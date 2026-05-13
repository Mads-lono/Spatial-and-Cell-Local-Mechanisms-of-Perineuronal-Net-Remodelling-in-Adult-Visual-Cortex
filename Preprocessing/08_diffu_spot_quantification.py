#!/usr/bin/env python3
"""
08_diffu_spot_quantification.py
================================
Region-level diffuse fluorescence quantification for PNN (WFA, C1) and PV
(C2) channels across all brain areas and injection zones.

For each animal × slice:
  1. Load C1 (PNN) and C2 (PV) hiRes TIFFs and the corresponding RGB atlas
     overlay PNG (produced by 02_flat_to_png.py from VisuAlign).
  2. Map atlas colours to region labels using the Rainbow 2017 JSON.
  3. Load tissue mask (std_mask from ilastik, 05_tissue_masker.py) and
     injection zone masks (MaskLoose/MaskStrict from 04_generate_injection_masks.py).
  4. Derive three injection zones:
       Core     = std_mask & MaskStrict
       Penumbra = std_mask & (MaskLoose & ~MaskStrict)
       Outside  = std_mask & ~MaskLoose
  5. Detect the hemisphere boundary from the QuickNII/VisuAlign anchoring
     vectors (atlas_x = ox + (x/w)*ux + (y/h)*vx; midline at 228 voxels).
  6. For each combination of (cell_type, zone, hemisphere) measure:
       areaPx, areaMm2, total_intensity, mean_intensity, mean_norm_btm20
     where mean_norm_btm20 is normalised to the 20th percentile of
     non-zero pixels in the Outside zone (same hemisphere, same stain).
  7. Save per-animal results to OUTPUT_PATH/{animal_id}_NESTED_RESULTS.csv
     and hemisphere midline QC to {animal_id}_MIDLINE_DATA.json.

INPUT:
    DATA/{animal}/hiRes/          — C1, C2 TIFFs (split by 00)
    DATA/{animal}/RGBatlas/       — atlas PNG overlays (from 02)
    DATA/{animal}/masks/          — tissue masks (from 05)
    DATA/{animal}/edited_masks/   — manually corrected masks (optional)
    /path/to/analysis_results/    — MaskLoose/MaskStrict TIFFs (from 04)
    /path/to/Visualign_json/      — per-animal VisuAlign JSON files

OUTPUT:
    /path/to/analysis_results/{animal_id}_NESTED_RESULTS.csv
    /path/to/analysis_results/{animal_id}_MIDLINE_DATA.json

Usage:
    python 08_diffu_spot_quantification.py

Requirements:
    pip install numpy pandas scikit-image scipy tqdm
"""

import os
import re
import json
import numpy as np
import pandas as pd
from skimage import io, transform
import scipy.ndimage as ndi
from tqdm.auto import tqdm
import time

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_PATH      = "/path/to/project"
DATA_ROOT      = os.path.join(BASE_PATH, "DATA")
OUTPUT_PATH    = os.path.join(BASE_PATH, "analysis_results")
JSON_LOOKUP_PATH = os.path.join(BASE_PATH, "rainbow_2017.json")
EXCLUSION_FILE = os.path.join(BASE_PATH, "do_not_use.csv")
INJECTION_MASK_DIR = "/path/to/analysis_results"

# VisuAlign JSON location pattern. {animal_id} is substituted at runtime.
VISUALIGN_JSON_PATTERN = "/path/to/Visualign_json/{animal_id}_visualign.json"

# Allen CCF midline (x-axis, voxels). The atlas is 456 voxels wide; midline = 228.
ATLAS_MIDLINE_VOXELS = 228

# Image calibration
MICRONS_PER_PIXEL_HIRES = 0.645
PIXEL_AREA_MM2 = (MICRONS_PER_PIXEL_HIRES / 1000) ** 2

# Normalisation: 20th percentile of outside-zone pixels (same hemisphere)
NORMALIZATION_PERCENTILE = 20

# Test mode — set TEST_MODE = False for a full run
TEST_MODE            = False
TEST_ANIMAL_IDS      = []      # e.g. ["ADAMTS4_1"]
TEST_SPECIFIC_SLICES = None    # e.g. ["s011"] or None (use TEST_SLICES count)
TEST_SLICES          = 3

# =============================================================================
# QUICKNII JSON PARSING
# =============================================================================

def load_quicknii_json(json_path: str) -> dict | None:
    """
    Load a QuickNII/VisuAlign JSON and return a dict mapping slice_id to
    {'anchoring': [ox,oy,oz,ux,uy,uz,vx,vy,vz], 'width': w, 'height': h}.
    """
    if not os.path.exists(json_path):
        return None
    with open(json_path) as f:
        data = json.load(f)

    result = {}
    for s in data.get("slices", []):
        fname     = s.get("filename", "")
        anchoring = s.get("anchoring", None)
        width     = s.get("width",    None)
        height    = s.get("height",   None)
        if fname and anchoring and len(anchoring) == 9:
            m = re.search(r"(s\d+)", fname)
            if m:
                result[m.group(1)] = {
                    "anchoring": anchoring,
                    "width": width,
                    "height": height,
                }
    return result or None


def find_quicknii_json(animal_path: str, mouse_id: str) -> str | None:
    """Resolve the VisuAlign JSON path for an animal using the configured pattern."""
    json_path = VISUALIGN_JSON_PATTERN.format(
        animal_path=animal_path, animal_id=mouse_id
    )
    if os.path.exists(json_path):
        return json_path
    print(f"  WARNING: JSON not found at {json_path}")
    return None


def create_hemisphere_masks(
    slice_info: dict, mask_shape: tuple, hires_shape: tuple
) -> tuple:
    """
    Create left/right hemisphere boolean masks at hiRes resolution using
    the QuickNII anchoring formula:
        atlas_x = ox + (x / w) * ux + (y / h) * vx
    Pixels with atlas_x < ATLAS_MIDLINE_VOXELS → left hemisphere.

    Returns (left_mask_hires, right_mask_hires, midline_info_dict).
    """
    ox, oy, oz, ux, uy, uz, vx, vy, vz = slice_info["anchoring"]
    json_w = slice_info.get("width",  mask_shape[1])
    json_h = slice_info.get("height", mask_shape[0])

    mask_h, mask_w = mask_shape
    y_coords, x_coords = np.mgrid[0:mask_h, 0:mask_w]
    atlas_x = ox + (x_coords / json_w) * ux + (y_coords / json_h) * vx

    left_mask_res  = (atlas_x <  ATLAS_MIDLINE_VOXELS).astype(np.uint8)
    right_mask_res = (atlas_x >= ATLAS_MIDLINE_VOXELS).astype(np.uint8)

    # Upsample to hiRes
    hires_h, hires_w = hires_shape
    scale_y = mask_h / hires_h
    scale_x = mask_w / hires_w
    hires_coords = np.indices(hires_shape)
    scaled = np.array([hires_coords[0] * scale_y, hires_coords[1] * scale_x])

    left_hires  = ndi.map_coordinates(left_mask_res,  scaled, order=0).astype(bool)
    right_hires = ndi.map_coordinates(right_mask_res, scaled, order=0).astype(bool)

    # Midline QC
    y_rows      = np.arange(mask_h) / json_h
    midline_x   = (
        (ATLAS_MIDLINE_VOXELS - ox - y_rows * vx) * json_w / ux
        if abs(ux) > 1e-10
        else np.full(mask_h, np.nan)
    )
    in_image     = (midline_x >= 0) & (midline_x < mask_w)
    left_pix     = left_mask_res.sum()
    right_pix    = right_mask_res.sum()
    total_pix    = left_pix + right_pix

    midline_info = {
        "method":          "quicknii_anchoring",
        "mean_x":          float(np.nanmean(midline_x)),
        "min_x":           float(np.nanmin(midline_x)),
        "max_x":           float(np.nanmax(midline_x)),
        "angle_deg":       float(np.degrees(np.arctan2(
            midline_x[-1] - midline_x[0] if not np.isnan(midline_x[0]) else 0,
            len(midline_x),
        ))),
        "valid_rows_pct":  float(100 * in_image.sum() / len(midline_x)),
        "left_pct":        float(100 * left_pix  / total_pix) if total_pix > 0 else 0,
        "right_pct":       float(100 * right_pix / total_pix) if total_pix > 0 else 0,
        "anchoring":       slice_info["anchoring"],
        "json_dims":       [json_h, json_w],
        "mask_dims":       list(mask_shape),
        "ux":              float(ux),
    }
    return left_hires, right_hires, midline_info

# =============================================================================
# UTILITIES
# =============================================================================

def build_lut_dict(region_lookup: list) -> tuple:
    """Build colour→label and label→name/regionID mappings from the JSON."""
    color_to_label, id_to_name, id_to_reg = {}, {}, {}
    for i, reg in enumerate(region_lookup, 1):
        if reg.get("name") == "empty":
            continue
        h = (int(reg["red"]) << 16) | (int(reg["green"]) << 8) | int(reg["blue"])
        color_to_label[h] = i
        id_to_name[i]     = reg["name"]
        id_to_reg[i]      = reg["index"]
    all_ids = np.array(list(id_to_name.keys()))
    return color_to_label, id_to_name, id_to_reg, all_ids


def find_file(directory: str, mouse_id: str, slice_id: str, suffix: str) -> str | None:
    """Find a file in directory whose name contains mouse_id, slice_id, and suffix."""
    if not os.path.exists(directory):
        return None
    for f in os.listdir(directory):
        if mouse_id in f and slice_id in f and f.endswith(suffix):
            return os.path.join(directory, f)
    return None

# =============================================================================
# SLICE PROCESSING
# =============================================================================

def process_slice(
    slice_id: str,
    mouse_id: str,
    animal_path: str,
    lut_dict: dict,
    all_label_ids: np.ndarray,
    quicknii_anchoring: dict,
) -> tuple:
    """
    Process one slice and return (results_DataFrame, midline_info_dict).
    Returns (None, None) on any failure.
    """
    try:
        hires_dir  = os.path.join(animal_path, "hiRes")
        atlas_dir  = os.path.join(animal_path, "RGBatlas")
        mask_dir   = os.path.join(animal_path, "masks")
        edited_dir = os.path.join(animal_path, "edited_masks")

        c1_path  = find_file(hires_dir, mouse_id, slice_id, "-C1.tif")
        c2_path  = find_file(hires_dir, mouse_id, slice_id, "-C2.tif")
        atl_path = find_file(atlas_dir, mouse_id, slice_id, "_nl.png") or \
                   find_file(atlas_dir, mouse_id, slice_id, ".png")

        if not c1_path or not atl_path:
            return None, None

        img_pnn = io.imread(c1_path)
        img_pv  = io.imread(c2_path) if c2_path else np.zeros_like(img_pnn)
        rgb_atl = io.imread(atl_path)[:, :, :3]

        h_shape = img_pnn.shape
        a_shape = rgb_atl.shape[:2]
        scale   = np.array(a_shape) / np.array(h_shape)
        coords  = np.indices(h_shape) * scale[:, np.newaxis, np.newaxis]

        def load_scaled_mask(p):
            if not p or not os.path.exists(p):
                return np.zeros(h_shape, dtype=bool)
            m = io.imread(p)
            if m.ndim > 2:
                m = m[:, :, 0]
            m_rs = transform.resize(
                m.astype(float), a_shape, order=0, preserve_range=True, anti_aliasing=False
            )
            return ndi.map_coordinates(m_rs, coords, order=0).astype(bool)

        std_p  = (find_file(edited_dir, mouse_id, slice_id, ".png") or
                  find_file(mask_dir,   mouse_id, slice_id, ".png"))
        std_m  = load_scaled_mask(std_p)

        loose_p  = os.path.join(INJECTION_MASK_DIR, f"MaskLoose_C3-{mouse_id}_{slice_id}.tif")
        strict_p = os.path.join(INJECTION_MASK_DIR, f"MaskStrict_C3-{mouse_id}_{slice_id}.tif")
        loose_m  = load_scaled_mask(loose_p)
        strict_m = load_scaled_mask(strict_p)

        core_m     = std_m & strict_m
        penumbra_m = std_m & (loose_m & ~strict_m)
        outside_m  = std_m & ~loose_m

        # Map atlas colours to region labels
        atl_h     = ((rgb_atl[:, :, 0].astype(np.int32) << 16) |
                     (rgb_atl[:, :, 1].astype(np.int32) <<  8) |
                      rgb_atl[:, :, 2].astype(np.int32))
        s_h_img   = ndi.map_coordinates(atl_h, coords, order=0).astype(np.int32)
        base_lab  = np.zeros(h_shape, dtype=np.int32)
        for c in np.unique(s_h_img):
            if c in lut_dict:
                base_lab[s_h_img == c] = lut_dict[c]

        # Hemisphere masks
        if slice_id not in quicknii_anchoring:
            print(f"  ERROR: no anchoring for {slice_id} — skipping")
            return None, None

        slice_info = quicknii_anchoring[slice_id]
        mask_path  = std_p
        if mask_path:
            raw_mask   = io.imread(mask_path)
            mask_shape = raw_mask.shape[:2] if raw_mask.ndim > 2 else raw_mask.shape
        else:
            mask_shape = (slice_info.get("height", a_shape[0]),
                          slice_info.get("width",  a_shape[1]))

        left_hemi, right_hemi, midline_info = create_hemisphere_masks(
            slice_info, mask_shape, h_shape
        )

        # Measurement function
        def measure(img, mask, lab, hemi_mask, hemi_name):
            combined  = mask & hemi_mask
            baseline  = outside_m & hemi_mask
            px_base   = img[baseline]
            norm      = (np.percentile(px_base[px_base > 0], NORMALIZATION_PERCENTILE)
                         if px_base.size > 0 else 1)
            l_map = base_lab.copy()
            l_map[~combined] = 0
            px = ndi.sum_labels(np.ones_like(l_map, dtype=np.uint64), l_map, index=all_label_ids)
            sm = ndi.sum_labels(img, l_map, index=all_label_ids)
            df = pd.DataFrame({
                "label_id":         all_label_ids,
                "areaPx":           px,
                "total_intensity":  sm,
                "cell_type":        lab,
                "slice_id":         slice_id,
                "hemisphere":       hemi_name,
                "midline_method":   midline_info["method"],
                "midline_angle_deg": midline_info["angle_deg"],
            })
            df["mean_intensity"]   = df["total_intensity"] / df["areaPx"].replace(0, np.nan)
            df["mean_norm_btm20"]  = df["mean_intensity"] / norm
            return df

        results = []
        for hemi_mask, hemi_name in [(left_hemi, "left"), (right_hemi, "right")]:
            results.extend([
                measure(img_pnn, std_m,     "PNN_Total",    hemi_mask, hemi_name),
                measure(img_pnn, core_m,    "PNN_Core",     hemi_mask, hemi_name),
                measure(img_pnn, penumbra_m,"PNN_Penumbra", hemi_mask, hemi_name),
                measure(img_pnn, outside_m, "PNN_Outside",  hemi_mask, hemi_name),
                measure(img_pv,  std_m,     "PV_Total",     hemi_mask, hemi_name),
                measure(img_pv,  core_m,    "PV_Core",      hemi_mask, hemi_name),
                measure(img_pv,  penumbra_m,"PV_Penumbra",  hemi_mask, hemi_name),
                measure(img_pv,  outside_m, "PV_Outside",   hemi_mask, hemi_name),
            ])

        return pd.concat(results).fillna(0), midline_info

    except Exception as exc:
        import traceback
        print(f"  Error on {slice_id}: {exc}")
        traceback.print_exc()
        return None, None

# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    print("=" * 60)
    print("DIFFU_SPOT QUANTIFICATION")
    print(f"Atlas midline: {ATLAS_MIDLINE_VOXELS} voxels")
    if TEST_MODE:
        print(f"*** TEST MODE — animals: {TEST_ANIMAL_IDS}, "
              f"slices: {TEST_SPECIFIC_SLICES or TEST_SLICES} ***")
    print("=" * 60)

    os.makedirs(OUTPUT_PATH, exist_ok=True)

    with open(JSON_LOOKUP_PATH) as f:
        lut_dict, id_name, id_reg, all_ids = build_lut_dict(json.load(f))

    exclude = []
    if os.path.exists(EXCLUSION_FILE):
        exclude = pd.read_csv(EXCLUSION_FILE).iloc[:, 0].astype(str).tolist()
        print(f"Exclusion file loaded — skipping {len(exclude)} animal(s).")

    all_folders = [d for d in os.listdir(DATA_ROOT)
                   if os.path.isdir(os.path.join(DATA_ROOT, d))]
    animals     = sorted([d for d in all_folders if d not in exclude])

    if TEST_MODE:
        animals = [a for a in animals if a in TEST_ANIMAL_IDS]

    print(f"Animals to process: {len(animals)}\n")

    for mouse in tqdm(animals, desc="Animals"):
        a_path  = os.path.join(DATA_ROOT, mouse)
        h_dir   = os.path.join(a_path, "hiRes")
        if not os.path.exists(h_dir):
            continue

        json_path = find_quicknii_json(a_path, mouse)
        if not json_path:
            print(f"  ERROR: no VisuAlign JSON for {mouse} — skipping")
            continue
        anchoring = load_quicknii_json(json_path)
        if not anchoring:
            print(f"  ERROR: cannot load JSON for {mouse} — skipping")
            continue
        print(f"  JSON loaded ({len(anchoring)} slices)")

        slice_ids = sorted({
            re.search(r"(s\d+)", f).group(1)
            for f in os.listdir(h_dir)
            if re.search(r"(s\d+)", f)
        })

        if TEST_MODE:
            if TEST_SPECIFIC_SLICES:
                slice_ids = [s for s in slice_ids if s in TEST_SPECIFIC_SLICES]
            else:
                slice_ids = slice_ids[:TEST_SLICES]

        results      = []
        midline_data = {}

        for sid in tqdm(slice_ids, desc=f"  {mouse}", leave=False):
            res, ml = process_slice(sid, mouse, a_path, lut_dict, all_ids, anchoring)
            if res is not None:
                results.append(res)
                midline_data[sid] = ml

        if not results:
            continue

        full_df = pd.concat(results)
        full_df["brain_area"] = full_df["label_id"].map(id_name)
        full_df["regionID"]   = full_df["label_id"].map(id_reg)
        full_df["areaMm2"]    = full_df["areaPx"] * PIXEL_AREA_MM2
        full_df["animal_id"]  = mouse
        full_df = full_df[full_df["areaPx"] > 0]
        full_df = full_df[full_df["brain_area"] != "root"]

        col_order = ["animal_id", "brain_area", "regionID", "label_id",
                     "hemisphere", "slice_id", "cell_type",
                     "areaPx", "areaMm2", "total_intensity",
                     "mean_intensity", "mean_norm_btm20",
                     "midline_method", "midline_angle_deg"]
        full_df = full_df[[c for c in col_order if c in full_df.columns]]

        out_csv  = os.path.join(OUTPUT_PATH, f"{mouse}_NESTED_RESULTS.csv")
        out_json = os.path.join(OUTPUT_PATH, f"{mouse}_MIDLINE_DATA.json")
        full_df.to_csv(out_csv, index=False)
        with open(out_json, "w") as f:
            json.dump(midline_data, f, indent=2)
        print(f"  Saved: {out_csv} ({len(full_df)} rows)")

    print("\n" + "=" * 60)
    print("COMPLETE")


if __name__ == "__main__":
    main()
