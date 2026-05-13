#!/usr/bin/env python3
"""
09_cell_fluorescence_analysis.py
=================================
Per-cell fluorescence intensity analysis for PNN and PV cells.

For each animal × slice:
  1. Load C1 (PNN) and C2 (PV) hiRes TIFFs, the RGB atlas overlay PNG,
     and tissue masks (regular and, where present, manually edited).
  2. Load cell coordinates from DATA/{animal}/counts/ (output of 07_csvsplitter.py).
  3. Assign each cell to a brain region via the region label map.
  4. Assign each cell to a hemisphere (left/right) from the QuickNII anchoring.
  5. For each valid cell (centroid inside tissue mask), extract a
     PATCH_SIZE × PATCH_SIZE pixel patch and a larger context window, then
     compute per-cell metrics:
       absolute         — sum of patch pixel values
       mean             — mean patch intensity
       normalized_btm20 — (96th percentile of patch − 20th percentile of context)
                          / 20th percentile of context
       normalized_median — same but normalised to median of context
       top4_mean/median  — mean/median of the top 4% of patch pixels

  6. Save per-animal results to OUTPUT_PATH/{animal_id}_cell_fluorescence_analysis.csv.

Edited masks can be used for specific animals/slices by populating the
SPECIAL_CASES dict. Most experiments can leave it empty.

INPUT:
    DATA/{animal}/hiRes/       — C1, C2 TIFFs
    DATA/{animal}/RGBatlas/    — atlas PNG overlays
    DATA/{animal}/masks/       — tissue masks
    DATA/{animal}/counts/      — per-slice coordinate CSVs (from 07)
    /path/to/Visualign_json/   — VisuAlign JSON files

OUTPUT:
    /path/to/analysis_results/{animal_id}_cell_fluorescence_analysis.csv

Usage:
    python 09_cell_fluorescence_analysis.py

Requirements:
    pip install numpy pandas scikit-image scipy tqdm
"""

import os
import re
import json
import time
from typing import List, Tuple, Dict, Optional, Any

import numpy as np
import pandas as pd
from skimage import io, transform
import scipy.ndimage as ndi
from tqdm.auto import tqdm

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_PATH   = "/path/to/project"
DATA_ROOT   = os.path.join(BASE_PATH, "DATA")
OUTPUT_PATH = os.path.join(BASE_PATH, "analysis_results")
JSON_LOOKUP_PATH = os.path.join(BASE_PATH, "rainbow_2017.json")

# VisuAlign JSON location pattern. {animal_id} is substituted at runtime.
VISUALIGN_JSON_PATTERN = "/path/to/Visualign_json/{animal_id}_visualign.json"

# Allen CCF midline (voxels)
ATLAS_MIDLINE_VOXELS = 228

# ── Analysis parameters ────────────────────────────────────────────────────────
PATCH_SIZE          = 45     # pixels, centred on each cell
CONTEXT_MULTIPLIER  = 6      # context window = PATCH_SIZE × CONTEXT_MULTIPLIER
TOP_PERCENTILE      = 96     # top 4% threshold
BOTTOM_PERCENTILE   = 20     # normalisation baseline percentile
MIN_CONTEXT_PIXELS  = 100    # minimum mask pixels needed in context for normalisation

# ── Image calibration ──────────────────────────────────────────────────────────
MICRONS_PER_PIXEL_HIRES = 0.645

# ── Output precision ──────────────────────────────────────────────────────────
FLOAT_PRECISION = 6

# ── Edited mask overrides ─────────────────────────────────────────────────────
# By default, PNN cells use the regular tissue mask and PV cells use the
# edited mask (if available). If a specific animal/slice pair should use
# the edited mask for BOTH channels, add it here:
#   { "animal_id": {"s001", "s003"} }  — specific slices
#   { "animal_id": "ALL" }             — all slices for that animal
# Leave empty for most experiments.
SPECIAL_CASES: Dict[str, Any] = {}

# ── Test mode ──────────────────────────────────────────────────────────────────
TEST_MODE            = False
TEST_ANIMAL_IDS      = []
TEST_SPECIFIC_SLICES = None
TEST_SLICES          = 2

# =============================================================================
# QUICKNII HELPERS  (shared with 08_diffu_spot_quantification.py)
# =============================================================================

def load_quicknii_json(json_path: str) -> Optional[Dict[str, Dict]]:
    if not os.path.exists(json_path):
        return None
    with open(json_path) as f:
        data = json.load(f)
    result = {}
    for s in data.get("slices", []):
        fname, anchoring = s.get("filename", ""), s.get("anchoring", None)
        width, height    = s.get("width"), s.get("height")
        if fname and anchoring and len(anchoring) == 9:
            m = re.search(r"(s\d+)", fname)
            if m:
                result[m.group(1)] = {
                    "anchoring": anchoring, "width": width, "height": height
                }
    return result or None


def find_quicknii_json(animal_path: str, mouse_id: str) -> Optional[str]:
    p = VISUALIGN_JSON_PATTERN.format(animal_path=animal_path, animal_id=mouse_id)
    return p if os.path.exists(p) else None


def create_hemisphere_masks(
    slice_info: Dict, mask_shape: Tuple, hires_shape: Tuple
) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """Identical logic to 08_diffu_spot_quantification.py — see that file."""
    ox, oy, oz, ux, uy, uz, vx, vy, vz = slice_info["anchoring"]
    json_w = slice_info.get("width",  mask_shape[1])
    json_h = slice_info.get("height", mask_shape[0])

    mask_h, mask_w = mask_shape
    y_c, x_c       = np.mgrid[0:mask_h, 0:mask_w]
    atlas_x         = ox + (x_c / json_w) * ux + (y_c / json_h) * vx

    l_res = (atlas_x <  ATLAS_MIDLINE_VOXELS).astype(np.uint8)
    r_res = (atlas_x >= ATLAS_MIDLINE_VOXELS).astype(np.uint8)

    hires_h, hires_w = hires_shape
    hc = np.indices(hires_shape)
    sc = np.array([hc[0] * (mask_h / hires_h), hc[1] * (mask_w / hires_w)])

    l_hires = ndi.map_coordinates(l_res, sc, order=0).astype(bool)
    r_hires = ndi.map_coordinates(r_res, sc, order=0).astype(bool)

    y_rows    = np.arange(mask_h) / json_h
    mx        = ((ATLAS_MIDLINE_VOXELS - ox - y_rows * vx) * json_w / ux
                 if abs(ux) > 1e-10 else np.full(mask_h, np.nan))
    lp, rp    = l_res.sum(), r_res.sum()
    total_pix = lp + rp

    midline_info = {
        "method":          "quicknii_anchoring",
        "mean_x":          float(np.nanmean(mx)),
        "angle_deg":       float(np.degrees(np.arctan2(
            mx[-1] - mx[0] if not np.isnan(mx[0]) else 0, len(mx)))),
        "left_pct":        float(100 * lp / total_pix) if total_pix > 0 else 0,
        "right_pct":       float(100 * rp / total_pix) if total_pix > 0 else 0,
        "anchoring":       slice_info["anchoring"],
        "json_dims":       [json_h, json_w],
        "mask_dims":       list(mask_shape),
        "ux":              float(ux),
    }
    return l_hires, r_hires, midline_info


def get_hemisphere(coords: np.ndarray, left: np.ndarray, right: np.ndarray) -> List[str]:
    """Return hemisphere label for each (x, y) coordinate pair."""
    h, w      = left.shape
    result    = []
    for x, y in coords:
        xi, yi = int(round(x)), int(round(y))
        if 0 <= xi < w and 0 <= yi < h:
            result.append("left" if left[yi, xi] else "right" if right[yi, xi] else "unknown")
        else:
            result.append("unknown")
    return result

# =============================================================================
# FILE UTILITIES
# =============================================================================

def find_file(directory: str, mouse_id: str, slice_id: str, suffix: str) -> Optional[str]:
    if not os.path.exists(directory):
        return None
    for fn in os.listdir(directory):
        if mouse_id in fn and slice_id in fn and fn.endswith(suffix):
            return os.path.join(directory, fn)
    return None


def find_mask(mask_dir: str, mouse_id: str, slice_id: str, prefer_edited: bool = False) -> Optional[str]:
    if not os.path.exists(mask_dir):
        return None
    matches = [f for f in os.listdir(mask_dir) if f.endswith(".png")
               and mouse_id in f and slice_id in f]
    if not matches:
        return None
    if len(matches) == 1:
        return os.path.join(mask_dir, matches[0])
    key    = "edited" if prefer_edited else ""
    ranked = [f for f in matches if key in f.lower()] or matches
    return os.path.join(mask_dir, sorted(ranked, key=len)[0])


def load_coords(counts_dir: str, mouse_id: str, slice_id: str, cell_type: str) -> np.ndarray:
    """Load (x, y) coordinates from a per-slice counts CSV."""
    if not os.path.exists(counts_dir):
        return np.array([])
    for fn in os.listdir(counts_dir):
        if not fn.endswith(".csv"):
            continue
        if mouse_id in fn and slice_id in fn:
            if (f"cells_{cell_type}" in fn or f"-{cell_type}.csv" in fn
                    or f"_{cell_type}.csv" in fn):
                df = pd.read_csv(os.path.join(counts_dir, fn))
                if "x" in df.columns and "y" in df.columns:
                    return df[["x", "y"]].values
                elif "X" in df.columns and "Y" in df.columns:
                    return df[["X", "Y"]].values
                elif df.shape[1] >= 2:
                    return df.iloc[:, :2].values
    return np.array([])


def build_lut(region_lookup: list) -> Tuple[np.ndarray, Dict, Dict]:
    max_val    = 256 ** 3
    lut        = np.zeros(max_val, dtype=np.int32)
    id_to_name = {}
    id_to_rid  = {}
    counter    = 1
    for reg in region_lookup:
        if reg.get("name") == "empty":
            continue
        ch = (reg["red"] << 16) | (reg["green"] << 8) | reg["blue"]
        lut[ch]           = counter
        id_to_name[counter] = reg["name"]
        id_to_rid[counter]  = reg["index"]
        counter += 1
    return lut, id_to_name, id_to_rid

# =============================================================================
# CELL METRICS
# =============================================================================

def extract_patch(img: np.ndarray, mask: np.ndarray, x: int, y: int, half: int):
    h, w = img.shape
    return (img [max(0,y-half):min(h,y+half), max(0,x-half):min(w,x+half)],
            mask[max(0,y-half):min(h,y+half), max(0,x-half):min(w,x+half)])


def cell_metrics(patch: np.ndarray, context: np.ndarray, mask_p: np.ndarray,
                 mask_c: np.ndarray) -> tuple:
    valid_p = patch[mask_p]
    valid_c = context[mask_c]
    if valid_p.size == 0:
        return 0.0, 0.0, np.nan, np.nan, np.nan, np.nan

    absolute   = float(valid_p.sum())
    mean_int   = float(valid_p.mean())
    t96        = float(np.percentile(valid_p, TOP_PERCENTILE))
    top4_px    = valid_p[valid_p >= t96]
    top4_mean  = float(top4_px.mean())   if top4_px.size > 0 else np.nan
    top4_med   = float(np.median(top4_px)) if top4_px.size > 0 else np.nan

    if valid_c.size >= MIN_CONTEXT_PIXELS:
        b20      = float(np.percentile(valid_c, BOTTOM_PERCENTILE))
        med      = float(np.median(valid_c))
        n_b20    = (t96 - b20) / abs(b20) if abs(b20) > 0.01 else t96 - b20
        n_med    = (t96 - med) / abs(med) if abs(med) > 0.01 else t96 - med
    else:
        n_b20 = n_med = np.nan

    return absolute, mean_int, n_b20, n_med, top4_mean, top4_med


def process_cells(
    coords: np.ndarray,
    cell_type: str,
    img: np.ndarray,
    mask: np.ndarray,
    label_map: np.ndarray,
    mouse_id: str,
    slice_id: str,
    left_hemi: Optional[np.ndarray],
    right_hemi: Optional[np.ndarray],
    midline_info: Optional[Dict],
) -> List[Dict]:
    """Compute per-cell metrics for all coordinates of one cell type."""
    results   = []
    half_p    = PATCH_SIZE // 2
    half_c    = (PATCH_SIZE * CONTEXT_MULTIPLIER) // 2
    coords_i  = np.round(coords).astype(int)

    valid_idx = [i for i, (x, y) in enumerate(coords_i)
                 if 0 <= x < mask.shape[1] and 0 <= y < mask.shape[0] and mask[y, x]]

    if not valid_idx:
        return results

    valid_coords = coords_i[valid_idx]
    hemispheres  = (get_hemisphere(valid_coords, left_hemi, right_hemi)
                    if left_hemi is not None else ["unknown"] * len(valid_idx))

    for i, idx in enumerate(valid_idx):
        x, y       = coords_i[idx]
        p,  mp     = extract_patch(img, mask, x, y, half_p)
        c,  mc     = extract_patch(img, mask, x, y, half_c)
        abs_, mn, nb20, nmed, t4m, t4md = cell_metrics(p, c, mp, mc)

        rec = {
            "mouse_id":          mouse_id,
            "slice_id":          slice_id,
            "cell_type":         cell_type,
            "hemisphere":        hemispheres[i],
            "x_hires":           x,
            "y_hires":           y,
            "label_id":          label_map[y, x],
            "absolute":          abs_,
            "mean":              mn,
            "normalized_btm20":  nb20,
            "normalized_median": nmed,
            "top4_mean":         t4m,
            "top4_median":       t4md,
        }
        if midline_info:
            rec["midline_method"]    = midline_info["method"]
            rec["midline_angle_deg"] = midline_info["angle_deg"]
        results.append(rec)

    return results

# =============================================================================
# SLICE PROCESSING
# =============================================================================

def use_edited_for_both(mouse_id: str, slice_id: str) -> bool:
    sc = SPECIAL_CASES.get(mouse_id)
    if sc is None:
        return False
    return sc == "ALL" or slice_id in sc


def process_slice(
    mouse_id: str,
    slice_id: str,
    animal_path: str,
    lut: np.ndarray,
    has_edited: bool,
    anchoring: Optional[Dict],
) -> Tuple[List[Dict], Optional[Dict]]:

    hires_dir  = os.path.join(animal_path, "hiRes")
    counts_dir = os.path.join(animal_path, "counts")
    masks_dir  = os.path.join(animal_path, "masks")
    edited_dir = os.path.join(animal_path, "edited_masks")
    atlas_dir  = os.path.join(animal_path, "RGBatlas")

    pnn_path = find_file(hires_dir, mouse_id, slice_id, "-C1.tif")
    pv_path  = find_file(hires_dir, mouse_id, slice_id, "-C2.tif")
    atl_path = (find_file(atlas_dir, mouse_id, slice_id, "_nl.png") or
                find_file(atlas_dir, mouse_id, slice_id, ".png"))

    if not all([pnn_path, pv_path, atl_path]):
        print(f"  Warning: missing images for {mouse_id} {slice_id}")
        return [], None

    img_pnn   = io.imread(pnn_path)
    img_pv    = io.imread(pv_path)
    rgb_atlas = io.imread(atl_path)[:, :, :3]
    hires_shape = img_pnn.shape

    # Mask selection
    force_edited = use_edited_for_both(mouse_id, slice_id)
    reg_mask_p   = find_mask(masks_dir,  mouse_id, slice_id)
    edit_mask_p  = find_mask(edited_dir, mouse_id, slice_id, prefer_edited=True) if has_edited else None

    if not reg_mask_p:
        print(f"  Warning: no mask for {mouse_id} {slice_id}")
        return [], None

    if force_edited and edit_mask_p:
        mask_pnn = mask_pv = io.imread(edit_mask_p)
    else:
        mask_pnn = io.imread(reg_mask_p)
        mask_pv  = io.imread(edit_mask_p) if edit_mask_p else mask_pnn

    # Scale atlas and masks to hiRes
    scaled_atlas  = transform.resize(rgb_atlas, hires_shape,
                                     order=0, preserve_range=True, anti_aliasing=False).astype(np.uint8)
    scaled_pnn_m  = transform.resize(mask_pnn, hires_shape,
                                     order=0, preserve_range=True, anti_aliasing=False) > 0
    scaled_pv_m   = transform.resize(mask_pv,  hires_shape,
                                     order=0, preserve_range=True, anti_aliasing=False) > 0

    rgb_int   = scaled_atlas.astype(np.int32)
    col_hash  = (rgb_int[:, :, 0] << 16) | (rgb_int[:, :, 1] << 8) | rgb_int[:, :, 2]
    label_map = lut[col_hash]

    # Hemisphere masks
    l_hemi = r_hemi = midline_info = None
    if anchoring and slice_id in anchoring:
        s_info        = anchoring[slice_id]
        raw_mask      = io.imread(reg_mask_p)
        mask_shape    = raw_mask.shape[:2] if raw_mask.ndim > 2 else raw_mask.shape
        l_hemi, r_hemi, midline_info = create_hemisphere_masks(
            s_info, mask_shape, hires_shape
        )

    pnn_coords = load_coords(counts_dir, mouse_id, slice_id, "C1")
    pv_coords  = load_coords(counts_dir, mouse_id, slice_id, "C2")

    results = []
    if pnn_coords.size:
        results.extend(process_cells(
            pnn_coords, "PNN", img_pnn, scaled_pnn_m,
            label_map, mouse_id, slice_id, l_hemi, r_hemi, midline_info
        ))
    if pv_coords.size:
        results.extend(process_cells(
            pv_coords, "PV", img_pv, scaled_pv_m,
            label_map, mouse_id, slice_id, l_hemi, r_hemi, midline_info
        ))

    print(f"    {slice_id}: {len(results)} cells processed")
    return results, midline_info

# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    start = time.time()
    print("=" * 60)
    print("CELL-SPECIFIC FLUORESCENCE ANALYSIS")
    if TEST_MODE:
        print(f"*** TEST MODE — {TEST_ANIMAL_IDS} ***")
    print("=" * 60)

    os.makedirs(OUTPUT_PATH, exist_ok=True)

    with open(JSON_LOOKUP_PATH) as f:
        region_lookup = json.load(f)
    lut, id_to_name, id_to_rid = build_lut(region_lookup)

    animals = sorted([
        d for d in os.listdir(DATA_ROOT)
        if os.path.isdir(os.path.join(DATA_ROOT, d))
    ])
    if TEST_MODE:
        animals = [a for a in animals if a in TEST_ANIMAL_IDS]

    print(f"Animals: {len(animals)}\n")

    for mouse_id in tqdm(animals, desc="Animals"):
        a_path    = os.path.join(DATA_ROOT, mouse_id)
        hires_dir = os.path.join(a_path, "hiRes")
        counts_dir = os.path.join(a_path, "counts")

        if not os.path.exists(counts_dir):
            print(f"  No counts folder for {mouse_id} — skipping")
            continue

        has_edited = os.path.exists(os.path.join(a_path, "edited_masks"))

        json_path = find_quicknii_json(a_path, mouse_id)
        anchoring = load_quicknii_json(json_path) if json_path else None
        if not anchoring:
            print(f"  No VisuAlign JSON for {mouse_id} — hemisphere will be 'unknown'")

        slice_ids = sorted({
            re.search(r"(s\d+)", f).group(1)
            for f in os.listdir(hires_dir)
            if re.search(r"(s\d+)", f) and f.endswith("-C1.tif")
        })

        if TEST_MODE:
            if TEST_SPECIFIC_SLICES:
                slice_ids = [s for s in slice_ids if s in TEST_SPECIFIC_SLICES]
            else:
                slice_ids = slice_ids[:TEST_SLICES]

        all_results  = []
        all_midlines = {}

        for sid in tqdm(slice_ids, desc=f"  {mouse_id}", leave=False):
            try:
                res, ml = process_slice(mouse_id, sid, a_path, lut, has_edited, anchoring)
                all_results.extend(res)
                if ml:
                    all_midlines[sid] = ml
            except Exception as exc:
                import traceback
                print(f"  Error {sid}: {exc}")
                traceback.print_exc()

        if not all_results:
            print(f"  No results for {mouse_id}")
            continue

        df = pd.DataFrame(all_results)
        df["brain_area"] = df["label_id"].map(id_to_name).fillna("Unknown")
        df["regionID"]   = df["label_id"].map(id_to_rid).fillna(0)

        final_cols = ["mouse_id", "slice_id", "cell_type", "hemisphere",
                      "x_hires", "y_hires", "regionID", "brain_area",
                      "absolute", "mean", "normalized_btm20", "normalized_median",
                      "top4_mean", "top4_median"]
        if "midline_method" in df.columns:
            final_cols += ["midline_method", "midline_angle_deg"]
        df = df[[c for c in final_cols if c in df.columns]]

        out_name = f"TEST_{mouse_id}_cell_fluorescence_analysis.csv" if TEST_MODE else \
                   f"{mouse_id}_cell_fluorescence_analysis.csv"
        df.to_csv(os.path.join(OUTPUT_PATH, out_name), index=False,
                  float_format=f"%.{FLOAT_PRECISION}f")
        print(f"  Saved {len(df)} cells → {out_name}")

    print(f"\nDone in {(time.time()-start)/60:.1f} min")


if __name__ == "__main__":
    main()
