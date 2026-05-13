#!/usr/bin/env python3
"""
11_assign_zones.py
==================
Assigns each detected cell (from _cell_fluorescence_analysis.csv files) to
an injection zone (Core / Penumbra / Outside) by looking up the corresponding
pixel in the MaskLoose and MaskStrict TIFFs.

Zone logic (per pixel at scaled coordinates):
    Core     = pixel is non-zero in MaskStrict
    Penumbra = pixel is non-zero in MaskLoose but zero in MaskStrict
    Outside  = pixel is zero in MaskLoose (or no mask exists for this slice)

Cell coordinates are in hiRes image space. The mask TIFFs are at reduced
resolution (~1005 × 1500 px); coordinates are scaled linearly before lookup.
The scale factor is derived from the ratio of hiRes image dimensions to the
known mask width (MASK_WIDTH_PX).

Label corrections (swap, rename, exclusion) are applied before saving.

INPUT:
    /path/to/analysis_results/*_cell_fluorescence_analysis.csv
    /path/to/analysis_results/MaskLoose_C3-{animal}_{slice}.tif
    /path/to/analysis_results/MaskStrict_C3-{animal}_{slice}.tif
    /path/to/DATA/                 (for hiRes image dimensions)

OUTPUT:
    /path/to/analysis_results/cells_with_zones.csv

Usage:
    python 11_assign_zones.py

Requirements:
    pip install numpy pandas tifffile tqdm
"""

import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from tifffile import imread
from tqdm import tqdm

warnings.filterwarnings("ignore")

# =============================================================================
# CONFIGURATION
# =============================================================================
# Directory containing *_cell_fluorescence_analysis.csv and Mask*.tif files
CELL_RESULTS_PATH = Path("/path/to/analysis_results")
INJECTION_MASK_DIR = Path("/path/to/analysis_results")
OUTPUT_PATH        = Path("/path/to/analysis_results")

# DATA directory (used to locate hiRes images for dimension inference)
DATA_PATH = Path("/path/to/DATA")

# ── Label corrections ──────────────────────────────────────────────────────────
# Confirmed label swaps: map animal_id → corrected treatment string.
# Leave empty ({}) if no swaps apply to your experiment.
SWAP_ANIMALS: dict = {}   # e.g. {"mScarlet_4": "ADAMTS4_MD", "ADAMTS4_MD_4": "mScarlet"}

# Confirmed injection failures: these animals are excluded from the output.
EXCLUDE_ANIMALS: list = []   # e.g. ["C6ST1_ADAMTS15_4"]

# Treatment name rename (old → new). Applied after swap correction.
TREATMENT_RENAME: dict = {}   # e.g. {"C6ST1_ADAMTS4": "C6ST1_ADAMTS15"}

# Apply the label swap to both the animal_id and treatment columns
APPLY_SWAP_CORRECTION = True

# Known mask resolution width (px). MaskLoose/MaskStrict TIFFs are always
# generated at this width by 04_generate_injection_masks.py.
MASK_WIDTH_PX = 1500

# =============================================================================
# HELPERS
# =============================================================================

def get_treatment(animal_id: str) -> str:
    """Infer treatment group from animal ID by longest-matching prefix."""
    PREFIXES = [
        "C6ST1_ADAMTS15", "C6ST1_ADAMTS4",
        "ADAMTS4_MD", "ADAMTS15", "ADAMTS4",
        "C6ST1", "mScarlet",
    ]
    for p in PREFIXES:
        if str(animal_id).startswith(p):
            return p
    return "Unknown"


def apply_swap(animal_id: str, treatment: str) -> str:
    """Return the corrected treatment label for an animal, if a swap applies."""
    num = str(animal_id).split("_")[-1]
    key = f"{treatment}_{num}"
    return SWAP_ANIMALS.get(key, treatment)


def load_mask(mask_dir: Path, animal_id: str, slice_id: str, mask_type: str) -> np.ndarray | None:
    """Load a MaskLoose or MaskStrict TIFF as a boolean array."""
    fname = mask_dir / f"{mask_type}_C3-{animal_id}_{slice_id}.tif"
    if not fname.exists():
        return None
    try:
        img = imread(str(fname))
        return img > 0
    except Exception as exc:
        print(f"  Warning: cannot load {fname.name}: {exc}")
        return None


def assign_zone(
    x_hires: np.ndarray,
    y_hires: np.ndarray,
    strict_mask: np.ndarray | None,
    loose_mask:  np.ndarray | None,
    scale_x: float,
    scale_y: float,
) -> np.ndarray:
    """
    Assign injection zones to all cells in a slice.

    Cell hiRes coordinates are scaled to mask resolution, then the mask pixel
    at that position determines the zone.

    Returns a string array with values 'Core', 'Penumbra', or 'Outside'.
    """
    n = len(x_hires)
    zones = np.full(n, "Outside", dtype=object)

    if strict_mask is None and loose_mask is None:
        return zones

    ref_mask = loose_mask if loose_mask is not None else strict_mask
    h, w     = ref_mask.shape

    mx = np.clip((x_hires * scale_x).astype(int), 0, w - 1)
    my = np.clip((y_hires * scale_y).astype(int), 0, h - 1)

    if loose_mask is not None:
        in_loose           = loose_mask[my, mx]
        zones[in_loose]    = "Penumbra"

    if strict_mask is not None:
        in_strict          = strict_mask[my, mx]
        zones[in_strict]   = "Core"

    return zones

# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    print("=" * 60)
    print("ZONE ASSIGNMENT → cells_with_zones.csv")
    print("=" * 60)

    cell_files = sorted(CELL_RESULTS_PATH.glob("*_cell_fluorescence_analysis.csv"))
    cell_files = [f for f in cell_files if not f.name.startswith("TEST_")]
    print(f"Found {len(cell_files)} animal file(s).\n")

    all_records = []

    for cell_file in tqdm(cell_files, desc="Animals"):
        animal_id = cell_file.stem.split("_cell_fluorescence")[0]

        if animal_id in EXCLUDE_ANIMALS:
            print(f"  EXCLUDED: {animal_id}")
            continue

        df = pd.read_csv(cell_file)
        if df.empty:
            continue

        # Apply swap correction
        if APPLY_SWAP_CORRECTION and SWAP_ANIMALS:
            raw_treatment = df["cell_type"].apply(lambda _: get_treatment(animal_id))
            df["treatment"] = raw_treatment.apply(lambda t: apply_swap(animal_id, t))
        else:
            df["treatment"] = get_treatment(animal_id)

        # Apply rename
        if TREATMENT_RENAME:
            df["treatment"] = df["treatment"].replace(TREATMENT_RENAME)

        # Infer hiRes image dimensions for scaling (from first available C1 file)
        hires_dir = DATA_PATH / animal_id / "hiRes"
        img_dims  = None
        if hires_dir.exists():
            for f in sorted(hires_dir.iterdir()):
                if f.name.endswith("-C1.tif"):
                    try:
                        img = imread(str(f))
                        img_dims = img.shape[:2]  # (H, W)
                        break
                    except Exception:
                        pass

        # Process each slice
        slice_records = []
        for slice_id, slice_df in df.groupby("slice_id"):
            loose  = load_mask(INJECTION_MASK_DIR, animal_id, slice_id, "MaskLoose")
            strict = load_mask(INJECTION_MASK_DIR, animal_id, slice_id, "MaskStrict")

            # Compute scale factors
            if img_dims and (loose is not None or strict is not None):
                ref   = loose if loose is not None else strict
                h, w  = ref.shape
                # Use actual image dims if available; fallback to mask dims
                img_h, img_w = img_dims
                scale_x = w / img_w
                scale_y = h / img_h
            else:
                scale_x = scale_y = 1.0

            x_arr = slice_df["x_hires"].values.astype(float)
            y_arr = slice_df["y_hires"].values.astype(float)

            zones = assign_zone(x_arr, y_arr, strict, loose, scale_x, scale_y)

            slice_df = slice_df.copy()
            slice_df["zone"]      = zones
            slice_df["mouse_id"]  = animal_id
            slice_records.append(slice_df)

        if slice_records:
            animal_df = pd.concat(slice_records, ignore_index=True)
            all_records.append(animal_df)
            print(f"  ✓ {animal_id}: {len(animal_df)} cells")

    if not all_records:
        print("No records to save.")
        return

    output_df = pd.concat(all_records, ignore_index=True)

    # Standardise column order
    priority_cols = ["mouse_id", "slice_id", "cell_type", "hemisphere", "zone",
                     "x_hires", "y_hires", "treatment", "regionID", "brain_area",
                     "absolute", "mean", "normalized_btm20", "normalized_median",
                     "top4_mean", "top4_median"]
    col_order = [c for c in priority_cols if c in output_df.columns]
    col_order += [c for c in output_df.columns if c not in col_order]
    output_df = output_df[col_order]

    out_path = OUTPUT_PATH / "cells_with_zones.csv"
    output_df.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}")
    print(f"  Total cells: {len(output_df)}")
    print(f"  Animals:     {output_df['mouse_id'].nunique()}")
    if "zone" in output_df.columns:
        print(f"  Zone counts:\n{output_df['zone'].value_counts().to_string()}")


if __name__ == "__main__":
    main()
