#!/usr/bin/env python3
"""
12_compute_zone_density.py
===========================
Computes PNN and PV cell density (cells/mm²) per animal × brain_area × zone
by looking up individual cell coordinates in the MaskLoose/MaskStrict TIFFs
and dividing cell counts by the corresponding areaMm2 from cells_with_zones.csv.

This fills the gap between the cell-coordinate data (cells_with_zones.csv) and
the zone-level LMM pipeline: cells_with_zones.csv contains areaMm2 per zone
but not a zone-level cell count or density. This script provides those.

Zone assignment per cell coordinate:
    Core     = pixel is non-zero in MaskStrict
    Penumbra = pixel is non-zero in MaskLoose but zero in MaskStrict
    Outside  = pixel is zero in MaskLoose (or no mask exists for that slice)

Coordinate scaling:
    Cell coordinates (x_hires, y_hires) are in full hiRes image space.
    Masks are at reduced resolution (MASK_WIDTH_PX wide).
    Scale is derived per slice from the ratio of hiRes image dimensions to
    MASK_WIDTH_PX (masks are always this width from script 04).

INPUT:
    /path/to/cells_with_zones.csv
    /path/to/analysis_results/MaskLoose_C3-{animal}_{slice}.tif
    /path/to/analysis_results/MaskStrict_C3-{animal}_{slice}.tif

OUTPUT:
    /path/to/results/zone_density/zone_density.csv
    Columns: animal_id, treatment, brain_area, zone, staining,
             cell_count, areaMm2, density_mm2

The output joins into the R pipeline on (animal_id, brain_area, zone, staining).

Usage:
    python 12_compute_zone_density.py

Requirements:
    pip install numpy pandas tifffile pillow
"""

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from tifffile import imread

warnings.filterwarnings("ignore")

# =============================================================================
# CONFIGURATION
# =============================================================================
# cells_with_zones.csv produced by script 11
CELLS_CSV = Path("/path/to/cells_with_zones.csv")

# Directory containing MaskLoose/MaskStrict TIFFs (from script 04)
MASK_DIR = Path("/path/to/analysis_results")

# Output directory
OUTPUT_DIR = Path("/path/to/results/zone_density")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Label corrections ──────────────────────────────────────────────────────────
# Apply swap and rename corrections to match the statistical pipeline.
SWAP: dict = {}               # e.g. {"mScarlet_4": "ADAMTS4_MD"}
TREATMENT_RENAME: dict = {}   # e.g. {"C6ST1_ADAMTS4": "C6ST1_ADAMTS15"}
EXCLUDE_ANIMALS: list = []    # e.g. ["C6ST1_ADAMTS15_4"]

# Hemisphere to use (injection side)
IPSI_HEMI = "left"

# Known mask width in pixels — masks are always generated at this width
MASK_WIDTH_PX = 1500

# =============================================================================
# HELPERS
# =============================================================================

TREATMENT_PREFIXES = [
    "C6ST1_ADAMTS15", "C6ST1_ADAMTS4",
    "ADAMTS4_MD", "ADAMTS15", "ADAMTS4",
    "C6ST1", "mScarlet",
]


def get_treatment(animal_id: str) -> str:
    for p in TREATMENT_PREFIXES:
        if str(animal_id).startswith(p):
            return p
    return "Unknown"


def apply_swap(animal_id: str, treatment: str) -> str:
    num = str(animal_id).split("_")[-1]
    key = f"{treatment}_{num}"
    return SWAP.get(key, treatment)


def load_mask(mask_dir: Path, animal_id: str, slice_id: str, mask_type: str):
    """Load a MaskLoose or MaskStrict TIFF as a boolean array, or None."""
    fname = mask_dir / f"{mask_type}_C3-{animal_id}_{slice_id}.tif"
    if not fname.exists():
        return None
    try:
        return imread(str(fname)) > 0
    except Exception as exc:
        warnings.warn(f"Cannot load {fname.name}: {exc}")
        return None


def assign_zones(
    x_hires: np.ndarray,
    y_hires: np.ndarray,
    strict_mask,
    loose_mask,
    scale_x: float,
    scale_y: float,
) -> np.ndarray:
    """
    Return zone labels ('Core', 'Penumbra', 'Outside') for each cell.

    Cell hiRes coordinates are scaled to mask resolution and the mask pixel
    at that position determines the zone.
    """
    n     = len(x_hires)
    zones = np.full(n, "Outside", dtype=object)

    if strict_mask is None and loose_mask is None:
        return zones

    ref  = loose_mask if loose_mask is not None else strict_mask
    h, w = ref.shape

    mx = np.clip((x_hires * scale_x).astype(int), 0, w - 1)
    my = np.clip((y_hires * scale_y).astype(int), 0, h - 1)

    if loose_mask is not None:
        zones[loose_mask[my, mx]]  = "Penumbra"
    if strict_mask is not None:
        zones[strict_mask[my, mx]] = "Core"

    return zones

# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    print("=" * 60)
    print("COMPUTE ZONE DENSITY")
    print("=" * 60)

    # ── Load cells_with_zones ─────────────────────────────────────────────────
    print(f"\nLoading {CELLS_CSV.name} ...")
    czv = pd.read_csv(CELLS_CSV)
    print(f"  {len(czv):,} rows, {czv['mouse_id'].nunique()} animals")

    # Apply corrections
    if EXCLUDE_ANIMALS:
        czv = czv[~czv["mouse_id"].isin(EXCLUDE_ANIMALS)]
    if SWAP:
        czv["treatment"] = czv.apply(
            lambda r: apply_swap(r["mouse_id"], get_treatment(r["mouse_id"])), axis=1
        )
    if TREATMENT_RENAME:
        czv["treatment"] = czv["treatment"].replace(TREATMENT_RENAME)

    # Restrict to ipsilateral hemisphere only
    czv = czv[czv["hemisphere"] == IPSI_HEMI]

    # ── Derive areaMm2 per animal × brain_area × zone ────────────────────────
    # Average across slices — areaMm2 in cells_with_zones is per-slice
    area_lut = (
        czv.groupby(["mouse_id", "brain_area", "zone"], as_index=False)["areaMm2"]
        .mean()
        .rename(columns={"mouse_id": "animal_id"})
    )

    # ── Infer hiRes image dimensions for scale factor ─────────────────────────
    # Use the maximum x/y coordinate per animal as a proxy for image width/height
    dims = (
        czv.groupby("mouse_id")
        .agg(img_w=("x_hires", "max"), img_h=("y_hires", "max"))
        .reset_index()
        .rename(columns={"mouse_id": "animal_id"})
    )
    # Add 10% margin to account for cells not reaching the image edge
    dims["img_w"] = (dims["img_w"] * 1.10).astype(int)
    dims["img_h"] = (dims["img_h"] * 1.10).astype(int)

    # ── Assign zones from mask lookup ─────────────────────────────────────────
    print("\nAssigning zones from mask lookup ...")
    all_records = []

    for animal_id, animal_df in czv.groupby("mouse_id"):
        treatment = get_treatment(animal_id)
        if SWAP:
            treatment = apply_swap(animal_id, treatment)
        if TREATMENT_RENAME:
            treatment = TREATMENT_RENAME.get(treatment, treatment)

        # Get image dimensions for this animal
        dim_row = dims[dims["animal_id"] == animal_id]
        img_w   = int(dim_row["img_w"].values[0]) if len(dim_row) else 18000
        img_h   = int(dim_row["img_h"].values[0]) if len(dim_row) else 13000

        for slice_id, slice_df in animal_df.groupby("slice_id"):
            loose  = load_mask(MASK_DIR, animal_id, slice_id, "MaskLoose")
            strict = load_mask(MASK_DIR, animal_id, slice_id, "MaskStrict")

            if loose is not None:
                h_mask, w_mask = loose.shape
            elif strict is not None:
                h_mask, w_mask = strict.shape
            else:
                h_mask, w_mask = MASK_WIDTH_PX * img_h // img_w, MASK_WIDTH_PX

            scale_x = w_mask / img_w
            scale_y = h_mask / img_h

            x = slice_df["x_hires"].values.astype(float)
            y = slice_df["y_hires"].values.astype(float)

            zones    = assign_zones(x, y, strict, loose, scale_x, scale_y)
            staining = slice_df["cell_type"].values  # PNN or PV

            for zone_label in ["Core", "Penumbra", "Outside"]:
                zone_mask = zones == zone_label
                for stain in ["PNN", "PV"]:
                    stain_mask = (staining == stain) & zone_mask
                    if not stain_mask.any():
                        continue

                    cells_in_zone = slice_df[stain_mask]
                    for brain_area, area_grp in cells_in_zone.groupby("brain_area"):
                        all_records.append({
                            "animal_id":  animal_id,
                            "treatment":  treatment,
                            "slice_id":   slice_id,
                            "brain_area": brain_area,
                            "zone":       zone_label,
                            "staining":   stain,
                            "cell_count": len(area_grp),
                        })

        print(f"  ✓ {animal_id}")

    if not all_records:
        print("No records generated.")
        return

    # ── Aggregate to animal × area × zone (sum across slices) ────────────────
    records_df = pd.DataFrame(all_records)
    density_df = (
        records_df
        .groupby(["animal_id", "treatment", "brain_area", "zone", "staining"],
                 as_index=False)["cell_count"]
        .sum()
    )

    # ── Join areaMm2 ──────────────────────────────────────────────────────────
    density_df = density_df.merge(
        area_lut, on=["animal_id", "brain_area", "zone"], how="left"
    )
    density_df["density_mm2"] = (
        density_df["cell_count"] / density_df["areaMm2"]
    ).replace([np.inf, -np.inf], np.nan)

    # ── Save ──────────────────────────────────────────────────────────────────
    out_path = OUTPUT_DIR / "zone_density.csv"
    density_df.to_csv(out_path, index=False)

    print(f"\nSaved: {out_path}")
    print(f"  Rows:    {len(density_df)}")
    print(f"  Animals: {density_df['animal_id'].nunique()}")
    print(f"  Zones:   {density_df['zone'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
