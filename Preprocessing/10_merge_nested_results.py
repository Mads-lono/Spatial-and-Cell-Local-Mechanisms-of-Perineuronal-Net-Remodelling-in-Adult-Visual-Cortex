#!/usr/bin/env python3
"""
10_merge_nested_results.py
===========================
Merges per-animal diffuse fluorescence (_NESTED_RESULTS.csv) and cell
fluorescence (_cell_fluorescence_analysis.csv) outputs into two combined
datasets:

  merged_dataset_[timestamp].csv       — Total-zone hemisphere-level data
  merged_dataset_zones_[timestamp].csv — Zone-level diffuse data (Core/Penumbra/Outside)

The merge is performed on (brain_area, regionID, animal_id, staining,
hemisphere). Brain region metadata (acronym, sphinx_id) is joined from the
Allen Brain Atlas structures JSON when available.

INPUT:
    /path/to/analysis_results/*_NESTED_RESULTS.csv
    /path/to/analysis_results/*_cell_fluorescence_analysis.csv
    /path/to/structures_with_hierarchy.json   (optional — for acronyms)

OUTPUT:
    /path/to/merged_datasets/merged_dataset_[timestamp].csv
    /path/to/merged_datasets/merged_dataset_zones_[timestamp].csv

Usage:
    python 10_merge_nested_results.py

Requirements:
    pip install pandas numpy
"""

import re
import json
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

warnings.filterwarnings("ignore")

# =============================================================================
# CONFIGURATION
# =============================================================================
DIFFUSE_DIR     = Path("/path/to/analysis_results")
CELL_DIR        = Path("/path/to/analysis_results")
STRUCTURES_FILE = Path("/path/to/structures_with_hierarchy.json")
OUTPUT_DIR      = Path("/path/to/merged_datasets")

# Set to True to combine left and right hemispheres into one row.
AGGREGATE_HEMISPHERES = False

# Zones to exclude from the output (empty list = include all).
EXCLUDE_ZONES: list = []

# Treatment group identifier. Adjust the prefix-matching logic below if your
# animal ID naming convention differs.
REFERENCE_GROUP = "mScarlet"

# =============================================================================
# TREATMENT EXTRACTION
# =============================================================================

# Map animal ID prefixes to canonical treatment names.
# Order matters — more specific prefixes must appear before shorter ones.
TREATMENT_PREFIXES = [
    "C6ST1_ADAMTS15",
    "C6ST1_ADAMTS4",
    "ADAMTS4_MD",
    "ADAMTS15",
    "ADAMTS4",
    "C6ST1",
    "mScarlet",
]


def get_treatment(animal_id: str) -> str:
    """Extract treatment group from animal ID by prefix matching."""
    for prefix in TREATMENT_PREFIXES:
        if str(animal_id).startswith(prefix):
            return prefix
    return "Unknown"

# =============================================================================
# BRAIN STRUCTURES LOOKUP
# =============================================================================

def load_structures() -> tuple:
    """Load brain structure metadata for acronym and sphinx_id lookup."""
    if not STRUCTURES_FILE.exists():
        print(f"  Warning: {STRUCTURES_FILE} not found — acronyms will be empty.")
        return {}, {}

    with open(STRUCTURES_FILE) as f:
        structures = json.load(f)

    by_name = {}
    by_id   = {}
    for s in structures:
        info = {
            "acronym":             s.get("acronym", ""),
            "sphinx_id":           s.get("sphinx_id"),
            "parent_structure_id": s.get("parent_structure_id"),
            "depth":               s.get("depth"),
        }
        by_name[s["name"]] = info
        by_id[s["id"]]     = {"name": s["name"], **info}

    print(f"  Loaded {len(by_id)} brain structures.")
    return by_name, by_id

# =============================================================================
# ZONE / MARKER PARSING
# =============================================================================

def extract_zone_marker(cell_type_str: str) -> tuple:
    """
    Split cell_type strings like 'PNN_Core' or 'PV_Total' into (zone, marker).
    Returns ('Total', cell_type_str) as a fallback.
    """
    if pd.isna(cell_type_str):
        return "Total", cell_type_str
    parts  = str(cell_type_str).split("_", 1)
    marker = parts[0]
    zone   = parts[1] if len(parts) > 1 else "Total"
    return zone, marker

# =============================================================================
# PROCESS DIFFUSE FILES
# =============================================================================

def process_diffuse_file(filepath: Path, structures_by_id: dict) -> pd.DataFrame | None:
    """
    Process one _NESTED_RESULTS.csv file.

    Aggregates per-slice measurements to animal × brain_area × hemisphere × zone
    using area-weighted mean intensity.
    """
    m = re.match(r"(.+?)_NESTED_RESULTS\.csv", filepath.name)
    if not m:
        print(f"  Warning: cannot parse {filepath.name}")
        return None

    animal_id = m.group(1)
    df        = pd.read_csv(filepath)
    if df.empty:
        return None

    df["zone"]   = df["cell_type"].apply(lambda x: extract_zone_marker(x)[0])
    df["marker"] = df["cell_type"].apply(lambda x: extract_zone_marker(x)[1])

    if EXCLUDE_ZONES:
        df = df[~df["zone"].isin(EXCLUDE_ZONES)]
    df = df[df["marker"].isin(["PNN", "PV"])].copy()
    if df.empty:
        return None

    group_cols = (["brain_area", "regionID", "marker", "zone"]
                  if AGGREGATE_HEMISPHERES
                  else ["brain_area", "regionID", "marker", "zone", "hemisphere"])

    results = []
    for group_vals, grp in df.groupby(group_cols):
        if AGGREGATE_HEMISPHERES:
            brain_area, regionID, marker, zone = group_vals
            hemisphere = "both"
        else:
            brain_area, regionID, marker, zone, hemisphere = group_vals

        total_px   = grp["areaPx"].sum()
        total_mm2  = grp["areaMm2"].sum()
        if total_px == 0:
            continue

        wt_mean_int  = (grp["mean_intensity"]   * grp["areaPx"]).sum() / total_px
        wt_norm_int  = (grp["mean_norm_btm20"]  * grp["areaPx"]).sum() / total_px
        struct       = structures_by_id.get(int(regionID), {})

        results.append({
            "brain_area":  brain_area,
            "acronym":     struct.get("acronym", ""),
            "regionID":    int(regionID),
            "sphinx_id":   struct.get("sphinx_id"),
            "hemisphere":  hemisphere,
            "zone":        zone,
            "marker":      marker,
            "areaPx":      total_px,
            "areaMm2":     total_mm2,
            "diffFluo":    wt_mean_int,
            "avgPxIntensity": wt_norm_int,
            "animal_id":   animal_id,
            "treatment":   get_treatment(animal_id),
            "n_slices":    grp["slice_id"].nunique(),
        })

    return pd.DataFrame(results) if results else None

# =============================================================================
# PROCESS CELL FLUORESCENCE FILES
# =============================================================================

def process_cell_file(filepath: Path, structures_by_id: dict) -> pd.DataFrame | None:
    """
    Process one _cell_fluorescence_analysis.csv file.

    Aggregates to animal × brain_area × hemisphere, reporting cell count and
    mean normalised intensity (normalized_btm20).
    """
    df = pd.read_csv(filepath)
    if df.empty:
        return None

    animal_id = (df["mouse_id"].iloc[0] if "mouse_id" in df.columns
                 else filepath.stem.split("_cell_")[0])
    df = df[df["cell_type"].isin(["PNN", "PV"])].copy()
    if df.empty:
        return None

    group_cols = (["brain_area", "regionID", "cell_type"]
                  if AGGREGATE_HEMISPHERES
                  else ["brain_area", "regionID", "cell_type", "hemisphere"])

    results = []
    for group_vals, grp in df.groupby(group_cols):
        if AGGREGATE_HEMISPHERES:
            brain_area, regionID, cell_type = group_vals
            hemisphere = "both"
        else:
            brain_area, regionID, cell_type, hemisphere = group_vals

        int_col   = "normalized_btm20" if "normalized_btm20" in grp.columns else "mean"
        struct    = structures_by_id.get(int(regionID), {})

        results.append({
            "brain_area": brain_area,
            "acronym":    struct.get("acronym", ""),
            "regionID":   int(regionID),
            "sphinx_id":  struct.get("sphinx_id"),
            "hemisphere": hemisphere,
            "marker":     cell_type,
            "cell_count": len(grp),
            "intensity":  grp[int_col].mean(),
            "animal_id":  animal_id,
            "treatment":  get_treatment(animal_id),
            "n_slices":   grp["slice_id"].nunique() if "slice_id" in grp.columns else 1,
        })

    return pd.DataFrame(results) if results else None

# =============================================================================
# VALIDATION
# =============================================================================

def validate(diffuse_dir: Path, cell_dir: Path) -> tuple:
    diffuse_files = sorted([f for f in diffuse_dir.glob("*_NESTED_RESULTS.csv")
                            if not f.name.startswith("TEST_")])
    cell_files    = sorted([f for f in cell_dir.glob("*_cell_fluorescence_analysis.csv")
                            if not f.name.startswith("TEST_")])

    def animal_id(f, pattern):
        m = re.match(pattern, f.name)
        return m.group(1) if m else None

    d_animals = {animal_id(f, r"(.+?)_NESTED_RESULTS\.csv")          for f in diffuse_files}
    c_animals = {animal_id(f, r"(.+?)_cell_fluorescence_analysis\.csv") for f in cell_files}
    d_animals.discard(None)
    c_animals.discard(None)

    print(f"  Diffuse files:           {len(diffuse_files)}")
    print(f"  Cell files:              {len(cell_files)}")
    print(f"  Animals with both:       {len(d_animals & c_animals)}")
    if d_animals - c_animals:
        print(f"  Diffuse only:            {sorted(d_animals - c_animals)}")
    if c_animals - d_animals:
        print(f"  Cell only:               {sorted(c_animals - d_animals)}")

    return diffuse_files, cell_files

# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("MERGE NESTED RESULTS")
    print("=" * 70)

    _, structures_by_id = load_structures()

    print("\nValidating files ...")
    diffuse_files, cell_files = validate(DIFFUSE_DIR, CELL_DIR)

    # ── Diffuse data ──────────────────────────────────────────────────────────
    print("\nProcessing diffuse files ...")
    all_diffuse = []
    for f in diffuse_files:
        result = process_diffuse_file(f, structures_by_id)
        if result is not None:
            all_diffuse.append(result)
            print(f"  ✓ {f.stem}: {len(result)} rows")
        else:
            print(f"  ✗ {f.stem}: no valid data")

    df_diffuse = pd.concat(all_diffuse, ignore_index=True) if all_diffuse else pd.DataFrame()
    print(f"Total diffuse rows: {len(df_diffuse)}")

    # ── Cell data ─────────────────────────────────────────────────────────────
    print("\nProcessing cell files ...")
    all_cells = []
    for f in cell_files:
        result = process_cell_file(f, structures_by_id)
        if result is not None:
            all_cells.append(result)
            print(f"  ✓ {f.stem}: {len(result)} rows, "
                  f"{result['cell_count'].sum()} cells")
        else:
            print(f"  ✗ {f.stem}: no valid data")

    df_cells = pd.concat(all_cells, ignore_index=True) if all_cells else pd.DataFrame()
    print(f"Total cell rows: {len(df_cells)}")

    # ── Merge ─────────────────────────────────────────────────────────────────
    print("\nMerging ...")
    if not df_diffuse.empty:
        df_diffuse["staining"]   = df_diffuse["marker"].map({"PNN": "WFA", "PV": "PV"})
        df_diffuse["cell_type"]  = df_diffuse["marker"]
    if not df_cells.empty:
        df_cells["staining"]     = df_cells["marker"].map({"PNN": "WFA", "PV": "PV"})
        df_cells["cell_type"]    = df_cells["marker"]

    merge_keys = (["brain_area", "regionID", "animal_id", "staining"]
                  if AGGREGATE_HEMISPHERES
                  else ["brain_area", "regionID", "animal_id", "staining", "hemisphere"])

    if not df_diffuse.empty and not df_cells.empty:
        df_total       = df_diffuse[df_diffuse["zone"] == "Total"].copy()
        df_zones       = df_diffuse[df_diffuse["zone"] != "Total"].copy()
        df_merged      = df_total.merge(
            df_cells[merge_keys + ["cell_count", "intensity"]],
            on=merge_keys, how="outer"
        )
    elif not df_diffuse.empty:
        df_merged = df_diffuse.copy()
        df_zones  = pd.DataFrame()
    elif not df_cells.empty:
        df_merged         = df_cells.copy()
        df_merged["zone"] = "Total"
        df_zones          = pd.DataFrame()
    else:
        print("ERROR: no data to merge.")
        return

    # Fill treatment from animal_id where missing
    df_merged["treatment"] = df_merged.apply(
        lambda r: get_treatment(r["animal_id"]) if pd.isna(r.get("treatment")) else r["treatment"],
        axis=1,
    )

    # Derived metrics
    if "areaMm2" in df_merged.columns:
        df_merged["density"] = (df_merged["cell_count"] / df_merged["areaMm2"]).replace(
            [np.inf, -np.inf], np.nan
        )
    if "intensity" in df_merged.columns and "density" in df_merged.columns:
        df_merged["energy"] = df_merged["density"] * df_merged["intensity"]

    # ── Save ──────────────────────────────────────────────────────────────────
    out_main  = OUTPUT_DIR / f"merged_dataset_{timestamp}.csv"
    out_zones = OUTPUT_DIR / f"merged_dataset_zones_{timestamp}.csv"

    df_merged.to_csv(out_main, index=False)
    print(f"\nSaved: {out_main}  ({len(df_merged)} rows)")

    if not df_zones.empty:
        df_zones["treatment"] = df_zones["animal_id"].apply(get_treatment)
        df_zones.to_csv(out_zones, index=False)
        print(f"Saved: {out_zones}  ({len(df_zones)} rows)")

    print("\nDone.")


if __name__ == "__main__":
    main()
