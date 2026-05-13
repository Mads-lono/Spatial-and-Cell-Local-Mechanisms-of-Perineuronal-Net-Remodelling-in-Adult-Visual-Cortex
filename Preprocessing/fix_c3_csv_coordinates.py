#!/usr/bin/env python3
"""
diagnostics/fix_c3_csv_coordinates.py
=======================================
One-time coordinate correction for C3 Cellpose detection CSVs.

Background:
    The C3 detections (Counts_C3-*.csv) were run on the original unflipped
    TIFFs. Slices whose Net_Orientation_Change == "Mirrored" in the
    FINAL_TIFF_TRANSFORMATION_LOG.csv were subsequently horizontally flipped
    (and renamed) so that all images in the pipeline share a consistent
    orientation. The mask TIFFs were also correctly flipped at that stage.
    However, the Counts_C3 CSV files were only renamed — the Global_X
    coordinates inside them were never corrected.

    This script applies:
        Global_X_corrected = IMAGE_WIDTH - 1 - Global_X

    for every CSV whose (Animal, Final_ID) key maps to "Mirrored" in the
    transformation log.

Key detail:
    The lookup uses Final_ID (the post-rename slice ID on disk), NOT
    Original_ID. The CSV files on disk have already been renamed, so
    matching by Original_ID would target the wrong files.

Safety:
    Originals are backed up to {C3_DIR}/coord_fix_backup/ before modification.
    Already-backed-up files are skipped, making the script safe to rerun.

Usage:
    python diagnostics/fix_c3_csv_coordinates.py

    Set DRY_RUN = True to preview which files would be changed without
    modifying anything.

Requirements:
    pip install pandas
"""

import re
import shutil
from pathlib import Path

import pandas as pd

# =============================================================================
# CONFIGURATION
# =============================================================================
# Transformation log CSV
TRANSFORM_LOG = Path("/path/to/FINAL_TIFF_TRANSFORMATION_LOG.csv")

# Directory containing Counts_C3-*.csv files
C3_DIR = Path("/path/to/analysis_results")

# Backup directory (created automatically)
BACKUP_DIR = C3_DIR / "coord_fix_backup"

# hiRes image width in pixels — used for coordinate flip
# Confirm by checking a SegLabel TIFF: shape[1] gives the width
IMAGE_WIDTH = 16792

# Set True to preview without modifying any files
DRY_RUN = False

# =============================================================================
# HELPERS
# =============================================================================

def load_flip_lookup(log_path: Path) -> dict:
    """
    Build a {(animal, final_id): needs_flip} lookup from the transformation log.

    needs_flip is True when Net_Orientation_Change == "Mirrored".
    """
    log = pd.read_csv(log_path)
    log.columns = log.columns.str.strip()
    lookup = {}
    for _, row in log.iterrows():
        animal    = str(row["Animal"]).strip()
        final_id  = str(row["Final_ID"]).strip()
        needs_flip = str(row["Net_Orientation_Change"]).strip() == "Mirrored"
        lookup[(animal, final_id)] = needs_flip
    return lookup


def parse_c3_csv_name(fname: str):
    """
    Extract (animal, slice_id) from a filename like:
        Counts_C3-ADAMTS15_1_s001.csv   →  ('ADAMTS15_1', 's001')
        Counts_C3-C6ST1_ADAMTS4_2_s007.csv → ('C6ST1_ADAMTS4_2', 's007')

    Animal names may contain underscores, so the slice ID is matched from
    the end of the stem.
    """
    stem  = Path(fname).stem   # e.g. Counts_C3-ADAMTS15_1_s001
    match = re.match(r"Counts_C3-(.+)_(s\d{3})$", stem)
    if match:
        return match.group(1), match.group(2)
    return None, None

# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    print("=" * 60)
    print("FIX C3 CSV COORDINATES")
    if DRY_RUN:
        print("*** DRY RUN — no files will be modified ***")
    print("=" * 60)

    if not TRANSFORM_LOG.exists():
        print(f"ERROR: transformation log not found: {TRANSFORM_LOG}")
        return

    # Load flip lookup
    flip_lookup = load_flip_lookup(TRANSFORM_LOG)
    n_mirrored  = sum(flip_lookup.values())
    print(f"\nLoaded {len(flip_lookup)} records — {n_mirrored} marked Mirrored.\n")

    # Find all C3 Counts CSVs
    csv_files = sorted(C3_DIR.glob("Counts_C3-*.csv"))
    print(f"Found {len(csv_files)} Counts_C3 CSV files.\n")

    if not DRY_RUN:
        BACKUP_DIR.mkdir(exist_ok=True)

    stats = {
        "flipped":           0,
        "unchanged":         0,
        "already_done":      0,
        "not_in_log":        0,
        "parse_error":       0,
        "missing_column":    0,
    }

    for csv_path in csv_files:
        animal, slice_id = parse_c3_csv_name(csv_path.name)

        if animal is None:
            print(f"  PARSE ERROR: {csv_path.name}")
            stats["parse_error"] += 1
            continue

        key = (animal, slice_id)
        if key not in flip_lookup:
            stats["not_in_log"] += 1
            continue

        if not flip_lookup[key]:
            # Original orientation — no change needed
            stats["unchanged"] += 1
            continue

        # Check if already corrected (backup exists)
        backup_path = BACKUP_DIR / csv_path.name
        if backup_path.exists():
            stats["already_done"] += 1
            continue

        if DRY_RUN:
            print(f"  WOULD FLIP: {csv_path.name}  ({animal}, {slice_id})")
            stats["flipped"] += 1
            continue

        # Load, check, flip, save
        df = pd.read_csv(csv_path)

        if "Global_X" not in df.columns:
            print(f"  WARNING: No Global_X column — {csv_path.name}")
            stats["missing_column"] += 1
            continue

        shutil.copy2(csv_path, backup_path)
        df["Global_X"] = IMAGE_WIDTH - 1 - df["Global_X"]
        df.to_csv(csv_path, index=False)
        print(f"  FLIPPED: {csv_path.name}")
        stats["flipped"] += 1

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for key, val in stats.items():
        print(f"  {key:<22}: {val}")

    if DRY_RUN:
        print("\nDRY RUN complete — no files were modified.")
        print("Set DRY_RUN = False to apply corrections.")
    else:
        print(f"\nOriginals backed up to: {BACKUP_DIR}")
        print("Safe to rerun — already-corrected files are skipped.")


if __name__ == "__main__":
    main()
