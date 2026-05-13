#!/usr/bin/env python3
"""
07_csvsplitter.py
=================
Filters CPN detection output (localizations.csv) and splits it into
per-slice coordinate CSV files for use by the Diffu_Spot pipeline.

Steps:
  1. Read a CPN localizations.csv file.
  2. Keep only rows with rescore >= RESCORE_THRESHOLD (default 0.5).
  3. For each unique image name:
       - Strip the .tif extension and "_enhanced" suffix.
       - Rename channel suffixes to match the DATA/counts/ naming convention:
           _C0 → _A1_1-cells_C2   (PNN channel from CPN)
           _C1 → _A1_1-cells_C1   (PV channel from CPN)
       - Save x, y coordinates only to OUTPUT_DIR/{savename}.csv.

This is run once per animal after CPN inference is complete.

INPUT:  /path/to/CPN_output/{animal}/localizations.csv
OUTPUT: /path/to/DATA/{animal}/counts/{slice_name}.csv

Usage:
    python 07_csvsplitter.py

Requirements:
    pip install pandas
"""

import os
import pandas as pd

# =============================================================================
# CONFIGURATION
# =============================================================================
# Path to the CPN localizations CSV for one animal
INPUT_CSV = "/path/to/CPN_output/ANIMAL_1/localizations.csv"

# Directory where per-slice coordinate CSVs will be saved
OUTPUT_DIR = "/path/to/DATA/ANIMAL_1/counts"

# Minimum rescore confidence to retain a detection
RESCORE_THRESHOLD = 0.5

# =============================================================================
# SPLITTING LOGIC
# =============================================================================

def split_csv(input_csv: str, output_dir: str, rescore_threshold: float = 0.5) -> None:
    """
    Filter localizations by rescore and split into per-slice coordinate CSVs.

    The channel renaming (_C0/_C1 → _A1_1-cells_C2/_C1) preserves the naming
    convention expected by the Diffu_Spot script when loading counts from
    DATA/{animal}/counts/.
    """
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(input_csv)
    df_filtered = df[df["rescore"] >= rescore_threshold]

    print(f"Loaded {len(df)} detections, {len(df_filtered)} pass rescore >= {rescore_threshold}")

    saved = 0
    for img_name in df_filtered["imgName"].unique():
        savename = img_name.replace(".tif", "").replace("_enhanced", "")

        if "_C0" in savename:
            savename = savename.replace("_C0", "_A1_1-cells_C2")
        elif "_C1" in savename:
            savename = savename.replace("_C1", "_A1_1-cells_C1")

        subset = (
            df_filtered[df_filtered["imgName"] == img_name][["X", "Y"]]
            .rename(columns={"X": "x", "Y": "y"})
        )

        out_path = os.path.join(output_dir, f"{savename}.csv")
        subset.to_csv(out_path, index=False)
        saved += 1

    print(f"Saved {saved} per-slice CSV files to {output_dir}")


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    split_csv(INPUT_CSV, OUTPUT_DIR, RESCORE_THRESHOLD)
