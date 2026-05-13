#!/usr/bin/env python3
"""
diagnostics/mask_qc.py
=======================
Identifies injection mask false positives by analysing the distribution of
MaskLoose white-pixel percentages across all slices for each animal.

A large gap in the sorted percentage distribution indicates a bimodal
distribution: slices with a genuine injection zone (small, consistent
percentages) and slices where the mask fires on background (much larger
percentages). The midpoint of the largest gap is suggested as a cutoff.

Animals with many slices above the cutoff are likely injection failures and
should be excluded from analysis (via the exclusion file read by script 08).

INPUT:  /path/to/analysis_results/MaskLoose_C3-{animal}_s*.tif
OUTPUT: Printed report to stdout; no files written.

Usage:
    python diagnostics/mask_qc.py

Requirements:
    pip install numpy tifffile
"""

import numpy as np
from pathlib import Path
from tifffile import imread

# =============================================================================
# CONFIGURATION
# =============================================================================
# Directory containing MaskLoose and MaskStrict TIFFs
MASK_DIR = Path("/path/to/analysis_results")

# Number of top gaps to display in the report
TOP_GAPS = 5

# Percentage cutoff above which a slice is flagged as a likely false positive.
# Set to None to use the auto-detected cutoff (midpoint of largest gap).
MANUAL_CUTOFF = None   # e.g. 36.0

# =============================================================================
# QC ANALYSIS
# =============================================================================

def analyse_animal(animal_id: str, mask_dir: Path) -> None:
    """
    Print the MaskLoose white-pixel percentage distribution for one animal
    and suggest a false-positive cutoff based on the largest gap.
    """
    mask_files = sorted(mask_dir.glob(f"MaskLoose_C3-{animal_id}_s*.tif"))
    if not mask_files:
        print(f"  No MaskLoose files found for {animal_id}")
        return

    print(f"\n{'=' * 60}")
    print(f"Animal: {animal_id}  ({len(mask_files)} slices)")
    print(f"{'=' * 60}")

    percentages = []
    for fpath in mask_files:
        img      = imread(fpath)
        pct      = (np.count_nonzero(img) / img.size) * 100
        percentages.append(pct)
        print(f"  {fpath.name}: {pct:.2f}%")

    sorted_pcts = sorted(percentages)
    print(f"\nSorted values: {[f'{p:.2f}' for p in sorted_pcts]}")

    # Find gaps
    gaps = []
    for i in range(len(sorted_pcts) - 1):
        gap = sorted_pcts[i + 1] - sorted_pcts[i]
        gaps.append((gap, sorted_pcts[i], sorted_pcts[i + 1]))

    gaps.sort(reverse=True)
    print(f"\nTop {TOP_GAPS} gaps:")
    for gap, low, high in gaps[:TOP_GAPS]:
        print(f"  Gap {gap:.2f}%  between {low:.2f}% and {high:.2f}%")

    cutoff = MANUAL_CUTOFF
    if cutoff is None and gaps:
        best_gap, low, high = gaps[0]
        cutoff = (low + high) / 2

    if cutoff is not None:
        n_real   = sum(1 for p in percentages if p <= cutoff)
        n_false  = sum(1 for p in percentages if p >  cutoff)
        print(f"\n>>> Cutoff: {cutoff:.2f}%")
        print(f"    Real slices (≤ cutoff):          {n_real}")
        print(f"    Likely false positives (> cutoff): {n_false}")
        if n_false > len(mask_files) / 2:
            print(f"    *** MOST SLICES ARE FALSE POSITIVES — consider excluding {animal_id} ***")


def run(mask_dir: Path) -> None:
    """Discover all unique animal IDs from MaskLoose files and analyse each."""
    if not mask_dir.exists():
        print(f"Mask directory not found: {mask_dir}")
        return

    # Extract unique animal IDs from filenames like MaskLoose_C3-ANIMAL_sXXX.tif
    all_files   = list(mask_dir.glob("MaskLoose_C3-*.tif"))
    animal_ids  = set()
    for f in all_files:
        # Strip "MaskLoose_C3-" prefix and "_sXXX.tif" suffix
        stem = f.stem.replace("MaskLoose_C3-", "")
        # Animal ID is everything before the last _sXXX portion
        parts = stem.rsplit("_", 1)
        if len(parts) == 2 and parts[1].lower().startswith("s"):
            animal_ids.add(parts[0])

    if not animal_ids:
        print("No MaskLoose files found. Check MASK_DIR.")
        return

    print(f"Found {len(animal_ids)} animal(s).\n")

    for animal_id in sorted(animal_ids):
        analyse_animal(animal_id, mask_dir)

    print("\n" + "=" * 60)
    print("QC complete. Animals with predominantly false-positive masks")
    print("should be added to the exclusion file used by script 08.")


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    run(MASK_DIR)
