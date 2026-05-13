#!/usr/bin/env python3
"""
collect_all_results.py

Finds every CSV in the results directory, excludes diagnostics and large
raw tables, and concatenates all remaining tables into a single CSV with
two prefix columns: source_dir and source_file.

Output: ALL_RESULTS_COMBINED.csv

Usage:
  python collect_all_results.py

Requirements:
  pip install pandas
"""

import pandas as pd
from pathlib import Path

BASE = Path("/path/to/results")
OUT  = BASE / "ALL_RESULTS_COMBINED.csv"

# Files whose names contain any of these strings are skipped.
EXCLUDE_PATTERNS = [
    "diagnostic",
    "section_data",                     # large raw section-level table
    "atlas_enwrapment",                 # large raw atlas table
    "slice_enwrapment",                 # large raw slice table
    "zone_enwrapment",                  # large raw zone table
    "cell_virus_tags",                  # per-cell tag table — millions of rows
    "pv_cells_viral",                   # per-cell viral load — millions of rows
    "animal_metrics_wide",
    "metric_correlation_matrix",
    "analysis2_section_heterogeneity",  # per-section raw data
    "composite_hemisphere_data",
    "composite_section_data",
    "composite_validation_animal",
    "analysis_merged",                  # per-animal merged table
    "scarlet_expression_animal",
    "binned_gradient_with_baseline",    # raw bin data
    "cell_level_distances",             # per-cell distance table
    "off_target_atlas_enwrapment",      # large atlas off-target raw
    "partA_triple_colocalization",      # raw counts table
    "section_virus_enwrapment",         # raw section table
    "partC_adamts4_section_data",       # raw section data
    "zone_density_all",                 # full zone density raw table
    "plot_summary",                     # figure summary, not model results
    "expression_enwrapment_merged",     # raw merged table
    "group_means_ratio",
    "group_means_all",
    "group_means_noswapped",
    "partA_cv_table",
    "partA_expression_descriptive",
    "partD_zone_gradient_by_treatment",
    "partD_core_outside_ratio",
    "animal_summary",
    "animal_zone_summary",
    "partG_animal_layer_summary",
    "per_area_contrasts",
    "fdr_comparison",
    "effect_size_comparison",
    "one_sample_tests",
    "contrasts_raw_ratio",
    "posterior_comparison",
]


def should_skip(fname: str) -> bool:
    fl = fname.lower()
    return any(pat in fl for pat in EXCLUDE_PATTERNS)


records = []
skipped = []

csv_files = sorted(BASE.rglob("*.csv"))
print(f"Found {len(csv_files)} CSV files total.\n")

for path in csv_files:
    rel   = path.relative_to(BASE)
    parts = rel.parts
    fname = path.name

    if should_skip(fname):
        skipped.append(str(rel))
        continue

    dir_parts  = list(parts[:-1])
    source_dir = "/".join(dir_parts)

    try:
        df = pd.read_csv(path, dtype=str)  # read all as str to avoid type clashes
        if df.empty:
            skipped.append(f"{rel} [empty]")
            continue
        df.insert(0, "source_file", fname)
        df.insert(0, "source_dir",  source_dir)
        records.append(df)
        print(f"  OK   {rel}  ({len(df)} rows, {len(df.columns)-2} cols)")
    except Exception as e:
        skipped.append(f"{rel} [ERROR: {e}]")
        print(f"  SKIP {rel}  [{e}]")

print(f"\n{'='*60}")
print(f"Loaded {len(records)} files.")
print(f"Skipped {len(skipped)} files.\n")

if records:
    combined = pd.concat(records, axis=0, ignore_index=True, sort=False)
    combined.to_csv(OUT, index=False)
    print(f"Written: {OUT}")
    print(f"  Total rows: {len(combined):,}")
    print(f"  Total columns: {len(combined.columns)}")
else:
    print("No records to write.")

print("\nSkipped files:")
for s in skipped:
    print(f"  {s}")
