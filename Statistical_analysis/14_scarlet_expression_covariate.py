
#!/usr/bin/env python3
"""
14_scarlet_expression_covariate.py
====================================
Uses Cellpose-segmented mScarlet (C3) cell detections to build a per-animal
transgene expression index, then tests whether expression level predicts
PNN outcome measures — separating enzyme activity from injection microenvironment.

Pipeline:
  1. Load all Counts_C3-*.csv files from results directory
  2. Apply do_not_use exclusions — convert Original_IDs to Final_IDs via log
  3. Load FINAL_TIFF_TRANSFORMATION_LOG to build Final_ID-keyed lookups
  4. Apply name corrections (C6ST1_ADAMTS4 → C6ST1_ADAMTS15)
  5. Apply swap correction (mScarlet_4 ↔ ADAMTS4_MD_4 treatment labels)
  6. Exclude injection failure animal
  7. Aggregate to per-animal expression metrics:
       - n_cells_core: mScarlet cell count in Core zone
       - mean_intensity_core: mean mScarlet intensity in Core zone
       - total_intensity_core: total mScarlet signal in Core zone
  8. Merge with enwrapment outcomes from cells_with_zones
  9. Three analytical parts:
       Part A - Descriptive: expression variability across animals
       Part B - Within-group dose-response: expression vs enwrapment
       Part C - Covariate LMM: does treatment effect survive adjustment
       Part D - mScarlet controls: injection microenvironment baseline
       Part E - ADAMTS4 vs ADAMTS15 expression comparison
"""

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import gc
from tqdm import tqdm
from joblib import Parallel, delayed

warnings.filterwarnings('ignore')

# ── Paths ─────────────────────────────────────────────────────────────────────
C3_RESULTS_DIR  = Path("/path/to/c3_results")
TRANSFORM_LOG   = Path("/path/to/FINAL_TIFF_TRANSFORMATION_LOG.csv")
DO_NOT_USE      = Path("/path/to/do_not_use.csv")
CELLS_CSV       = Path("/path/to/analysis_results/cells_with_zones.csv")
MERGED_ZONES    = Path("/path/to/merged_datasets/merged_dataset_zones.csv")
OUT_DIR         = Path("/path/to/results/14_scarlet_expression")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Constants ─────────────────────────────────────────────────────────────────
IMG_WIDTH       = 16792   # image width in pixels
TREATMENT_ORDER = ["mScarlet","ADAMTS4","ADAMTS4_MD","ADAMTS15","C6ST1","C6ST1_ADAMTS15"]
EXCLUDE         = ""   # any mouse_id containing this string will be excluded (e.g. injection failure)
ALPHA           = 0.05

PALETTE = {
    "mScarlet":        "#888888",
    "ADAMTS4":         "#4e9af1",
    "ADAMTS4_MD":      "#f17c4e",
    "ADAMTS15":        "#4ef196",
    "C6ST1":           "#c44ef1",
    "C6ST1_ADAMTS15":  "#f1c44e",
}

# ── Helpers ───────────────────────────────────────────────────────────────────
def parse_c3_filename(fname):
    """Parse Counts_C3-{treatment}_{animal}_{slice}.csv
    Returns (treatment_raw, animal_num, final_slice).
    Slice extracted from filename IS the Final_ID — files were already renamed
    by apply_transformations_to_masks.ipynb cell 9.
    """
    m = re.match(r"Counts_C3-(.+)_(\d+)_(s\d+)\.csv", fname)
    if m:
        return m.group(1), m.group(2), m.group(3)
    return None, None, None


# =============================================================================
# STEP 1: Load transformation log and do_not_use
# =============================================================================

print("Loading transformation log and exclusion list...")
tlog = pd.read_csv(TRANSFORM_LOG)
tlog.columns = tlog.columns.str.strip()

# Build orig→final mapping per animal for do_not_use conversion
orig_to_final = {}  # (animal_id, orig_slice) → final_slice
final_valid   = set()  # (animal_id, final_slice) — all valid final slices

for _, row in tlog.iterrows():
    animal_id   = row["Animal"]
    orig_slice  = row["Original_ID"]
    final_slice = row["Final_ID"]
    orig_to_final[(animal_id, orig_slice)] = final_slice
    final_valid.add((animal_id, final_slice))

# do_not_use contains filenames like MaskLoose_C3-ADAMTS15_1_s051_thumb.png
# These use Original_IDs — convert to Final_IDs for correct exclusion
dnu = pd.read_csv(DO_NOT_USE)
dnu.columns = ["filename"]
dnu["filename"] = dnu["filename"].str.strip()

def parse_dnu(fname):
    m = re.match(r"MaskLoose_C3-(.+)_(\d+)_(s\d+)_thumb\.png", fname)
    if m:
        return f"{m.group(1)}_{m.group(2)}", m.group(3)
    return None, None

dnu[["animal_id", "orig_slice"]] = dnu["filename"].apply(
    lambda x: pd.Series(parse_dnu(x)))
dnu = dnu.dropna(subset=["animal_id"])

# Convert do_not_use from Original_ID → Final_ID
exclude_set = set()
for _, row in dnu.iterrows():
    key = (row["animal_id"], row["orig_slice"])
    final = orig_to_final.get(key)
    if final:
        exclude_set.add((row["animal_id"], final))
    else:
        # Slice not in log — exclude by original ID as fallback
        exclude_set.add((row["animal_id"], row["orig_slice"]))

print(f"  Exclusion list: {len(exclude_set)} slice/animal pairs (Final_IDs)")

# =============================================================================
# STEP 2: Build Final_ID-keyed lookup (no flip needed — done on disk)
# =============================================================================
# Coordinates in CSV files are already corrected. No X flip applied here.

# =============================================================================
# STEP 3: Load MaskStrict and MaskLoose per slice for zone assignment
# =============================================================================
# We need to know, for each mScarlet cell, whether it falls in Core, Penumbra,
# or Outside. Core = inside MaskStrict, Penumbra = MaskLoose only, Outside = neither.
# We load masks on demand per slice to avoid memory overload.

def get_zone(x, y, animal_id, final_slice, img_width=IMG_WIDTH):
    """
    Load MaskStrict and MaskLoose for a given animal/final_slice and return
    zone label for a cell at (x, y). Coordinates are already in corrected
    (post-flip) space. Mask files are named with Final_IDs.
    Returns 'Core', 'Penumbra', or 'Outside'.
    Caches masks in module-level dict to avoid reloading.
    """
    key = (animal_id, final_slice)
    if key not in _mask_cache:
        strict_path = C3_RESULTS_DIR / f"MaskStrict_C3-{animal_id}_{final_slice}.tif"
        loose_path  = C3_RESULTS_DIR / f"MaskLoose_C3-{animal_id}_{final_slice}.tif"
        try:
            _mask_cache[key] = (
                tifffile.imread(str(strict_path)).astype(bool),
                tifffile.imread(str(loose_path)).astype(bool)
            )
        except Exception:
            _mask_cache[key] = (None, None)
    strict_mask, loose_mask = _mask_cache[key]
    if strict_mask is None:
        return "Outside"
    h, w = strict_mask.shape
    yi, xi = int(np.clip(y, 0, h-1)), int(np.clip(x, 0, w-1))
    if strict_mask[yi, xi]:
        return "Core"
    elif loose_mask[yi, xi]:
        return "Penumbra"
    return "Outside"

_mask_cache = {}

# =============================================================================
# STEP 4: Load all Counts CSVs
# =============================================================================
print("Loading all Counts CSVs (this may take a minute)...")
all_records = []
skipped_dnu    = 0
skipped_no_map = 0

count_files = sorted(C3_RESULTS_DIR.glob("Counts_C3-*.csv"))

for fpath in count_files:
    treatment_raw, animal_num, final_slice = parse_c3_filename(fpath.name)
    if treatment_raw is None:
        continue

    animal_id_raw = f"{treatment_raw}_{animal_num}"

    # do_not_use check — keyed on Final_ID
    if (animal_id_raw, final_slice) in exclude_set:
        skipped_dnu += 1
        continue

    # Verify slice is in transformation log
    if (animal_id_raw, final_slice) not in final_valid:
        skipped_no_map += 1
        continue

    try:
        df = pd.read_csv(fpath)
    except Exception:
        continue

    if df.empty or "Mean_Intensity" not in df.columns:
        continue

    # NOTE: No X flip applied here — Global_X already corrected in files
    # by fix_c3_csv_coordinates.py

    df["animal_id_raw"] = animal_id_raw
    df["final_slice"]   = final_slice
    df["treatment_raw"] = treatment_raw
    df["animal_num"]    = animal_num

    all_records.append(df)

cells_c3 = pd.concat(all_records, ignore_index=True)
print(f"  Loaded: {len(cells_c3)} mScarlet cells")
print(f"  Excluded (do_not_use): {skipped_dnu} slices")
print(f"  Skipped (not in log):  {skipped_no_map} slices")

# =============================================================================
# STEP 5: Apply corrections
# =============================================================================
# Name correction: C6ST1_ADAMTS4 → C6ST1_ADAMTS15
cells_c3["treatment_raw"] = cells_c3["treatment_raw"].str.replace(
    "C6ST1_ADAMTS4", "C6ST1_ADAMTS15", regex=False)

# Build mouse_id
cells_c3["mouse_id"] = cells_c3["treatment_raw"] + "_" + cells_c3["animal_num"].astype(str)

# Swap correction
cells_c3.loc[cells_c3["mouse_id"] == "mScarlet_4",   "treatment_raw"] = "ADAMTS4_MD"
cells_c3.loc[cells_c3["mouse_id"] == "ADAMTS4_MD_4", "treatment_raw"] = "mScarlet"

# Exclude injection failure
cells_c3 = cells_c3[cells_c3["mouse_id"] != EXCLUDE]

# Rename for clarity
cells_c3 = cells_c3.rename(columns={
    "treatment_raw": "treatment",
    "Mean_Intensity": "scarlet_intensity",
    "Global_X": "x",
    "Global_Y": "y"
})

print(f"\nAfter corrections: {len(cells_c3)} cells")
print("Animals:", sorted(cells_c3["mouse_id"].unique()))



# =============================================================================
# STEP 6: Parallelized Zone Assignment with Progress Bar
# =============================================================================
print(f"\nAssigning zones using 12 CPU cores...")

def process_single_slice(group_keys, indices, x_vals, y_vals):
    """
    Function to be run in parallel worker processes.
    Loads MaskStrict/MaskLoose using Final_ID — masks are already renamed.
    """
    animal_id, final_slice = group_keys
    strict_path = C3_RESULTS_DIR / f"MaskStrict_C3-{animal_id}_{final_slice}.tif"
    loose_path  = C3_RESULTS_DIR / f"MaskLoose_C3-{animal_id}_{final_slice}.tif"
    
    # Initialize with default
    labels = np.full(len(indices), "Outside", dtype=object)
    
    if not strict_path.exists() or not loose_path.exists():
        return indices, labels

    try:
        # Load masks as bools to save space
        # Use TiffFile context manager for clean file closing
        with tifffile.TiffFile(str(strict_path)) as tif:
            strict_mask = tif.asarray().astype(bool)
        with tifffile.TiffFile(str(loose_path)) as tif:
            loose_mask = tif.asarray().astype(bool)
        
        h, w = strict_mask.shape
        
        # Vectorized coordinate lookup
        xs = np.clip(x_vals.astype(int), 0, w - 1)
        ys = np.clip(y_vals.astype(int), 0, h - 1)
        
        # Priority: Strict (Core) > Loose (Penumbra) > Outside
        is_strict = strict_mask[ys, xs]
        is_loose  = loose_mask[ys, xs]
        
        labels[is_loose]  = "Penumbra"
        labels[is_strict] = "Core"
        
        # --- Critical Memory Cleanup ---
        del strict_mask
        del loose_mask
        gc.collect() 
        # -------------------------------
        
        return indices, labels
        
    except Exception:
        return indices, labels

# 1. Prepare tasks
slice_groups = cells_c3.groupby(["animal_id_raw", "final_slice"])
tasks = [
    delayed(process_single_slice)(
        name, idx, cells_c3.loc[idx, "x"].values, cells_c3.loc[idx, "y"].values
    ) 
    for name, idx in slice_groups.groups.items()
]

# 2. Execute with progress bar
# n_jobs=12 to utilize your hardware while leaving overhead for the OS
results = Parallel(n_jobs=12)(
    tqdm(tasks, desc="Processing Slices", total=len(tasks))
)

# 3. Update main DataFrame
cells_c3["zone_c3"] = "Outside"
for indices, labels in results:
    cells_c3.loc[indices, "zone_c3"] = labels

# Final clear
del results
gc.collect()

print(f"\nProcessing complete. Total mScarlet cells assigned: {len(cells_c3)}")

# =============================================================================
# STEP 7: Aggregate to animal level — expression indices
# =============================================================================
print("\nAggregating to animal level...")

# Per animal × zone
agg_zone = cells_c3.groupby(["mouse_id","treatment","zone_c3"]).agg(
    n_cells     = ("scarlet_intensity", "count"),
    mean_intens = ("scarlet_intensity", "mean"),
    total_intens= ("scarlet_intensity", "sum"),
    median_intens=("scarlet_intensity","median")
).reset_index()

# Focus on Core zone as primary expression proxy
core_expr = agg_zone[agg_zone["zone_c3"] == "Core"].copy()
core_expr = core_expr.rename(columns={
    "n_cells":      "n_scarlet_core",
    "mean_intens":  "mean_scarlet_core",
    "total_intens": "total_scarlet_core",
    "median_intens":"median_scarlet_core"
})
core_expr = core_expr.drop(columns=["zone_c3"])

# Overall per animal
agg_animal = cells_c3.groupby(["mouse_id","treatment"]).agg(
    n_scarlet_total      = ("scarlet_intensity","count"),
    mean_scarlet_total   = ("scarlet_intensity","mean"),
    total_scarlet_total  = ("scarlet_intensity","sum")
).reset_index()

expr_animal = agg_animal.merge(core_expr, on=["mouse_id","treatment"], how="left")
expr_animal["treatment"] = pd.Categorical(
    expr_animal["treatment"], categories=TREATMENT_ORDER, ordered=True)
expr_animal = expr_animal.sort_values(["treatment","mouse_id"])

expr_animal.to_csv(OUT_DIR / "scarlet_expression_animal.csv", index=False)
print("\nExpression summary per animal:")
print(expr_animal[["mouse_id","treatment","n_scarlet_core",
                    "mean_scarlet_core","total_scarlet_core"]].to_string(index=False))

# =============================================================================
# STEP 8: Load PNN outcomes
# =============================================================================
print("\nLoading PNN outcome data...")

# Enwrapment from cells_with_zones
cells_pnn = pd.read_csv(CELLS_CSV)
# Corrections
cells_pnn["mouse_id"]  = cells_pnn["mouse_id"].str.replace(
    "C6ST1_ADAMTS4_", "C6ST1_ADAMTS15_", regex=False)
cells_pnn["treatment"] = cells_pnn["treatment"].str.replace(
    "C6ST1_ADAMTS4$", "C6ST1_ADAMTS15", regex=True)
cells_pnn.loc[cells_pnn["mouse_id"] == "mScarlet_4",   "treatment"] = "ADAMTS4_MD"
cells_pnn.loc[cells_pnn["mouse_id"] == "ADAMTS4_MD_4", "treatment"] = "mScarlet"
cells_pnn = cells_pnn[cells_pnn["mouse_id"] != EXCLUDE]

# Compute enwrapment per animal (ipsilateral, visual cortex, all zones)
from scipy.spatial import cKDTree
COLOC_THRESH = 30

enwrap_records = []
for (mid, treat, sl), grp in cells_pnn[cells_pnn["hemisphere"] == "left"].groupby(
        ["mouse_id","treatment","slice_id"]):
    pv  = grp[grp["cell_type"] == "PV"][["x_hires","y_hires"]].dropna()
    pnn = grp[grp["cell_type"] == "PNN"][["x_hires","y_hires"]].dropna()
    if len(pv) < 5 or len(pnn) < 5:
        continue
    tree = cKDTree(pnn.values)
    dists, _ = tree.query(pv.values, k=1)
    frac = (dists <= COLOC_THRESH).mean()
    enwrap_records.append({
        "mouse_id": mid, "treatment": treat,
        "slice_id": sl, "frac_enwrapped": frac, "n_pv": len(pv)
    })

enwrap_df = pd.DataFrame(enwrap_records)
enwrap_animal = enwrap_df.groupby(["mouse_id","treatment"])["frac_enwrapped"].mean().reset_index()
enwrap_animal.columns = ["mouse_id","treatment","frac_enwrapped"]

# WFA intensity from merged_zones
mz = pd.read_csv(MERGED_ZONES)
mz = mz.rename(columns={"animal_id":"mouse_id"})
mz["mouse_id"]  = mz["mouse_id"].str.replace("C6ST1_ADAMTS4_","C6ST1_ADAMTS15_",regex=False)
mz["treatment"] = mz["treatment"].str.replace("C6ST1_ADAMTS4$","C6ST1_ADAMTS15",regex=True)
mz.loc[mz["mouse_id"] == "mScarlet_4",   "treatment"] = "ADAMTS4_MD"
mz.loc[mz["mouse_id"] == "ADAMTS4_MD_4", "treatment"] = "mScarlet"
mz = mz[mz["mouse_id"] != EXCLUDE]

is_vis = mz["brain_area"].str.lower().str.contains("visual|vis", na=False)
wfa_vis = mz[(mz["staining"] == "WFA") & is_vis]
avgpx_animal = wfa_vis.groupby(["mouse_id","treatment"])["avgPxIntensity"].mean().reset_index()

# Merge all outcomes with expression
outcomes = enwrap_animal.merge(avgpx_animal, on=["mouse_id","treatment"], how="outer")
analysis_df = outcomes.merge(expr_animal, on=["mouse_id","treatment"], how="inner")
analysis_df["treatment"] = pd.Categorical(
    analysis_df["treatment"], categories=TREATMENT_ORDER, ordered=True)
analysis_df.to_csv(OUT_DIR / "analysis_merged.csv", index=False)
print(f"\nMerged analysis table: {len(analysis_df)} animals")

# =============================================================================
# PART A: Descriptive — expression variability
# =============================================================================
print("\n" + "="*70)
print("PART A: Expression variability across animals")
print("="*70)

desc = analysis_df.groupby("treatment")[
    ["n_scarlet_core","mean_scarlet_core","frac_enwrapped"]
].agg(["mean","std","min","max"]).round(3)
print(desc.to_string())
desc.to_csv(OUT_DIR / "partA_expression_descriptive.csv")

# CV within treatment groups
cv_table = analysis_df.groupby("treatment")["n_scarlet_core"].agg(
    mean="mean", std="std").assign(
    cv=lambda x: x["std"] / x["mean"] * 100).round(2)
print("\nCoefficient of variation — Core mScarlet cell count:")
print(cv_table.to_string())
cv_table.to_csv(OUT_DIR / "partA_cv_table.csv")

# =============================================================================
# PART B: Within-group dose-response
# =============================================================================
print("\n" + "="*70)
print("PART B: Within-group dose-response (expression vs enwrapment)")
print("="*70)

results_B = []
for treat in TREATMENT_ORDER:
    sub = analysis_df[analysis_df["treatment"] == treat].dropna(
        subset=["n_scarlet_core","frac_enwrapped"])
    if len(sub) < 3:
        continue
    r, p = stats.pearsonr(sub["n_scarlet_core"], sub["frac_enwrapped"])
    r_int, p_int = stats.pearsonr(sub["mean_scarlet_core"], sub["frac_enwrapped"])
    results_B.append({
        "treatment": treat, "n": len(sub),
        "r_count_vs_enwrap": round(r, 3),
        "p_count_vs_enwrap": round(p, 4),
        "r_intensity_vs_enwrap": round(r_int, 3),
        "p_intensity_vs_enwrap": round(p_int, 4)
    })
    print(f"  {treat}: r(count)={r:.3f} p={p:.4f} | r(intensity)={r_int:.3f} p={p_int:.4f} | n={len(sub)}")

results_B_df = pd.DataFrame(results_B)
results_B_df.to_csv(OUT_DIR / "partB_within_group_correlations.csv", index=False)

# =============================================================================
# PART C: Covariate LMM — does treatment effect survive expression adjustment?
# =============================================================================
print("\n" + "="*70)
print("PART C: Treatment effect with expression as covariate")
print("="*70)
print("(Animal-level lm — n too small for LMM here)")

try:
    import statsmodels.formula.api as smf

    dat = analysis_df.dropna(subset=["frac_enwrapped","n_scarlet_core"]).copy()
    dat["treatment_cat"] = pd.Categorical(dat["treatment"],
                                           categories=TREATMENT_ORDER)
    dat["log_scarlet"] = np.log1p(dat["n_scarlet_core"])

    # Model 1: treatment only
    m1 = smf.ols("frac_enwrapped ~ C(treatment, Treatment('mScarlet'))",
                  data=dat).fit()
    # Model 2: treatment + expression covariate
    m2 = smf.ols("frac_enwrapped ~ C(treatment, Treatment('mScarlet')) + log_scarlet",
                  data=dat).fit()

    print("\n-- Model 1: treatment only --")
    print(m1.summary().tables[1])
    print("\n-- Model 2: treatment + log(scarlet cells) --")
    print(m2.summary().tables[1])
    print(f"\nR² change: {m1.rsquared:.3f} → {m2.rsquared:.3f}")
    print(f"Scarlet covariate p = {m2.pvalues.get('log_scarlet', float('nan')):.4f}")

    # Compare treatment coefficients before/after covariate
    treat_terms = [t for t in m1.params.index if "treatment" in t]
    coef_compare = pd.DataFrame({
        "term": treat_terms,
        "beta_no_cov":   [m1.params[t] for t in treat_terms],
        "beta_with_cov": [m2.params.get(t, np.nan) for t in treat_terms],
        "p_no_cov":      [m1.pvalues[t] for t in treat_terms],
        "p_with_cov":    [m2.pvalues.get(t, np.nan) for t in treat_terms]
    })
    coef_compare["beta_shift_pct"] = (
        (coef_compare["beta_with_cov"] - coef_compare["beta_no_cov"]) /
        coef_compare["beta_no_cov"].abs() * 100).round(1)
    print("\n-- Coefficient stability --")
    print(coef_compare.round(4).to_string(index=False))
    coef_compare.to_csv(OUT_DIR / "partC_covariate_stability.csv", index=False)

    # Expression covariate model within treated groups only
    treated = dat[dat["treatment"].isin(["ADAMTS4","ADAMTS4_MD","ADAMTS15",
                                          "C6ST1","C6ST1_ADAMTS15"])].copy()
    m3 = smf.ols("frac_enwrapped ~ log_scarlet", data=treated).fit()
    print(f"\n-- Across treated groups: enwrap ~ log(scarlet) --")
    print(f"  r² = {m3.rsquared:.3f}, p(scarlet) = {m3.pvalues['log_scarlet']:.4f}")
    m3_res = pd.DataFrame({
        "term": m3.params.index,
        "coef": m3.params.values,
        "p":    m3.pvalues.values
    })
    m3_res.to_csv(OUT_DIR / "partC_treated_expression_lm.csv", index=False)

except ImportError:
    print("statsmodels not available — skipping Part C")

# =============================================================================
# PART D: mScarlet controls — injection microenvironment baseline
# =============================================================================
print("\n" + "="*70)
print("PART D: mScarlet zone gradient (injection effect baseline)")
print("="*70)

# For mScarlet animals, compute enwrapment by zone
mscarlet_cells = cells_pnn[
    (cells_pnn["treatment"] == "mScarlet") &
    (cells_pnn["hemisphere"] == "left")
].copy()

is_vis_pnn = mscarlet_cells["brain_area"].str.lower().str.contains("visual|vis", na=False)
mscarlet_vis = mscarlet_cells[is_vis_pnn]

ms_enwrap = []
for (mid, sl, z), grp in mscarlet_vis.groupby(["mouse_id","slice_id","zone"]):
    pv  = grp[grp["cell_type"] == "PV"][["x_hires","y_hires"]].dropna()
    pnn = grp[grp["cell_type"] == "PNN"][["x_hires","y_hires"]].dropna()
    if len(pv) < 5 or len(pnn) < 5:
        continue
    tree = cKDTree(pnn.values)
    dists, _ = tree.query(pv.values, k=1)
    frac = (dists <= COLOC_THRESH).mean()
    ms_enwrap.append({"mouse_id": mid, "zone": z, "frac_enwrapped": frac})

ms_zone = pd.DataFrame(ms_enwrap).groupby(
    ["mouse_id","zone"])["frac_enwrapped"].mean().reset_index()
ms_zone_mean = ms_zone.groupby("zone")["frac_enwrapped"].agg(
    ["mean","std","count"]).round(4)
print("\nmScarlet zone gradient (Core/Penumbra/Outside enwrapment):")
print(ms_zone_mean.to_string())

# Compare with treated groups
zone_by_treat = []
is_vis_all = cells_pnn["brain_area"].str.lower().str.contains("visual|vis", na=False)
vis_cells = cells_pnn[is_vis_all & (cells_pnn["hemisphere"] == "left")]

for (mid, treat, sl, z), grp in vis_cells.groupby(
        ["mouse_id","treatment","slice_id","zone"]):
    pv  = grp[grp["cell_type"] == "PV"][["x_hires","y_hires"]].dropna()
    pnn = grp[grp["cell_type"] == "PNN"][["x_hires","y_hires"]].dropna()
    if len(pv) < 5 or len(pnn) < 5:
        continue
    tree = cKDTree(pnn.values)
    dists, _ = tree.query(pv.values, k=1)
    frac = (dists <= COLOC_THRESH).mean()
    zone_by_treat.append({
        "mouse_id": mid, "treatment": treat,
        "zone": z, "frac_enwrapped": frac
    })

zone_df = pd.DataFrame(zone_by_treat)
zone_means = zone_df.groupby(["treatment","zone"])["frac_enwrapped"].agg(
    ["mean","std"]).round(4)
print("\nZone enwrapment by treatment (Core/Outside gradient):")
print(zone_means.to_string())
zone_means.to_csv(OUT_DIR / "partD_zone_gradient_by_treatment.csv")

# Core/Outside ratio per animal (removes baseline variability)
zone_pivot = zone_df.groupby(["mouse_id","treatment","zone"])[
    "frac_enwrapped"].mean().unstack("zone").reset_index()
if "Core" in zone_pivot.columns and "Outside" in zone_pivot.columns:
    zone_pivot["core_outside_ratio"] = zone_pivot["Core"] / zone_pivot["Outside"]
    ratio_summary = zone_pivot.groupby("treatment")[
        "core_outside_ratio"].agg(["mean","std"]).round(4)
    print("\nCore/Outside enwrapment ratio per treatment:")
    print(ratio_summary.to_string())
    zone_pivot.to_csv(OUT_DIR / "partD_core_outside_ratio.csv", index=False)

# =============================================================================
# PART E: ADAMTS4 vs ADAMTS15 expression comparison
# =============================================================================
print("\n" + "="*70)
print("PART E: ADAMTS4 vs ADAMTS15 — expression levels compared")
print("="*70)

a4  = analysis_df[analysis_df["treatment"] == "ADAMTS4"]
a15 = analysis_df[analysis_df["treatment"] == "ADAMTS15"]

for metric in ["n_scarlet_core","mean_scarlet_core","total_scarlet_core"]:
    a4_vals  = a4[metric].dropna()
    a15_vals = a15[metric].dropna()
    if len(a4_vals) < 2 or len(a15_vals) < 2:
        continue
    t, p = stats.ttest_ind(a4_vals, a15_vals)
    print(f"\n  {metric}:")
    print(f"    ADAMTS4:  mean={a4_vals.mean():.1f} ± {a4_vals.std():.1f}")
    print(f"    ADAMTS15: mean={a15_vals.mean():.1f} ± {a15_vals.std():.1f}")
    print(f"    t={t:.3f}, p={p:.4f}")

# =============================================================================
# FIGURES
# =============================================================================
print("\nGenerating figures...")

fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# Fig 1: Core mScarlet cell count by treatment
ax = axes[0, 0]
treat_order = [t for t in TREATMENT_ORDER if t in analysis_df["treatment"].values]
means = analysis_df.groupby("treatment")["n_scarlet_core"].mean()[treat_order]
sems  = analysis_df.groupby("treatment")["n_scarlet_core"].sem()[treat_order]
bars  = ax.bar(treat_order, means,
               color=[PALETTE[t] for t in treat_order], alpha=0.85)
ax.errorbar(range(len(treat_order)), means, yerr=sems,
            fmt="none", color="black", capsize=4)
for i, (treat, row) in enumerate(analysis_df.groupby("treatment")):
    vals = row["n_scarlet_core"].dropna()
    ax.scatter([i]*len(vals), vals, color="white", edgecolors="black",
               s=40, zorder=3)
ax.set_title("Core mScarlet cell count\n(proxy for transgene expression)")
ax.set_ylabel("Cells in Core zone")
ax.set_xticklabels(treat_order, rotation=35, ha="right", fontsize=9)

# Fig 2: Expression vs enwrapment scatter (all groups)
ax = axes[0, 1]
for treat in TREATMENT_ORDER:
    sub = analysis_df[analysis_df["treatment"] == treat].dropna(
        subset=["n_scarlet_core","frac_enwrapped"])
    ax.scatter(sub["n_scarlet_core"], sub["frac_enwrapped"],
               color=PALETTE.get(treat,"grey"), label=treat, s=60, alpha=0.85)
r_all, p_all = stats.pearsonr(
    analysis_df["n_scarlet_core"].dropna(),
    analysis_df.loc[analysis_df["n_scarlet_core"].notna(), "frac_enwrapped"])
ax.set_xlabel("Core mScarlet cells")
ax.set_ylabel("Fraction enwrapped")
ax.set_title(f"Expression vs enwrapment (all groups)\nr={r_all:.3f}, p={p_all:.4f}")
ax.legend(fontsize=7, loc="upper right")

# Fig 3: Expression vs enwrapment — active treatments only
ax = axes[0, 2]
active = ["ADAMTS15","C6ST1_ADAMTS15","ADAMTS4","ADAMTS4_MD"]
sub_act = analysis_df[analysis_df["treatment"].isin(active)].dropna(
    subset=["n_scarlet_core","frac_enwrapped"])
for treat in active:
    s = sub_act[sub_act["treatment"] == treat]
    ax.scatter(s["n_scarlet_core"], s["frac_enwrapped"],
               color=PALETTE.get(treat,"grey"), label=treat, s=60, alpha=0.85)
if len(sub_act) >= 3:
    r_act, p_act = stats.pearsonr(sub_act["n_scarlet_core"],
                                   sub_act["frac_enwrapped"])
    x_line = np.linspace(sub_act["n_scarlet_core"].min(),
                         sub_act["n_scarlet_core"].max(), 50)
    slope, intercept, *_ = stats.linregress(sub_act["n_scarlet_core"],
                                             sub_act["frac_enwrapped"])
    ax.plot(x_line, slope*x_line + intercept, "k--", lw=1.2, alpha=0.6)
    ax.set_title(f"Expression vs enwrapment (active treatments)\nr={r_act:.3f}, p={p_act:.4f}")
ax.set_xlabel("Core mScarlet cells")
ax.set_ylabel("Fraction enwrapped")
ax.legend(fontsize=8)

# Fig 4: Core/Outside ratio by treatment
ax = axes[1, 0]
if "core_outside_ratio" in zone_pivot.columns:
    ratio_by_treat = zone_pivot.groupby("treatment")["core_outside_ratio"].mean()
    ratio_se       = zone_pivot.groupby("treatment")["core_outside_ratio"].sem()
    t_list = [t for t in TREATMENT_ORDER if t in ratio_by_treat.index]
    ax.bar(t_list, ratio_by_treat[t_list],
           color=[PALETTE[t] for t in t_list], alpha=0.85)
    ax.errorbar(range(len(t_list)), ratio_by_treat[t_list],
                yerr=ratio_se[t_list], fmt="none", color="black", capsize=4)
    ax.axhline(1.0, color="grey", linestyle="--", lw=0.8)
    ax.set_title("Core/Outside enwrapment ratio\n(< 1 = Core more reduced)")
    ax.set_ylabel("Ratio")
    ax.set_xticklabels(t_list, rotation=35, ha="right", fontsize=9)

# Fig 5: Zone gradient — treated vs mScarlet
ax = axes[1, 1]
zone_line = zone_df.groupby(["treatment","zone"])["frac_enwrapped"].mean().reset_index()
zone_order_plot = ["Core","Penumbra","Outside"]
for treat in TREATMENT_ORDER:
    sub = zone_line[zone_line["treatment"] == treat]
    sub = sub.set_index("zone").reindex(zone_order_plot)
    lw  = 2.5 if treat == "mScarlet" else 1.2
    ls  = "--" if treat == "mScarlet" else "-"
    ax.plot(zone_order_plot, sub["frac_enwrapped"].values,
            color=PALETTE.get(treat,"grey"), lw=lw, ls=ls,
            marker="o", markersize=5, label=treat)
ax.set_title("Enwrapment zone gradient by treatment")
ax.set_ylabel("Mean fraction enwrapped")
ax.legend(fontsize=7)

# Fig 6: ADAMTS4 vs ADAMTS15 expression
ax = axes[1, 2]
for treat in ["ADAMTS4","ADAMTS15"]:
    sub = analysis_df[analysis_df["treatment"] == treat]["n_scarlet_core"].dropna()
    ax.bar([treat], [sub.mean()], color=PALETTE[treat], alpha=0.85)
    ax.errorbar([treat], [sub.mean()], yerr=[sub.sem()],
                fmt="none", color="black", capsize=4)
    ax.scatter([treat]*len(sub), sub, color="white",
               edgecolors="black", s=50, zorder=3)
ax.set_title("ADAMTS4 vs ADAMTS15\nCore mScarlet cells")
ax.set_ylabel("Core mScarlet cells")

plt.tight_layout()
plt.savefig(OUT_DIR / "fig_scarlet_expression_analysis.pdf", bbox_inches="tight")
plt.savefig(OUT_DIR / "fig_scarlet_expression_analysis.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved fig_scarlet_expression_analysis.pdf")

print(f"\nScript 14 complete. Outputs: {OUT_DIR}")