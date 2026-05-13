#!/usr/bin/env python3
"""
unified_protease_pipeline.py
=============================
Pre-LMM statistical pipeline for PNN enwrapment and fluorescence data.
Runs five complementary analyses before the primary LMM pipeline.

Sections:
  S1 — Virus+ vs Virus− cell intensity (one-sample t-tests per group)
  S2 — PV–PNN colocalization (Welch's t-tests vs reference)
  S3 — Ipsilateral vs Contralateral hemisphere (paired within-animal tests)
  S4 — Core/Outside zone gradient (Welch's t-tests vs reference)
  S5 — Power analysis (post-hoc power from observed effect sizes)

Also includes a metric correlation diagnostic (Cell 27) at the end.

All BH-FDR corrections are applied within test families (stain × zone or
zone × metric), not globally.

Usage:
  python unified_protease_pipeline.py

Requirements:
  pip install numpy pandas scipy statsmodels matplotlib seaborn
"""

import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns

# =============================================================================
# SECTION 0 — SHARED SETUP
# =============================================================================

# ── Paths ─────────────────────────────────────────────────────────────────────
ANALYSIS_DIR  = Path("/path/to/analysis_results")
ZONES_CSV     = Path("/path/to/cells_with_zones.csv")
MERGED_DIR    = Path("/path/to/merged_datasets")
OUTPUT_DIR    = Path("/path/to/results/unified_pipeline")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Experimental design ───────────────────────────────────────────────────────
TREATMENTS  = ["ADAMTS4", "ADAMTS4_MD", "ADAMTS15", "C6ST1", "C6ST1_ADAMTS15", "mScarlet"]
TREAT_ORDER = ["mScarlet", "ADAMTS4", "ADAMTS4_MD", "ADAMTS15", "C6ST1", "C6ST1_ADAMTS15"]
ALPHA       = 0.05

# ── Confirmed label corrections ───────────────────────────────────────────────
# If animal label swaps were confirmed diagnostically, list them here as
# {animal_id: corrected_treatment}. Leave empty ({}) if no swaps apply.
SWAP = {}   # e.g. {"ANIMAL_X": "GROUP_A", "ANIMAL_Y": "GROUP_B"}

# Treatment name aliases — rename any legacy treatment labels to canonical names.
# Keys are old names, values are new names. Leave empty ({}) if not needed.
TREATMENT_RENAME = {}   # e.g. {"OLD_LABEL": "NEW_LABEL"}

# Hemisphere convention
IPSI_HEMI   = "left"
CONTRA_HEMI = "right"

# ── Helpers ───────────────────────────────────────────────────────────────────
def cohens_d(a, b):
    """Cohen's d for independent groups."""
    pooled_sd = np.sqrt((np.std(a, ddof=1)**2 + np.std(b, ddof=1)**2) / 2)
    return (np.mean(a) - np.mean(b)) / pooled_sd if pooled_sd > 0 else np.nan

def cohens_d_one_sample(diffs):
    """Cohen's d for one-sample test (diffs vs 0)."""
    sd = np.std(diffs, ddof=1)
    return np.mean(diffs) / sd if sd > 0 else np.nan

def apply_fdr(pvals):
    """BH-FDR correction. Returns adjusted p-values."""
    if len(pvals) == 0:
        return np.array([])
    _, adj, _, _ = multipletests(pvals, method="fdr_bh")
    return adj

def sig_label(p_adj):
    if p_adj < 0.001: return "***"
    if p_adj < 0.01:  return "**"
    if p_adj < 0.05:  return "*"
    return ""

print(f"Output directory: {OUTPUT_DIR}")

# =============================================================================
# SECTION 1 — VIRUS+ vs VIRUS− CELL INTENSITY
# =============================================================================
# Question: Do cells within the virus injection radius have lower WFA/PV
#   fluorescence intensity than cells outside the injection radius?
# Design: Per-animal difference (Virus+ mean − Virus− mean) in each
#   stain × zone. One-sample t-test vs 0 within each treatment group.
#   BH-FDR per stain × zone family.
# Input: cells_with_zones.csv + per-animal Virus CSV files

VIRUS_RADIUS = 30  # pixels
MIN_CELLS    = 5

df = pd.read_csv(ZONES_CSV, low_memory=False)

# Apply swap and rename corrections
df["treatment"] = df.apply(
    lambda r: SWAP.get(r["mouse_id"], r["treatment"]), axis=1)
if TREATMENT_RENAME:
    df["treatment"] = df["treatment"].replace(TREATMENT_RENAME)

print(f"\nLoaded {len(df):,} cells")

# Tag virus-positive cells via KDTree if not already present
if "virus_positive" not in df.columns:
    print("\nvirus_positive column not found — computing from Virus CSV files...")
    df["virus_positive"] = False

    for mouse_id, animal_df in df.groupby("mouse_id"):
        virus_path = ANALYSIS_DIR / f"{mouse_id}_cell_fluorescence_analysis_Virus.csv"
        if not virus_path.exists():
            print(f"  WARNING: no virus file for {mouse_id}, skipping")
            continue

        virus = pd.read_csv(virus_path, usecols=["slice_id", "x_hires", "y_hires"])

        for slice_id, slice_virus in virus.groupby("slice_id"):
            v_coords = slice_virus[["x_hires", "y_hires"]].dropna().values
            if len(v_coords) == 0:
                continue
            mask = (df["mouse_id"] == mouse_id) & (df["slice_id"] == slice_id)
            cell_coords = df.loc[mask, ["x_hires", "y_hires"]].values
            if len(cell_coords) == 0:
                continue
            tree = cKDTree(v_coords)
            dists, _ = tree.query(cell_coords, k=1)
            df.loc[mask, "virus_positive"] = dists <= VIRUS_RADIUS

    n_pos = df["virus_positive"].sum()
    print(f"Tagged {n_pos:,} virus-positive cells ({n_pos/len(df)*100:.1f}%)")
else:
    print("virus_positive column found — using existing values")

print(f"\nTreatment counts (post-correction):")
print(df.drop_duplicates("mouse_id")[["mouse_id", "treatment"]].sort_values("treatment")
        .to_string(index=False))

# Compute per-animal Virus+ vs Virus− intensity (ipsilateral hemisphere only)
records = []
for mouse_id, animal_df in df.groupby("mouse_id"):
    treatment = animal_df["treatment"].iloc[0]
    ipsi = animal_df[animal_df["hemisphere"] == IPSI_HEMI]

    for stain in ["PNN", "PV"]:
        stain_df = ipsi[ipsi["cell_type"] == stain]
        for zone in ["Core", "Penumbra"]:
            zone_df = stain_df[stain_df["zone"] == zone]
            vpos = zone_df[zone_df["virus_positive"] == True]["normalized_btm20"].dropna()
            vneg = zone_df[zone_df["virus_positive"] == False]["normalized_btm20"].dropna()
            if len(vpos) >= MIN_CELLS and len(vneg) >= MIN_CELLS:
                records.append({
                    "mouse_id": mouse_id,
                    "treatment": treatment,
                    "stain": stain,
                    "zone": zone,
                    "vpos_mean": vpos.mean(),
                    "vneg_mean": vneg.mean(),
                    "diff": vpos.mean() - vneg.mean(),
                    "n_vpos": len(vpos),
                    "n_vneg": len(vneg),
                })

intensity = pd.DataFrame(records)
print(f"\nAnimal × stain × zone records: {len(intensity)}")

# One-sample t-tests per treatment × stain × zone
all_tests, families = [], {}

for (stain, zone), grp in intensity.groupby(["stain", "zone"]):
    family_key = f"{stain}_{zone}"
    family_tests = []

    for treatment in [t for t in TREAT_ORDER if t != "mScarlet"]:
        tdata = grp[grp["treatment"] == treatment]["diff"].dropna()
        if len(tdata) < 3:
            continue
        t_stat, p_raw = stats.ttest_1samp(tdata, popmean=0)
        d = cohens_d_one_sample(tdata.values)
        family_tests.append({
            "stain": stain, "zone": zone, "treatment": treatment,
            "n": len(tdata), "mean_diff": tdata.mean(),
            "t": t_stat, "df": len(tdata) - 1,
            "p_raw": p_raw, "cohen_d": d, "family": family_key
        })
    families[family_key] = family_tests

for fam, tests in families.items():
    if not tests:
        continue
    pvals = [t["p_raw"] for t in tests]
    padj  = apply_fdr(pvals)
    for t, pa in zip(tests, padj):
        t["p_adj"] = pa
        t["sig"] = sig_label(pa)
        all_tests.append(t)

s1_df = pd.DataFrame(all_tests)
s1_df.to_csv(OUTPUT_DIR / "S1_virus_intensity.csv", index=False)

print("\n=== SECTION 1 RESULTS: Virus+ vs Virus− intensity ===")
print(s1_df[["stain", "zone", "treatment", "n", "mean_diff", "t", "df",
             "p_raw", "p_adj", "sig", "cohen_d"]]
      .sort_values("p_raw").to_string(index=False))

# Section 1 figure: Cohen's d heatmap
fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for ax, stain in zip(axes, ["PNN", "PV"]):
    pivot = s1_df[s1_df["stain"] == stain].pivot(
        index="treatment", columns="zone", values="cohen_d"
    ).reindex([t for t in TREAT_ORDER if t != "mScarlet"])
    sig_pivot = s1_df[s1_df["stain"] == stain].pivot(
        index="treatment", columns="zone", values="sig"
    ).reindex([t for t in TREAT_ORDER if t != "mScarlet"]).fillna("")

    vals = pivot.values.ravel()
    finite_vals = vals[np.isfinite(vals)]
    vmax = np.max(np.abs(finite_vals)) if len(finite_vals) > 0 else 1.0

    sns.heatmap(pivot, ax=ax, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
                annot=True, fmt=".2f", linewidths=0.5, annot_kws={"size": 9})

    for ri, row_label in enumerate(pivot.index):
        for ci, col_label in enumerate(pivot.columns):
            if col_label in sig_pivot.columns and row_label in sig_pivot.index:
                s = sig_pivot.loc[row_label, col_label]
                if s:
                    ax.text(ci + 0.5, ri + 0.12, s, ha="center", va="top",
                            fontsize=10, fontweight="bold")
    ax.set_title(f"S1: Virus+/- intensity — {stain}\n(Cohen's d; * FDR<0.05)")

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "S1_cohend_heatmap.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved S1_cohend_heatmap.png")

# =============================================================================
# SECTION 2 — PV–PNN COLOCALIZATION (ENWRAPMENT)
# =============================================================================
# Question: Do treatment groups show reduced PNN enwrapment of PV cells?
# Design: Per-animal colocalisation metrics per zone, then Welch's t-test
#   vs reference group. BH-FDR per zone × metric family.
# Metrics: frac_pv_with_pnn, frac_pnn_with_pv, n_coloc_pairs

COLOC_RADIUS = 30  # pixels

coloc_records = []
for mouse_id, animal_df in df.groupby("mouse_id"):
    treatment = animal_df["treatment"].iloc[0]

    for zone in ["Core", "Penumbra", "Outside", "Contralateral"]:
        zone_df = animal_df[animal_df["zone"] == zone]
        pv_z    = zone_df[zone_df["cell_type"] == "PV"][["x_hires", "y_hires"]].dropna()
        pnn_z   = zone_df[zone_df["cell_type"] == "PNN"][["x_hires", "y_hires"]].dropna()

        if len(pv_z) < 3 or len(pnn_z) < 3:
            continue

        pv_coords  = pv_z.values
        pnn_coords = pnn_z.values

        pnn_tree = cKDTree(pnn_coords)
        pv_tree  = cKDTree(pv_coords)

        pv_has_pnn  = pnn_tree.query_ball_point(pv_coords,  COLOC_RADIUS)
        pnn_has_pv  = pv_tree.query_ball_point(pnn_coords, COLOC_RADIUS)

        frac_pv_pnn = np.mean([len(x) > 0 for x in pv_has_pnn])
        frac_pnn_pv = np.mean([len(x) > 0 for x in pnn_has_pv])
        n_pairs     = sum(len(x) > 0 for x in pv_has_pnn)

        coloc_records.append({
            "mouse_id": mouse_id, "treatment": treatment, "zone": zone,
            "frac_pv_with_pnn": frac_pv_pnn,
            "frac_pnn_with_pv": frac_pnn_pv,
            "n_coloc_pairs": n_pairs,
            "n_pv": len(pv_z), "n_pnn": len(pnn_z),
        })

coloc       = pd.DataFrame(coloc_records)
animal_coloc = (coloc
    .groupby(["mouse_id", "treatment", "zone"])
    [["frac_pv_with_pnn", "frac_pnn_with_pv", "n_coloc_pairs"]]
    .mean()
    .reset_index())

print(f"\nAnimal × zone records: {len(animal_coloc)}")

# Welch's t-tests: each treatment vs reference
METRICS_S2 = ["frac_pv_with_pnn", "frac_pnn_with_pv", "n_coloc_pairs"]
families_s2 = {}

for zone in animal_coloc["zone"].unique():
    for metric in METRICS_S2:
        ctrl = animal_coloc[
            (animal_coloc["treatment"] == "mScarlet") &
            (animal_coloc["zone"] == zone)
        ][metric].dropna()
        if len(ctrl) < 3:
            continue
        fam_key   = f"{zone}_{metric}"
        fam_tests = []
        for treatment in [t for t in TREAT_ORDER if t != "mScarlet"]:
            treat = animal_coloc[
                (animal_coloc["treatment"] == treatment) &
                (animal_coloc["zone"] == zone)
            ][metric].dropna()
            if len(treat) < 3:
                continue
            t_stat, p_raw = stats.ttest_ind(treat, ctrl, equal_var=False)
            d = cohens_d(treat.values, ctrl.values)
            fam_tests.append({
                "zone": zone, "metric": metric, "treatment": treatment,
                "n_treat": len(treat), "n_ctrl": len(ctrl),
                "mean_treat": treat.mean(), "mean_ctrl": ctrl.mean(),
                "mean_diff": treat.mean() - ctrl.mean(),
                "t": t_stat, "df": len(treat) + len(ctrl) - 2,
                "p_raw": p_raw, "cohen_d": d, "family": fam_key,
            })
        families_s2[fam_key] = fam_tests

all_s2 = []
for fam, tests in families_s2.items():
    if not tests:
        continue
    pvals = [t["p_raw"] for t in tests]
    padj  = apply_fdr(pvals)
    for t, pa in zip(tests, padj):
        t["p_adj"] = pa
        t["sig"] = sig_label(pa)
        all_s2.append(t)

s2_df = pd.DataFrame(all_s2)
s2_df.to_csv(OUTPUT_DIR / "S2_colocalization.csv", index=False)

print("\n=== SECTION 2 RESULTS: PV-PNN Colocalization ===")
sig_s2 = s2_df[s2_df["sig"] != ""]
print(f"Significant hits: {len(sig_s2)}")
print(s2_df[["zone", "metric", "treatment", "mean_treat", "mean_ctrl",
             "t", "p_raw", "p_adj", "sig", "cohen_d"]]
      .sort_values("p_raw").head(20).to_string(index=False))

# Section 2 figure
metric_plot = "frac_pv_with_pnn"
zones_plot  = ["Core", "Penumbra", "Outside", "Contralateral"]
fig, axes   = plt.subplots(1, len(zones_plot), figsize=(14, 4), sharey=True)

for ax, zone in zip(axes, zones_plot):
    sub = (animal_coloc[animal_coloc["zone"] == zone]
           .groupby("treatment")[metric_plot]
           .agg(["mean", "std"])
           .reindex(TREAT_ORDER))
    colors = ["#555555"] + ["#2E86AB"] * (len(TREAT_ORDER) - 1)
    ax.bar(range(len(TREAT_ORDER)), sub["mean"], yerr=sub["std"],
           color=colors, capsize=4, error_kw={"linewidth": 1.2})
    for i, treat in enumerate(TREAT_ORDER):
        if treat == "mScarlet":
            continue
        hits = s2_df[
            (s2_df["zone"] == zone) &
            (s2_df["metric"] == metric_plot) &
            (s2_df["treatment"] == treat) &
            (s2_df["sig"] != "")
        ]
        if len(hits) and not np.isnan(sub["mean"].iloc[i]):
            ypos = (sub["mean"].iloc[i] +
                    (sub["std"].iloc[i] if not np.isnan(sub["std"].iloc[i]) else 0) + 0.01)
            ax.text(i, ypos, hits["sig"].values[0],
                    ha="center", fontsize=10, fontweight="bold")
    ax.set_title(f"{zone}", fontsize=9)
    ax.set_xticks(range(len(TREAT_ORDER)))
    ax.set_xticklabels(TREAT_ORDER, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("frac_pv_with_pnn" if zone == zones_plot[0] else "")

plt.suptitle("S2: PV-PNN colocalization by zone", fontsize=10)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "S2_colocalization.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved S2_colocalization.png")

# =============================================================================
# SECTION 3 — IPSILATERAL vs CONTRALATERAL HEMISPHERE
# =============================================================================
# Question: Does the injected hemisphere differ from the uninjected hemisphere?
# Design: Within-animal paired comparison on region-level metrics.
#   Each animal contributes one difference score (mean ipsi − mean contra
#   across all brain areas). One-sample t-test: H₀: mean(diff) = 0.
# Input: merged_dataset_*.csv (auto-detects most recent file in MERGED_DIR)

print("\n=== SECTION 3: Ipsilateral vs Contralateral ===")
merged_files = sorted(MERGED_DIR.glob("merged_dataset_*.csv"))
if not merged_files:
    print("  No merged_dataset file found in MERGED_DIR — skipping Section 3.")
    s3_df = pd.DataFrame()
else:
    merged_path = merged_files[-1]
    print(f"Loading: {merged_path.name}")
    merged = pd.read_csv(merged_path, low_memory=False)

    # Apply swap and rename corrections (merged uses 'animal_id' column)
    SWAP_MERGED = SWAP.copy()
    merged["treatment"] = merged.apply(
        lambda r: SWAP_MERGED.get(r["animal_id"], r["treatment"]), axis=1)
    if TREATMENT_RENAME:
        merged["treatment"] = merged["treatment"].replace(TREATMENT_RENAME)

    print("Treatments after correction:", sorted(merged["treatment"].unique()))

    METRICS_S3  = ["diffFluo", "avgPxIntensity"]
    families_s3 = {}

    for staining in merged["staining"].unique():
        for metric in METRICS_S3:
            if metric not in merged.columns:
                continue
            fam_key   = f"{staining}_{metric}"
            fam_tests = []

            for treatment in [t for t in TREAT_ORDER if t != "mScarlet"]:
                animal_diffs = []
                treat_animals = merged[
                    (merged["treatment"] == treatment) &
                    (merged["staining"] == staining)
                ]["animal_id"].unique()

                for animal_id in treat_animals:
                    animal_df   = merged[(merged["animal_id"] == animal_id) &
                                         (merged["staining"] == staining)]
                    ipsi_mean   = animal_df[animal_df["hemisphere"] == IPSI_HEMI][metric].mean()
                    contra_mean = animal_df[animal_df["hemisphere"] == CONTRA_HEMI][metric].mean()
                    if not np.isnan(ipsi_mean) and not np.isnan(contra_mean):
                        animal_diffs.append(ipsi_mean - contra_mean)

                if len(animal_diffs) < 3:
                    continue
                animal_diffs = np.array(animal_diffs)
                t_stat, p_raw = stats.ttest_1samp(animal_diffs, popmean=0)
                d = cohens_d_one_sample(animal_diffs)
                fam_tests.append({
                    "staining": staining, "metric": metric, "treatment": treatment,
                    "n_animals": len(animal_diffs),
                    "mean_diff": animal_diffs.mean(),
                    "sd_diff": animal_diffs.std(ddof=1),
                    "t": t_stat, "df": len(animal_diffs) - 1,
                    "p_raw": p_raw, "cohen_d": d, "family": fam_key,
                })
            families_s3[fam_key] = fam_tests

    all_s3 = []
    for fam, tests in families_s3.items():
        if not tests:
            continue
        pvals = [t["p_raw"] for t in tests]
        padj  = apply_fdr(pvals)
        for t, pa in zip(tests, padj):
            t["p_adj"] = pa
            t["sig"] = sig_label(pa)
            all_s3.append(t)

    s3_df = pd.DataFrame(all_s3) if all_s3 else pd.DataFrame()
    if not s3_df.empty:
        s3_df.to_csv(OUTPUT_DIR / "S3_hemisphere.csv", index=False)
        print("\n=== SECTION 3 RESULTS ===")
        print(s3_df[["staining", "metric", "treatment", "n_animals", "mean_diff",
                      "sd_diff", "t", "p_raw", "p_adj", "sig", "cohen_d"]]
              .sort_values("p_raw").to_string(index=False))
    else:
        print("  No results — check MERGED_DIR path and column names.")

# =============================================================================
# SECTION 4 — CORE/OUTSIDE ZONE GRADIENT
# =============================================================================
# Question: Does enwrapment show a Core < Outside gradient in treatment groups?
# Design: Core/Outside ratio per animal from Section 2 colocalization data.
#   Welch's t-test vs reference group. BH-FDR across treatments.

metric_s4 = "frac_pv_with_pnn"
co_records = []

for mouse_id, animal_df in animal_coloc.groupby("mouse_id"):
    treatment   = animal_df["treatment"].iloc[0]
    core_val    = animal_df[animal_df["zone"] == "Core"][metric_s4].mean()
    outside_val = animal_df[animal_df["zone"] == "Outside"][metric_s4].mean()
    if not np.isnan(core_val) and not np.isnan(outside_val) and outside_val > 0:
        co_records.append({
            "mouse_id": mouse_id, "treatment": treatment,
            "core": core_val, "outside": outside_val,
            "ratio": core_val / outside_val,
            "diff":  core_val - outside_val,
        })

co_df = pd.DataFrame(co_records)
print("\nCore/Outside ratio by treatment (frac_pv_with_pnn):")
print(co_df.groupby("treatment")[["core", "outside", "ratio"]].mean().round(3))

ctrl_ratio = co_df[co_df["treatment"] == "mScarlet"]["ratio"].dropna()
s4_tests   = []

for treatment in [t for t in TREAT_ORDER if t != "mScarlet"]:
    treat_ratio = co_df[co_df["treatment"] == treatment]["ratio"].dropna()
    if len(treat_ratio) < 3 or len(ctrl_ratio) < 3:
        continue
    t_stat, p_raw = stats.ttest_ind(treat_ratio, ctrl_ratio, equal_var=False)
    d = cohens_d(treat_ratio.values, ctrl_ratio.values)
    s4_tests.append({
        "treatment": treatment,
        "n": len(treat_ratio),
        "mean_ratio_treat": treat_ratio.mean(),
        "mean_ratio_ctrl": ctrl_ratio.mean(),
        "ratio_diff": treat_ratio.mean() - ctrl_ratio.mean(),
        "t": t_stat, "df": len(treat_ratio) + len(ctrl_ratio) - 2,
        "p_raw": p_raw, "cohen_d": d,
    })

pvals_s4 = [t["p_raw"] for t in s4_tests]
padj_s4  = apply_fdr(pvals_s4)
for t, pa in zip(s4_tests, padj_s4):
    t["p_adj"] = pa
    t["sig"] = sig_label(pa)

s4_df = pd.DataFrame(s4_tests)
s4_df.to_csv(OUTPUT_DIR / "S4_zone_gradient.csv", index=False)

print("\n=== SECTION 4 RESULTS: Core/Outside Zone Gradient ===")
print(s4_df[["treatment", "n", "mean_ratio_treat", "mean_ratio_ctrl",
             "ratio_diff", "t", "p_raw", "p_adj", "sig", "cohen_d"]]
      .to_string(index=False))

# Section 4 figure
zones_ord = ["Core", "Penumbra", "Outside"]
fig, ax   = plt.subplots(figsize=(7, 5))

palette = {t: c for t, c in zip(
    TREAT_ORDER,
    ["#555555", "#4C9BE8", "#7FBFFF", "#E84C4C", "#E8A04C", "#A04CE8"]
)}

for treatment in TREAT_ORDER:
    treat_data = animal_coloc[animal_coloc["treatment"] == treatment]
    means = [treat_data[treat_data["zone"] == z]["frac_pv_with_pnn"].mean()
             for z in zones_ord]
    sems  = [treat_data[treat_data["zone"] == z]["frac_pv_with_pnn"].sem()
             for z in zones_ord]
    lw = 2.5 if treatment == "mScarlet" else 1.5
    ls = "--" if treatment == "mScarlet" else "-"
    ax.errorbar(range(3), means, yerr=sems, label=treatment,
                color=palette.get(treatment, "grey"), lw=lw, ls=ls,
                marker="o", capsize=4)

ax.set_xticks(range(3))
ax.set_xticklabels(zones_ord)
ax.set_ylabel("Fraction PV cells with adjacent PNN")
ax.set_title("S4: Zone gradient — PV-PNN enwrapment\n(mean ± SEM per zone)")
ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "S4_zone_gradient.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved S4_zone_gradient.png")

# =============================================================================
# SECTION 5 — POWER ANALYSIS
# =============================================================================
# Purpose: Quantify statistical power at the current n per group using
#   the effect sizes actually observed in Sections 1–4.
# Method: Post-hoc power using t-distribution approximation.

from scipy.stats import t as t_dist

def compute_power(n, d, alpha=0.05, two_tailed=True):
    """Approximate power for independent two-sample t-test."""
    df_pw = 2 * n - 2
    ncp   = d * np.sqrt(n / 2)
    t_crit = t_dist.ppf(1 - alpha / 2, df_pw) if two_tailed else t_dist.ppf(1 - alpha, df_pw)
    power  = 1 - t_dist.cdf(t_crit, df_pw, ncp) + t_dist.cdf(-t_crit, df_pw, ncp)
    return max(0.0, min(1.0, power))

def n_for_power(d, target_power=0.80, alpha=0.05):
    """Minimum n per group to achieve target power for given d."""
    for n in range(3, 200):
        if compute_power(n, abs(d), alpha) >= target_power:
            return n
    return ">200"

# Collect effect sizes from Sections 1–4
effect_sizes = []

if not s1_df.empty:
    for _, row in s1_df.iterrows():
        if not np.isnan(row["cohen_d"]):
            effect_sizes.append({
                "section": "S1: Virus+/- intensity",
                "label":   f"{row['stain']} {row['zone']} {row['treatment']}",
                "cohen_d": abs(row["cohen_d"]),
                "n": row["n"], "p_raw": row["p_raw"], "sig": row["sig"],
            })

if not s2_df.empty:
    for _, row in s2_df[s2_df["metric"] == "frac_pv_with_pnn"].iterrows():
        if not np.isnan(row["cohen_d"]):
            effect_sizes.append({
                "section": "S2: PV-PNN enwrapment",
                "label":   f"{row['zone']} {row['treatment']}",
                "cohen_d": abs(row["cohen_d"]),
                "n": row["n_treat"], "p_raw": row["p_raw"], "sig": row["sig"],
            })

if not s4_df.empty:
    for _, row in s4_df.iterrows():
        if not np.isnan(row["cohen_d"]):
            effect_sizes.append({
                "section": "S4: Zone gradient",
                "label":   row["treatment"],
                "cohen_d": abs(row["cohen_d"]),
                "n": row["n"], "p_raw": row["p_raw"], "sig": row["sig"],
            })

es_df = pd.DataFrame(effect_sizes)
if not es_df.empty:
    es_df["power_at_n4"]   = es_df["cohen_d"].apply(lambda d: compute_power(4, d))
    es_df["n_for_80pct"]   = es_df["cohen_d"].apply(lambda d: n_for_power(d, 0.80))
    es_df["n_for_90pct"]   = es_df["cohen_d"].apply(lambda d: n_for_power(d, 0.90))
    es_df.to_csv(OUTPUT_DIR / "S5_power_analysis.csv", index=False)

    print("\n=== SECTION 5: Power at current n for all observed effect sizes ===")
    print(es_df[["section", "label", "cohen_d", "power_at_n4", "n_for_80pct", "p_raw", "sig"]]
          .sort_values("power_at_n4").to_string(index=False))

    # Figure 1: power curves
    n_range = np.arange(3, 40)
    d_vals  = sorted(es_df["cohen_d"].dropna().unique())
    indices = np.round(np.linspace(0, len(d_vals) - 1, min(5, len(d_vals)))).astype(int)
    d_plot  = [d_vals[i] for i in indices]

    fig, ax = plt.subplots(figsize=(8, 5))
    cmap    = plt.cm.viridis(np.linspace(0.1, 0.9, len(d_plot)))
    for d, color in zip(d_plot, cmap):
        powers = [compute_power(n, d) for n in n_range]
        ax.plot(n_range, powers, color=color, lw=2, label=f"d = {d:.2f}")
    ax.axhline(0.80, color="red", ls="--", lw=1.5, label="80% power threshold")
    ax.axvline(4,    color="black", ls=":", lw=1.5, label="n=4 (current study)")
    ax.set_xlabel("n per group")
    ax.set_ylabel("Statistical power (1 - β)")
    ax.set_title("S5: Power curves for observed effect sizes\n(two-sided t-test, α=0.05)")
    ax.legend(fontsize=8, loc="lower right")
    ax.set_xlim(3, 40)
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(1.0))
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "S5_power_curves.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Figure 2: power per test bar chart
    fig, ax = plt.subplots(figsize=(10, 5))
    colors  = ["#2E86AB" if p >= 0.80 else "#E84C4C" for p in es_df["power_at_n4"]]
    ax.barh(range(len(es_df)), es_df["power_at_n4"], color=colors)
    ax.axvline(0.80, color="red", ls="--", lw=1.5, label="80% threshold")
    ax.set_yticks(range(len(es_df)))
    ax.set_yticklabels(
        [f"{r['section'][:6]} | {r['label']}" for _, r in es_df.iterrows()],
        fontsize=7)
    ax.set_xlabel("Power at current n")
    ax.set_title("S5: Statistical power per test\nBlue = adequately powered (≥80%),  Red = underpowered")
    ax.legend()
    ax.xaxis.set_major_formatter(mticker.PercentFormatter(1.0))
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "S5_power_per_test.png", dpi=150, bbox_inches="tight")
    plt.close()

    print(f"\nSaved S5 figures.")

    underpowered = es_df[es_df["power_at_n4"] < 0.80]
    print(f"\n{len(underpowered)} of {len(es_df)} tests are underpowered (<80%).")
    print("Tests with adequate power only for large effects (Cohen's d ≥ ~2.1).")
    print("LMMs using slice-level data with animal as random intercept recover")
    print("statistical power while preserving the animal as the biological unit.")

# =============================================================================
# METRIC CORRELATION DIAGNOSTIC
# =============================================================================
# Computes animal-level means for all candidate PNN integrity metrics and
# produces a correlation matrix + scatter grid.
#
# Candidate metrics:
#   frac_enwrapped   — from script 00
#   normalized_btm20 — cell-level WFA intensity
#   avgPxIntensity   — region-level mean WFA pixel intensity
#   diffFluo         — region-level total WFA diffuse fluorescence
#   density          — WFA cells per mm² (from zone_density)
#
# All metrics: WFA staining, ipsilateral hemisphere, visual cortex, animal level.

print("\n=== METRIC CORRELATION DIAGNOSTIC ===\n")

CORR_CELLS_CSV    = ZONES_CSV
CORR_MERGED_ZONES = next(MERGED_DIR.glob("merged_dataset_zones_*.csv"), None)
CORR_MERGED_MAIN  = next(MERGED_DIR.glob("merged_dataset_[!z]*.csv"), None)
CORR_OUT_DIR      = OUTPUT_DIR / "metric_correlations"
CORR_OUT_DIR.mkdir(parents=True, exist_ok=True)

COLOC_RADIUS_CORR = 30
EXCLUDE_ANIMAL    = None  # set to animal_id string to exclude, or None

VIS_KEYWORDS = ["visual", "VIS"]

PALETTE_CORR = {t: c for t, c in zip(TREAT_ORDER,
    ["#888888", "#4e9af1", "#f17c4e", "#4ef196", "#c44ef1", "#f1c44e"])}

def apply_corrections_corr(df_in, id_col="mouse_id"):
    """Apply swap and rename corrections."""
    df_out = df_in.copy()
    if id_col in df_out.columns and "treatment" in df_out.columns:
        df_out["treatment"] = df_out.apply(
            lambda r: SWAP.get(r[id_col], r["treatment"]), axis=1)
    if TREATMENT_RENAME and "treatment" in df_out.columns:
        df_out["treatment"] = df_out["treatment"].replace(TREATMENT_RENAME)
    return df_out

# ── Metric 1: frac_enwrapped from cells_with_zones ───────────────────────────
try:
    cells_corr = pd.read_csv(
        CORR_CELLS_CSV,
        usecols=["mouse_id", "slice_id", "cell_type", "hemisphere",
                 "x_hires", "y_hires", "brain_area", "treatment"]
    )
    cells_corr = apply_corrections_corr(cells_corr)
    if EXCLUDE_ANIMAL:
        cells_corr = cells_corr[cells_corr["mouse_id"] != EXCLUDE_ANIMAL]

    vis_mask = cells_corr["brain_area"].str.contains("|".join(VIS_KEYWORDS),
                                                      case=False, na=False)
    cells_corr = cells_corr[vis_mask]
    ipsi_cells = cells_corr[cells_corr["hemisphere"] == IPSI_HEMI]

    pv_c  = ipsi_cells[ipsi_cells["cell_type"] == "PV"]
    pnn_c = ipsi_cells[ipsi_cells["cell_type"] == "PNN"]

    enw_records = []
    for mouse_id, a_pv in pv_c.groupby("mouse_id"):
        a_pnn = pnn_c[pnn_c["mouse_id"] == mouse_id]
        if len(a_pnn) == 0:
            continue
        pv_xy  = a_pv[["x_hires",  "y_hires"]].values
        pnn_xy = a_pnn[["x_hires", "y_hires"]].values
        tree   = cKDTree(pnn_xy)
        dists, _ = tree.query(pv_xy, k=1)
        treatment = a_pv["treatment"].iloc[0]
        enw_records.append({
            "mouse_id": mouse_id, "treatment": treatment,
            "frac_enwrapped": (dists <= COLOC_RADIUS_CORR).mean()
        })

    enw_df = pd.DataFrame(enw_records)
    print(f"frac_enwrapped: {len(enw_df)} animals")
except Exception as e:
    print(f"  frac_enwrapped failed: {e}")
    enw_df = pd.DataFrame()

# ── Metric 2: normalised intensity from cells_with_zones ─────────────────────
try:
    int_cells = pd.read_csv(
        CORR_CELLS_CSV,
        usecols=["mouse_id", "cell_type", "hemisphere", "brain_area",
                 "treatment", "normalized_btm20"]
    )
    int_cells = apply_corrections_corr(int_cells)
    if EXCLUDE_ANIMAL:
        int_cells = int_cells[int_cells["mouse_id"] != EXCLUDE_ANIMAL]
    vis_mask2 = int_cells["brain_area"].str.contains("|".join(VIS_KEYWORDS),
                                                      case=False, na=False)
    int_cells = int_cells[vis_mask2 &
                          (int_cells["hemisphere"] == IPSI_HEMI) &
                          (int_cells["cell_type"] == "PNN")]
    int_df = (int_cells.groupby(["mouse_id", "treatment"])["normalized_btm20"]
              .mean().reset_index()
              .rename(columns={"normalized_btm20": "cell_wfa_intensity"}))
    print(f"cell_wfa_intensity: {len(int_df)} animals")
except Exception as e:
    print(f"  cell_wfa_intensity failed: {e}")
    int_df = pd.DataFrame()

# ── Metrics 3 & 4: avgPxIntensity and diffFluo from merged_main ──────────────
px_df, diff_df = pd.DataFrame(), pd.DataFrame()
if CORR_MERGED_MAIN and CORR_MERGED_MAIN.exists():
    try:
        mdf = pd.read_csv(CORR_MERGED_MAIN, low_memory=False)
        mdf = apply_corrections_corr(mdf, id_col="animal_id")
        if EXCLUDE_ANIMAL:
            mdf = mdf[mdf["animal_id"] != EXCLUDE_ANIMAL]
        mdf_wfa = mdf[
            (mdf["staining"].str.upper() == "WFA") &
            (mdf["hemisphere"] == IPSI_HEMI)
        ]
        if "avgPxIntensity" in mdf.columns:
            px_df = (mdf_wfa.groupby(["animal_id", "treatment"])["avgPxIntensity"]
                     .mean().reset_index()
                     .rename(columns={"animal_id": "mouse_id",
                                      "avgPxIntensity": "avg_px_intensity"}))
            print(f"avg_px_intensity: {len(px_df)} animals")
        if "diffFluo" in mdf.columns:
            diff_df = (mdf_wfa.groupby(["animal_id", "treatment"])["diffFluo"]
                       .mean().reset_index()
                       .rename(columns={"animal_id": "mouse_id",
                                        "diffFluo": "diff_fluo"}))
            print(f"diff_fluo: {len(diff_df)} animals")
    except Exception as e:
        print(f"  merged_main metrics failed: {e}")

# ── Merge all metrics ─────────────────────────────────────────────────────────
wide = enw_df.copy() if not enw_df.empty else pd.DataFrame()

for other, col in [(int_df,  "cell_wfa_intensity"),
                   (px_df,   "avg_px_intensity"),
                   (diff_df, "diff_fluo")]:
    if not other.empty and not wide.empty:
        wide = wide.merge(other[["mouse_id", col]], on="mouse_id", how="outer")
    elif not wide.empty:
        wide[col] = np.nan

metric_cols = [c for c in ["frac_enwrapped", "cell_wfa_intensity",
                             "avg_px_intensity", "diff_fluo"]
               if c in wide.columns]

if len(metric_cols) >= 2 and len(wide) > 3:
    wide.to_csv(CORR_OUT_DIR / "animal_metrics_wide.csv", index=False)

    # Correlation matrix
    corr_data  = wide[metric_cols].dropna()
    corr_pairs = []
    for i, m1 in enumerate(metric_cols):
        for j, m2 in enumerate(metric_cols):
            if j <= i:
                continue
            x = corr_data[m1].values
            y = corr_data[m2].values
            if len(x) < 3:
                continue
            r, p = stats.pearsonr(x, y)
            corr_pairs.append({"metric_1": m1, "metric_2": m2,
                                "pearson_r": r, "p": p,
                                "n": len(x)})

    corr_df = pd.DataFrame(corr_pairs)
    corr_df.to_csv(CORR_OUT_DIR / "metric_correlation_matrix.csv", index=False)
    print(f"\nMetric correlation matrix ({len(corr_data)} animals):")
    print(corr_df.to_string(index=False))

    # Heatmap
    r_matrix = corr_data.corr()
    fig, ax  = plt.subplots(figsize=(6, 5))
    sns.heatmap(r_matrix, annot=True, fmt=".2f", cmap="RdBu_r", center=0,
                vmin=-1, vmax=1, ax=ax, square=True)
    ax.set_title("Metric correlations (Pearson r)\nWFA, ipsilateral, visual cortex")
    plt.tight_layout()
    plt.savefig(CORR_OUT_DIR / "metric_correlation_heatmap.png",
                dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved metric_correlation_heatmap.png")
else:
    print("  Insufficient data for correlation matrix.")

# =============================================================================
# PIPELINE SUMMARY
# =============================================================================
print("\n=== PIPELINE SUMMARY ===")
print(f"  S1: Virus+/- intensity   → {OUTPUT_DIR / 'S1_virus_intensity.csv'}")
print(f"  S2: Colocalization        → {OUTPUT_DIR / 'S2_colocalization.csv'}")
print(f"  S3: Hemisphere comparison → {OUTPUT_DIR / 'S3_hemisphere.csv'}")
print(f"  S4: Zone gradient         → {OUTPUT_DIR / 'S4_zone_gradient.csv'}")
if not es_df.empty:
    print(f"  S5: Power analysis        → {OUTPUT_DIR / 'S5_power_analysis.csv'}")
print(f"  Correlations              → {CORR_OUT_DIR}")
print("\nAll pre-LMM analyses complete. Proceed to R pipeline (scripts 01–06).")
