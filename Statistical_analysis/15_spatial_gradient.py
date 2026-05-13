#!/usr/bin/env python3
"""
15_spatial_gradient.py
=======================
Models enwrapment as a continuous function of distance from the injection
centre, fitting a logistic distance-decay model per treatment group.

The injection centre per slice is estimated as the centroid of all Core-zone
PNN cells in that slice. Distance is computed from each PV cell to that
centroid. PV cells without a valid Core centroid in their slice are excluded.

Model (fitted per treatment group using GLMM):
    logit(P(enwrapped)) ~ distance_mm + (1|mouse_id)

The slope (beta_distance) is the spatial decay constant:
    - More negative = enwrapment falls off faster with distance from injection
    - mScarlet slope provides the injection-microenvironment baseline
    - Treatment slopes steeper than mScarlet = protease-driven spatial effect

Additionally fits a continuous dose-response curve within the Core zone only:
    logit(P(enwrapped)) ~ distance_mm + (1|mouse_id)  [Core only]

Outputs:
    decay_constants.csv         — slope per treatment (primary result)
    decay_constants_core.csv    — slope within Core zone only
    per_slice_gradient.csv      — slice-level summaries
    animal_gradient.csv         — animal-level decay estimates
    fig_spatial_gradient.pdf    — main figure

EXCLUSIONS / CORRECTIONS: standard pipeline corrections applied.
"""

import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.spatial import cKDTree
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
CELLS_CSV  = Path("/path/to/analysis_results/cells_with_zones.csv")
OUT_DIR    = Path("/path/to/project/"
                  "Analysis/Analyse_2/results/15_spatial_gradient")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Constants ─────────────────────────────────────────────────────────────────
TREATMENT_ORDER = ["mScarlet","ADAMTS4","ADAMTS4_MD",
                   "ADAMTS15","C6ST1","C6ST1_ADAMTS15"]
EXCLUDE         = ""      # any mouse_id containing this string will be excluded (e.g. injection failure)
COLOC_THRESH    = 30      # pixels — same as pipeline
PX_TO_MM        = 0.325   # microns per pixel at hires → mm
MIN_PV_SLICE    = 5       # minimum PV cells per slice to include
MIN_CORE_PNN    = 3       # minimum Core PNN cells for centroid estimate

PALETTE = {
    "mScarlet":        "#888888",
    "ADAMTS4":         "#4e9af1",
    "ADAMTS4_MD":      "#f17c4e",
    "ADAMTS15":        "#4ef196",
    "C6ST1":           "#c44ef1",
    "C6ST1_ADAMTS15":  "#f1c44e",
}

# ── Load and correct data ─────────────────────────────────────────────────────
print("Loading cells_with_zones...")
cells = pd.read_csv(CELLS_CSV, low_memory=False)

cells = cells[cells["mouse_id"] != EXCLUDE].copy()

# Restrict to visual cortex, ipsilateral hemisphere
is_vis  = cells["brain_area"].str.lower().str.contains("visual|vis", na=False)
cells   = cells[is_vis & (cells["hemisphere"] == "left")].copy()

print(f"  {len(cells)} cells retained (visual cortex, ipsilateral)")
print(f"  Treatments: {sorted(cells['treatment'].unique())}")

# ── Per-slice: compute injection centre and distances ─────────────────────────
print("\nComputing injection centres and distances...")

pv_rows  = cells[cells["cell_type"] == "PV"].copy()
pnn_rows = cells[cells["cell_type"] == "PNN"].copy()

# Injection centre = centroid of Core-zone PNN cells per slice
core_pnn = pnn_rows[pnn_rows["zone"] == "Core"]
centres  = (core_pnn
            .groupby(["mouse_id","treatment","slice_id"])
            .apply(lambda g: pd.Series({
                "cx": g["x_hires"].mean(),
                "cy": g["y_hires"].mean(),
                "n_core_pnn": len(g)
            }))
            .reset_index())
centres  = centres[centres["n_core_pnn"] >= MIN_CORE_PNN]
print(f"  Slices with valid Core centroid: {len(centres)}")

# ── Compute enwrapment and distance per PV cell ───────────────────────────────
print("Computing enwrapment and distances per PV cell...")
records = []

for (mid, treat, sl), pv_grp in pv_rows.groupby(
        ["mouse_id","treatment","slice_id"]):

    if len(pv_grp) < MIN_PV_SLICE:
        continue

    # Get injection centre for this slice
    c_row = centres[
        (centres["mouse_id"] == mid) &
        (centres["slice_id"] == sl)
    ]
    if len(c_row) == 0:
        continue
    cx, cy = c_row.iloc[0]["cx"], c_row.iloc[0]["cy"]

    # Get PNN coordinates for this slice
    pnn_grp = pnn_rows[
        (pnn_rows["mouse_id"] == mid) &
        (pnn_rows["slice_id"] == sl)
    ][["x_hires","y_hires"]].dropna()

    if len(pnn_grp) < MIN_PV_SLICE:
        continue

    # Enwrapment: nearest PNN within threshold
    tree = cKDTree(pnn_grp.values)
    pv_xy = pv_grp[["x_hires","y_hires"]].dropna()
    if len(pv_xy) == 0:
        continue
    dists_pnn, _ = tree.query(pv_xy.values, k=1)
    enwrapped = (dists_pnn <= COLOC_THRESH).astype(int)

    # Distance from each PV cell to injection centre (in mm)
    dx = (pv_xy["x_hires"].values - cx) * PX_TO_MM / 1000
    dy = (pv_xy["y_hires"].values - cy) * PX_TO_MM / 1000
    dist_mm = np.sqrt(dx**2 + dy**2)

    for i in range(len(pv_xy)):
        records.append({
            "mouse_id":   mid,
            "treatment":  treat,
            "slice_id":   sl,
            "zone":       pv_grp.iloc[i]["zone"],
            "enwrapped":  enwrapped[i],
            "dist_mm":    dist_mm[i],
            "cx": cx, "cy": cy
        })

cell_df = pd.DataFrame(records)
cell_df["treatment"] = pd.Categorical(
    cell_df["treatment"], categories=TREATMENT_ORDER, ordered=True)
cell_df.to_csv(OUT_DIR / "cell_level_distances.csv", index=False)
print(f"  {len(cell_df)} PV cells with distance computed")

# ── Bin analysis: enwrapment vs distance bins ─────────────────────────────────
print("\nComputing binned gradient...")

N_BINS   = 12
MAX_DIST = cell_df["dist_mm"].quantile(0.95)  # exclude extreme outliers
cell_trimmed = cell_df[cell_df["dist_mm"] <= MAX_DIST].copy()

cell_trimmed["dist_bin"] = pd.cut(
    cell_trimmed["dist_mm"],
    bins=N_BINS,
    labels=False
)
cell_trimmed["dist_bin_mid"] = cell_trimmed["dist_mm"].apply(
    lambda x: (np.floor(x / (MAX_DIST / N_BINS)) + 0.5) * (MAX_DIST / N_BINS)
)

binned = (cell_trimmed
          .groupby(["treatment","dist_bin_mid"])
          .agg(
              frac_enwrapped = ("enwrapped","mean"),
              n_cells        = ("enwrapped","count")
          )
          .reset_index())
binned.to_csv(OUT_DIR / "binned_gradient.csv", index=False)

# ── Fit logistic decay per treatment ─────────────────────────────────────────
print("\nFitting logistic decay curves...")

try:
    import statsmodels.formula.api as smf
    USE_GLMM = True
except ImportError:
    USE_GLMM = False
    print("  statsmodels not available — using scipy logistic regression")

decay_results = []

for treat in TREATMENT_ORDER:
    sub = cell_trimmed[cell_trimmed["treatment"] == treat].copy()
    if len(sub) < 50:
        continue

    if USE_GLMM:
        try:
            m = smf.mixedlm(
                "enwrapped ~ dist_mm",
                data=sub,
                groups=sub["mouse_id"]
            ).fit(reml=True, method="bfgs")
            beta_dist = m.fe_params["dist_mm"]
            se_dist   = m.bse["dist_mm"]
            p_dist    = m.pvalues["dist_mm"]
            beta_int  = m.fe_params["Intercept"]
            method    = "GLMM_linear"
        except Exception as e:
            print(f"  GLMM failed for {treat}: {e} — falling back to OLS")
            USE_GLMM = False

    if not USE_GLMM:
        # OLS fallback: regress frac_enwrapped on distance at slice level
        sl_agg = sub.groupby(["mouse_id","slice_id","dist_bin_mid"])[
            "enwrapped"].mean().reset_index()
        sl_agg.columns = ["mouse_id","slice_id","dist_mm","frac"]
        if len(sl_agg) < 4:
            continue
        slope, intercept, r, p_dist, se = stats.linregress(
            sl_agg["dist_mm"], sl_agg["frac"])
        beta_dist = slope
        se_dist   = se
        beta_int  = intercept
        p_dist    = p_dist
        method    = "OLS_slice"

    # Decay constant: half-distance (where enwrapment drops to 50% of intercept)
    # For linear model: dist_50 = -0.5 * beta_int / beta_dist
    if beta_dist != 0:
        dist_50 = -0.5 * beta_int / beta_dist
    else:
        dist_50 = np.nan

    decay_results.append({
        "treatment":   treat,
        "beta_dist":   round(beta_dist, 5),
        "se_dist":     round(se_dist, 5),
        "p_dist":      round(p_dist, 4),
        "beta_int":    round(beta_int, 4),
        "dist_50_mm":  round(dist_50, 3),
        "n_cells":     len(sub),
        "n_animals":   sub["mouse_id"].nunique(),
        "method":      method
    })
    print(f"  {treat}: β={beta_dist:.5f} (p={p_dist:.4f}), "
          f"dist_50={dist_50:.2f}mm, n={len(sub)}")

decay_df = pd.DataFrame(decay_results)
decay_df.to_csv(OUT_DIR / "decay_constants.csv", index=False)

# ── Core-zone only: how fast does enwrapment recover from the injection site? ─
print("\nCore-zone only decay analysis...")
core_decay = []
for treat in TREATMENT_ORDER:
    sub = cell_trimmed[
        (cell_trimmed["treatment"] == treat) &
        (cell_trimmed["zone"] == "Core")
    ].copy()
    if len(sub) < 20:
        continue
    sl_agg = sub.groupby(["mouse_id","slice_id","dist_bin_mid"])[
        "enwrapped"].mean().reset_index()
    sl_agg.columns = ["mouse_id","slice_id","dist_mm","frac"]
    if len(sl_agg) < 4:
        continue
    slope, intercept, r, p, se = stats.linregress(
        sl_agg["dist_mm"], sl_agg["frac"])
    core_decay.append({
        "treatment": treat,
        "beta_dist_core": round(slope, 5),
        "p_core": round(p, 4),
        "n_cells": len(sub)
    })
    print(f"  {treat} [Core]: β={slope:.5f} (p={p:.4f}), n={len(sub)}")

core_decay_df = pd.DataFrame(core_decay)
core_decay_df.to_csv(OUT_DIR / "decay_constants_core.csv", index=False)

# ── Per-animal decay ──────────────────────────────────────────────────────────
print("\nPer-animal decay estimates...")
animal_decay = []
for (mid, treat), sub in cell_trimmed.groupby(["mouse_id","treatment"]):
    if len(sub) < 20:
        continue
    sl_agg = sub.groupby(["slice_id","dist_bin_mid"])[
        "enwrapped"].mean().reset_index()
    sl_agg.columns = ["slice_id","dist_mm","frac"]
    if len(sl_agg) < 4:
        continue
    slope, intercept, r, p, se = stats.linregress(
        sl_agg["dist_mm"], sl_agg["frac"])
    animal_decay.append({
        "mouse_id": mid, "treatment": treat,
        "beta_dist": round(slope, 5),
        "p": round(p, 4), "r": round(r, 3),
        "n_cells": len(sub)
    })

animal_decay_df = pd.DataFrame(animal_decay)
animal_decay_df.to_csv(OUT_DIR / "animal_gradient.csv", index=False)

# ── Compare mScarlet slope to treated slopes ──────────────────────────────────
print("\nSlope comparison vs mScarlet:")
ms_slope = decay_df.loc[
    decay_df["treatment"] == "mScarlet", "beta_dist"].values
if len(ms_slope) > 0:
    ms_slope = ms_slope[0]
    for _, row in decay_df[decay_df["treatment"] != "mScarlet"].iterrows():
        diff = row["beta_dist"] - ms_slope
        print(f"  {row['treatment']}: β={row['beta_dist']:.5f} "
              f"(vs mScarlet {ms_slope:.5f}, diff={diff:+.5f})")

# Per-animal t-tests of slopes vs mScarlet
ms_slopes = animal_decay_df[
    animal_decay_df["treatment"] == "mScarlet"]["beta_dist"].values

print("\nOne-sample t-test: each treatment slope vs mScarlet mean:")
for treat in TREATMENT_ORDER[1:]:
    t_slopes = animal_decay_df[
        animal_decay_df["treatment"] == treat]["beta_dist"].values
    if len(t_slopes) < 2 or len(ms_slopes) < 2:
        continue
    t, p = stats.ttest_ind(t_slopes, ms_slopes)
    print(f"  {treat}: slopes={np.round(t_slopes,5)}, "
          f"t={t:.3f}, p={p:.4f}")

# ── Figures ───────────────────────────────────────────────────────────────────
print("\nGenerating figures...")

fig = plt.figure(figsize=(16, 12))
gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.35)

# ── Fig 1: Binned gradient curves, all treatments ─────────────────────────────
ax1 = fig.add_subplot(gs[0, :2])
for treat in TREATMENT_ORDER:
    sub = binned[binned["treatment"] == treat].sort_values("dist_bin_mid")
    if len(sub) == 0:
        continue
    lw = 2.5 if treat in ["mScarlet","ADAMTS15","C6ST1_ADAMTS15"] else 1.2
    ls = "--" if treat == "mScarlet" else "-"
    ax1.plot(sub["dist_bin_mid"], sub["frac_enwrapped"],
             color=PALETTE[treat], lw=lw, ls=ls,
             marker="o", markersize=4, label=treat)

ax1.set_xlabel("Distance from injection centre (mm)", fontsize=11)
ax1.set_ylabel("Fraction of PV cells enwrapped", fontsize=11)
ax1.set_title("Enwrapment vs distance from injection centre\n"
              "(visual cortex, ipsilateral, all zones)", fontsize=11)
ax1.legend(fontsize=8, loc="upper right")
ax1.axhline(0, color="grey", lw=0.5, ls=":")

# ── Fig 2: Decay constants bar chart ─────────────────────────────────────────
ax2 = fig.add_subplot(gs[0, 2])
treats_plot = [t for t in TREATMENT_ORDER if t in decay_df["treatment"].values]
betas = decay_df.set_index("treatment").loc[treats_plot, "beta_dist"]
ses   = decay_df.set_index("treatment").loc[treats_plot, "se_dist"]
ps    = decay_df.set_index("treatment").loc[treats_plot, "p_dist"]

bars = ax2.bar(range(len(treats_plot)),
               betas.values,
               color=[PALETTE[t] for t in treats_plot],
               alpha=0.85)
ax2.errorbar(range(len(treats_plot)), betas.values,
             yerr=ses.values, fmt="none", color="black", capsize=4)
for i, (t, p_val) in enumerate(zip(treats_plot, ps.values)):
    sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else \
          "*" if p_val < 0.05 else "." if p_val < 0.1 else ""
    if sig:
        ax2.text(i, betas.values[i] + ses.values[i] + 0.001,
                 sig, ha="center", fontsize=10, fontweight="bold")

ax2.axhline(0, color="grey", lw=0.5, ls="--")
ax2.set_xticks(range(len(treats_plot)))
ax2.set_xticklabels(treats_plot, rotation=40, ha="right", fontsize=8)
ax2.set_ylabel("Decay slope β (enwrap / mm)", fontsize=10)
ax2.set_title("Spatial decay constants\n(more negative = faster decay)", fontsize=10)

# ── Fig 3: Per-animal slopes ──────────────────────────────────────────────────
ax3 = fig.add_subplot(gs[1, 0])
for i, treat in enumerate(TREATMENT_ORDER):
    sub = animal_decay_df[animal_decay_df["treatment"] == treat]
    ax3.scatter([i]*len(sub), sub["beta_dist"],
                color=PALETTE[treat], s=60, alpha=0.85, zorder=3)
    if len(sub) >= 2:
        ax3.errorbar(i, sub["beta_dist"].mean(),
                     yerr=sub["beta_dist"].sem(),
                     fmt="o", color=PALETTE[treat],
                     markersize=10, markeredgecolor="black",
                     capsize=4, lw=2, zorder=4)

ax3.axhline(0, color="grey", lw=0.5, ls="--")
ax3.set_xticks(range(len(TREATMENT_ORDER)))
ax3.set_xticklabels(TREATMENT_ORDER, rotation=40, ha="right", fontsize=8)
ax3.set_ylabel("Per-animal slope β", fontsize=10)
ax3.set_title("Per-animal decay slopes\n(dots = animals, large = mean±SEM)", fontsize=10)

# ── Fig 4: Core-zone gradient ADAMTS15 vs mScarlet ───────────────────────────
ax4 = fig.add_subplot(gs[1, 1])
for treat in ["mScarlet","ADAMTS15","C6ST1_ADAMTS15"]:
    sub = cell_trimmed[
        (cell_trimmed["treatment"] == treat) &
        (cell_trimmed["zone"] == "Core")
    ]
    binned_core = sub.groupby(
        pd.cut(sub["dist_mm"], bins=8))["enwrapped"].agg(
        ["mean","count"]).reset_index()
    binned_core["dist_mid"] = binned_core["dist_mm"].apply(
        lambda x: x.mid if hasattr(x, "mid") else np.nan)
    binned_core = binned_core.dropna(subset=["dist_mid"])
    lw = 2.0 if treat in ["ADAMTS15","C6ST1_ADAMTS15"] else 1.2
    ls = "--" if treat == "mScarlet" else "-"
    ax4.plot(binned_core["dist_mid"], binned_core["mean"],
             color=PALETTE[treat], lw=lw, ls=ls,
             marker="o", markersize=4, label=treat)

ax4.set_xlabel("Distance from centre (mm)", fontsize=10)
ax4.set_ylabel("Fraction enwrapped (Core only)", fontsize=10)
ax4.set_title("Core-zone gradient:\nactive treatments vs mScarlet", fontsize=10)
ax4.legend(fontsize=8)

# ── Fig 5: Zone × distance interaction ───────────────────────────────────────
ax5 = fig.add_subplot(gs[1, 2])
for zone, ls in [("Core","--"), ("Penumbra","-."), ("Outside","-")]:
    for treat in ["mScarlet","ADAMTS15","C6ST1_ADAMTS15"]:
        sub = cell_trimmed[
            (cell_trimmed["treatment"] == treat) &
            (cell_trimmed["zone"] == zone)
        ]
        if len(sub) < 20:
            continue
        zone_binned = sub.groupby(
            pd.cut(sub["dist_mm"], bins=6))["enwrapped"].mean().reset_index()
        zone_binned["dist_mid"] = zone_binned["dist_mm"].apply(
            lambda x: x.mid if hasattr(x, "mid") else np.nan)
        zone_binned = zone_binned.dropna(subset=["dist_mid"])
        alpha = 0.9 if treat != "mScarlet" else 0.4
        ax5.plot(zone_binned["dist_mid"], zone_binned["enwrapped"],
                 color=PALETTE[treat], lw=1.5, ls=ls, alpha=alpha)

# Legend
legend_elements = [
    Line2D([0],[0], color=PALETTE["mScarlet"],    lw=2, label="mScarlet"),
    Line2D([0],[0], color=PALETTE["ADAMTS15"],     lw=2, label="ADAMTS15"),
    Line2D([0],[0], color=PALETTE["C6ST1_ADAMTS15"], lw=2, label="C6ST1_ADAMTS15"),
    Line2D([0],[0], color="grey", lw=2, ls="--",  label="Core"),
    Line2D([0],[0], color="grey", lw=2, ls="-.",  label="Penumbra"),
    Line2D([0],[0], color="grey", lw=2, ls="-",   label="Outside"),
]
ax5.legend(handles=legend_elements, fontsize=6, loc="upper left")
ax5.set_xlabel("Distance from centre (mm)", fontsize=10)
ax5.set_ylabel("Fraction enwrapped", fontsize=10)
ax5.set_title("Distance decay by zone\n(active treatments vs mScarlet)", fontsize=10)

plt.savefig(OUT_DIR / "fig_spatial_gradient.pdf", bbox_inches="tight")
plt.savefig(OUT_DIR / "fig_spatial_gradient.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved fig_spatial_gradient.pdf")

print(f"\nScript 15 complete. Outputs: {OUT_DIR}")
