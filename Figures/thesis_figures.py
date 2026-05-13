#!/usr/bin/env python3
"""
Figure Generation Script — All Thesis Sections

Representative microscopy figures (§3.1, §3.4, §3.5, appendix D):
Six-group V1 crop composite (2x3 panels)
Same crops, mScarlet-only greyscale (supplementary)
Whole-section overview composite, V1 fill + atlas + crop box
Primary findings — overview + crop pairs (Control / A15 / C6+A15)
C6ST-1 interaction pairs (Control / C6 / C6+A15)
Monocular deprivation pairs (Control / A4 / A4-MD)

Main figures:
Primary LMM forest plot (log-transformed)
Animal-level enwrapment ratios
Layer-stratified forest plots (4 panels)
Distance-decay enwrapment curves
Factorial interaction (model prediction)
Zone profile — ADAMTS-4 vs ADAMTS-4-MD
Virus-positive vs virus-negative enwrapment (Bayesian)
Distance-bin enwrapment (line plot)
Virus-status enwrapment comparison
 Viral load dose-response curves

Supplementary / appendix figures:
WFA cell density (hemisphere-level LMM)
Off-target specificity (Bayesian)
Viral load dose-response slopes
Hemisphere WFA and PV fluorescence (preliminary)
Zone-stratified virus+/virus- fluorescence (preliminary)
Power analysis — Cohen's d per comparison

The microscopy figures require additional inputs (multi-channel TIFFs, atlas
overlay PNGs, VisuAlign JSON files, cells_with_zones.csv). If any of these
are missing, the microscopy section is skipped and the data-plot figures
still run normally. Configure paths in the MICROSCOPY CONFIGURATION block.

Usage:
  python thesis_figures.py

Requirements:
  pip install pandas numpy matplotlib scipy tifffile pillow
"""

import os
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy import stats as sp_stats

warnings.filterwarnings("ignore", category=UserWarning)

# =============================================================================
# CONFIGURATION
# =============================================================================
DATA_PATH   = "/path/to/ALL_RESULTS_COMBINED.csv"
PRE_LMM_DIR = "/path/to/pre_lmm_results"
OUTPUT_DIR  = "/path/to/output/figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# FIGURE DIMENSIONS
# =============================================================================
FIG_W, FIG_H = 6.3, 5.0
FIG_WIDTH    = 6.3
FIG_HEIGHT   = 5

# =============================================================================
# COLORBLIND-ACCESSIBLE PALETTE
# =============================================================================
TX_COLORS = {
    "mScarlet":       "#888888",
    "ADAMTS4":        "#9ECAE1",
    "ADAMTS4_MD":     "#2171B5",
    "ADAMTS15":       "#FDAE6B",
    "C6ST1":          "#74C476",
    "C6ST1_ADAMTS15": "#CB181D",
}
treatment_colors = {
    "mScarlet":        "#888888",
    "ADAMTS4":         "#9ECAE1",
    "ADAMTS4_MD":      "#2171B5",
    "ADAMTS15":        "#FDAE6B",
    "C6ST1":           "#74C476",
    "C6ST1_ADAMTS15":  "#CB181D",
}

TX_ORDER = ["mScarlet", "ADAMTS4", "ADAMTS4_MD", "ADAMTS15", "C6ST1", "C6ST1_ADAMTS15"]
treatment_order = [
    "mScarlet", "ADAMTS4", "ADAMTS4_MD",
    "ADAMTS15", "C6ST1", "C6ST1_ADAMTS15",
]

treatment_markers = {
    "mScarlet":        "o",
    "ADAMTS4":         "s",
    "ADAMTS4_MD":      "D",
    "ADAMTS15":        "^",
    "C6ST1":           "v",
    "C6ST1_ADAMTS15":  "P",
}
TX_MARKERS = {
    "mScarlet": "o", "ADAMTS4": "s", "ADAMTS4_MD": "D",
    "ADAMTS15": "^", "C6ST1":   "v", "C6ST1_ADAMTS15": "P",
}
TX5 = [t for t in TX_ORDER if t != "mScarlet"]

# =============================================================================
# PLOT SETTINGS
# =============================================================================
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "DejaVu Sans"],
    "font.size":         11,
    "axes.linewidth":    1,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        300,
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
    "savefig.facecolor": "white",
})

# =============================================================================
# LOAD DATA
# =============================================================================
print("Loading data...")
D = pd.read_csv(DATA_PATH, low_memory=False)
D["source_dir"] = D["source_dir"].str.replace(r"^results/", "", regex=True)

for c in ["beta", "se", "SE", "ci_lo", "ci_hi", "p", "p_adj", "ratio",
          "post_mean", "post_sd", "ci_89_lo", "ci_89_hi",
          "frac_enwrapped", "n_cells", "dist_bin_mid", "pct_change",
          "beta_dist", "p_dist", "slope", "OR_per_unit",
          "post_mean_diff", "p_viruspos_lower",
          "post_mean_visual", "post_mean_offtarget", "p_visual_gt_offtarget",
          "F value", "Pr(>F)", "Chisq", "Pr(>Chisq)",
          "estimate", "t.ratio", "p.value"]:
    if c in D.columns:
        D[c] = pd.to_numeric(D[c], errors="coerce")

def q(sdir, sfile):
    return D[(D["source_dir"] == sdir) & (D["source_file"] == sfile)].copy()

ANIM = q("00_enwrapment", "animal_enwrapment.csv")[
    ["treatment", "mouse_id", "ratio"]].dropna()
MS_LOGMEAN = np.log(ANIM[ANIM["treatment"] == "mScarlet"]["ratio"].values).mean()

print(f"  Loaded {len(D)} rows")

# =============================================================================
# HELPERS
# =============================================================================
def stars(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return ""


def _filled(ax, x, y, mkr, color, sig, s=130):
    if sig:
        ax.scatter(x, y, s=s, marker=mkr, c=color,
                   edgecolors="black", linewidths=1.1, zorder=4)
    else:
        ax.scatter(x, y, s=s, marker=mkr, c="white",
                   edgecolors=color, linewidths=1.8, zorder=4)

# =============================================================================
# BRACKET HELPERS (Option B: single vertical rail)
# =============================================================================
BRACKET_LW = 1.5
BRACKET_FS = 13
BRACKET_TW = 0.04


def bracket_B_h(ax, y_pos, sig_dict, ref_y, rail_x,
                tick_w=BRACKET_TW, fontsize=BRACKET_FS, step=0.18):
    items = [(t, s) for t, s in sig_dict.items() if s and t in y_pos]
    if not items:
        return
    max_stars = max(len(s) for _, s in items)
    effective_step = step * (1 + 0.2 * (max_stars - 1))
    x_positions = [rail_x + idx * effective_step for idx in range(len(items))]

    for idx, (tx, st) in enumerate(items):
        t_y  = y_pos[tx]
        x_br = x_positions[idx]
        ax.plot([x_br, x_br], [ref_y, t_y],
                color="black", lw=BRACKET_LW, clip_on=False, zorder=5)
        ax.plot([x_br - tick_w, x_br], [ref_y, ref_y],
                color="black", lw=BRACKET_LW, clip_on=False, zorder=5)
        ax.plot([x_br - tick_w, x_br], [t_y, t_y],
                color="black", lw=BRACKET_LW, clip_on=False, zorder=5)
        mid_y = (ref_y + t_y) / 2
        if idx < len(x_positions) - 1:
            mid_x = (x_positions[idx] + x_positions[idx + 1]) / 2
        else:
            mid_x = x_positions[idx] + effective_step / 2
        ax.text(mid_x, mid_y, st,
                ha="center", va="center", fontsize=fontsize,
                fontweight="bold", clip_on=False, zorder=5)


def bracket_B_v(ax, x_pos, sig_dict, ref_x, rail_y,
                tick_h=BRACKET_TW, fontsize=BRACKET_FS):
    items = [(t, s) for t, s in sig_dict.items() if s and t in x_pos]
    if not items:
        return
    sig_xs = [x_pos[t] for t, _ in items]
    x_lo   = min([ref_x] + sig_xs)
    x_hi   = max([ref_x] + sig_xs)

    ax.plot([x_lo, x_hi], [rail_y, rail_y],
            color="black", lw=BRACKET_LW, clip_on=False, zorder=5)
    ax.plot([ref_x, ref_x], [rail_y, rail_y + tick_h],
            color="black", lw=BRACKET_LW, clip_on=False, zorder=5)

    for tx, st in items:
        t_x = x_pos[tx]
        ax.plot([t_x, t_x], [rail_y, rail_y + tick_h],
                color="black", lw=BRACKET_LW, clip_on=False, zorder=5)
        ax.text(t_x, rail_y - tick_h * 2.0, st,
                ha="center", va="top", fontsize=fontsize,
                fontweight="bold", clip_on=False, zorder=5)

# =============================================================================
# SECTION 3.2
# =============================================================================

# LMM forest plot with Option B brackets
def fig_3_2_1():
    lmm = q("01_primary_LMM", "lmm_log_results.csv")[
        ["treatment", "beta", "ci_lo", "ci_hi", "p"]].dropna(subset=["beta"])
    order = [t for t in TX5 if t in lmm["treatment"].values]
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.axvline(0, color="gray", ls="--", lw=1, zorder=1)
    rng = np.random.default_rng(42)
    ref_y = -0.5
    _filled(ax, 0, ref_y, "o", TX_COLORS["mScarlet"], False, s=120)
    y_pos, sig_dict = {}, {}
    for i, tx in enumerate(order):
        row  = lmm[lmm["treatment"] == tx].iloc[0]
        c, m = TX_COLORS[tx], TX_MARKERS[tx]
        rats = ANIM[ANIM["treatment"] == tx]["ratio"].values
        if len(rats):
            lb  = np.log(rats) - MS_LOGMEAN
            jit = rng.uniform(-0.15, 0.15, len(lb))
            ax.scatter(lb, i + jit, s=38, c=c, alpha=0.75,
                       edgecolors="white", linewidths=0.5, marker=m, zorder=2)
        ax.errorbar(row["beta"], i,
                    xerr=[[row["beta"] - row["ci_lo"]], [row["ci_hi"] - row["beta"]]],
                    fmt="none", color="black",
                    capsize=4, capthick=1.5, linewidth=1.5, zorder=3)
        _filled(ax, row["beta"], i, m, c, row["p"] < 0.05)
        y_pos[tx] = i
        if row["p"] < 0.05:
            sig_dict[tx] = stars(row["p"])
    rail_x = max(lmm[lmm["treatment"].isin(order)]["ci_hi"].max(), 0.05) + 0.06
    bracket_B_h(ax, y_pos, sig_dict, ref_y=ref_y, rail_x=rail_x, step=0.10)
    ax.set_yticks(list(range(len(order))))
    ax.set_yticklabels([""] * len(order))
    ax.set_ylim(ref_y - 0.3, len(order) - 0.4)
    ax.invert_yaxis()
    ax.set_xlabel("")
    ax.set_title("")
    plt.tight_layout()
    return fig

print("\nGenerating fig_3_2_1 ...")
fig_3_2_1().savefig(f"{OUTPUT_DIR}/fig_3_2_1.png")
plt.close("all")


# Animal-level ratios with Option B brackets
def fig_3_2_2():
    ms_rats    = ANIM[ANIM["treatment"] == "mScarlet"]["ratio"].values
    ms_logmean = np.log(ms_rats).mean()
    fig, ax    = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.axvline(1, color="gray", ls="--", lw=1, zorder=1)
    rng = np.random.default_rng(42)
    ref_y = -0.5
    y_pos, sig_dict = {}, {}
    for i, tx in enumerate(TX_ORDER):
        vals = ANIM[ANIM["treatment"] == tx]["ratio"].values
        c, m = TX_COLORS[tx], TX_MARKERS[tx]
        if tx == "mScarlet":
            _filled(ax, np.mean(vals), i, "o", c, False, s=160)
            continue
        log_vals = np.log(vals)
        log_ms   = np.log(ms_rats)
        diff     = log_vals.mean() - log_ms.mean()
        se_d     = np.sqrt(log_vals.var(ddof=1) / len(log_vals) +
                           log_ms.var(ddof=1) / len(log_ms))
        df_w     = len(log_vals) + len(log_ms) - 2
        t_stat   = diff / se_d if se_d > 0 else 0
        p_val    = 2 * sp_stats.t.sf(abs(t_stat), df_w)
        est      = np.exp(ms_logmean + diff)
        ci_lo    = np.exp(ms_logmean + diff - sp_stats.t.ppf(0.975, df_w) * se_d)
        ci_hi    = np.exp(ms_logmean + diff + sp_stats.t.ppf(0.975, df_w) * se_d)
        ax.errorbar(est, i, xerr=[[est - ci_lo], [ci_hi - est]],
                    fmt="none", color="black",
                    capsize=5, capthick=2, linewidth=2, zorder=2)
        _filled(ax, est, i, m, c, p_val < 0.05, s=150)
        y_pos[tx] = i
        if p_val < 0.05:
            sig_dict[tx] = stars(p_val)
        jit = rng.uniform(-0.15, 0.15, len(vals))
        ax.scatter(vals, i + jit, s=60, c=c, alpha=0.9,
                   marker=m, edgecolors="black", linewidths=0.5, zorder=4)
    rail_x = ANIM["ratio"].max() + 0.08
    bracket_B_h(ax, y_pos, sig_dict, ref_y=ref_y, rail_x=rail_x, step=0.05)
    ax.set_yticks(range(len(TX_ORDER)))
    ax.set_yticklabels([""] * len(TX_ORDER))
    ax.set_ylim(-0.6, len(TX_ORDER) - 0.4)
    ax.invert_yaxis()
    ax.set_xlabel("")
    ax.set_title("")
    plt.tight_layout()
    return fig

print("Generating fig_3_2_2 ...")
fig_3_2_2().savefig(f"{OUTPUT_DIR}/fig_3_2_2.png")
plt.close("all")


# Layer-stratified (4 panels, Option B brackets)
def fig_3_2_3():
    ld = q("03_atlas_stratified", "layer_contrasts.csv")[
        ["treatment", "layer_model", "beta", "ci_lo", "ci_hi", "p", "p_adj"]
    ].dropna(subset=["treatment", "beta"])
    layers = ["layer 1", "layer 2/3", "layer 5", "layer 6"]
    fig, axes = plt.subplots(1, 4, figsize=(12, FIG_H), sharey=True,
                              gridspec_kw={"wspace": 0.06})
    for ax, layer in zip(axes, layers):
        ldf   = ld[ld["layer_model"] == layer]
        ax.axvline(0, color="gray", ls="--", lw=1, zorder=1)
        ref_y = -0.5
        _filled(ax, 0, ref_y, "o", TX_COLORS["mScarlet"], False, s=80)
        y_pos, sig_dict = {}, {}
        for i, tx in enumerate(TX5):
            row = ldf[ldf["treatment"] == tx]
            if len(row) == 0:
                continue
            row    = row.iloc[0]
            c, m   = TX_COLORS[tx], TX_MARKERS[tx]
            is_sig = pd.notna(row["p_adj"]) and row["p_adj"] < 0.05
            ax.errorbar(row["beta"], i,
                        xerr=[[row["beta"] - row["ci_lo"]], [row["ci_hi"] - row["beta"]]],
                        fmt="none", color="black",
                        capsize=4, capthick=1.5, linewidth=1.5, zorder=2)
            _filled(ax, row["beta"], i, m, c, is_sig, s=100)
            y_pos[tx] = i
            if is_sig:
                sig_dict[tx] = stars(row["p_adj"])
        bracket_B_h(ax, y_pos, sig_dict, ref_y=ref_y,
                    rail_x=0.85, step=0.50, tick_w=0.05, fontsize=11)
        ax.set_yticks(range(len(TX5)))
        ax.set_yticklabels([""] * len(TX5))
        ax.set_xlim(-2.5, 2.0)
        ax.set_ylim(ref_y - 0.3, len(TX5) - 0.4)
        ax.invert_yaxis()
        ax.set_title("")
    plt.tight_layout()
    return fig

print("Generating fig_3_2_3 ...")
fig_3_2_3().savefig(f"{OUTPUT_DIR}/fig_3_2_3.png")
plt.close("all")

# =============================================================================
# SECTION 3.3
# =============================================================================

# Distance-decay curves (line plot)
def fig_3_3_1():
    bd = q("15_spatial_gradient", "binned_gradient.csv")[
        ["treatment", "dist_bin_mid", "frac_enwrapped", "n_cells"]].dropna()
    fig, ax = plt.subplots(figsize=(FIG_W, 4.8))
    for tx in TX_ORDER:
        c, m = TX_COLORS[tx], TX_MARKERS[tx]
        grp  = bd[(bd["treatment"] == tx) &
                  (bd["dist_bin_mid"] <= 0.8)].sort_values("dist_bin_mid")
        if len(grp) < 4:
            continue
        spl  = UnivariateSpline(grp["dist_bin_mid"], grp["frac_enwrapped"],
                                s=0.05, k=3)
        xmin, xmax = grp["dist_bin_mid"].min(), grp["dist_bin_mid"].max()
        xs = np.linspace(xmin, xmax, 200)
        ax.plot(xs, spl(xs), color=c, lw=1.8,
                ls="--" if tx == "mScarlet" else "-", zorder=3)
        xm = np.linspace(xmin, xmax, 6)
        ax.scatter(xm, spl(xm), s=52, marker=m, facecolor=c,
                   edgecolor="black", linewidth=0.8, zorder=4)
    ax.set_ylim(0, 0.55)
    ax.set_xlim(0, 0.82)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("")
    plt.tight_layout()
    return fig

print("Generating fig_3_3_1 ...")
fig_3_3_1().savefig(f"{OUTPUT_DIR}/fig_3_3_1.png")
plt.close("all")


# Off-target specificity (Bayesian)
def fig_3_3_s1():
    spec = q("13_bayesian_lmm", "m6_specificity.csv")[
        ["treatment", "post_mean_visual", "post_mean_offtarget",
         "p_visual_gt_offtarget"]].dropna(subset=["post_mean_visual"])
    order = list(reversed(TX5))

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.axvline(0, color="gray", ls="--", lw=1, zorder=1)
    OFF = 0.18

    y_pos, sig_dict = {}, {}
    for i, tx in enumerate(order):
        row = spec[spec["treatment"] == tx]
        if len(row) == 0:
            continue
        row    = row.iloc[0]
        c, m   = TX_COLORS[tx], TX_MARKERS[tx]
        is_sig = row["p_visual_gt_offtarget"] >= 0.95

        _filled(ax, row["post_mean_visual"],   i + OFF, m, c, is_sig, s=120)
        ax.scatter(row["post_mean_offtarget"], i - OFF,
                   s=120, c="white", marker=m,
                   edgecolors="#AAAAAA", linewidths=1.4, zorder=4)
        ax.plot([row["post_mean_offtarget"], row["post_mean_visual"]],
                [i - OFF, i + OFF],
                color=c, lw=1.0, alpha=0.5, zorder=2)

        y_pos[tx] = i + OFF
        if is_sig:
            extreme_p = 1 - row["p_visual_gt_offtarget"]
            sig_dict[tx] = stars(extreme_p) if extreme_p < 0.05 else "†"

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([""] * len(order))
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.set_xlabel("")
    ax.set_title("")
    plt.tight_layout()
    return fig

print("Generating fig_3_3_s1 ...")
fig_3_3_s1().savefig(f"{OUTPUT_DIR}/fig_3_3_s1.png")
plt.close("all")

# =============================================================================
# SECTION 3.4
# =============================================================================

# Factorial interaction (model prediction)
def fig_3_4_1():
    m8 = q("23_bayesian_additive", "m8_factorial_result.csv")[
        ["parameter", "post_mean"]].dropna(subset=["post_mean"])

    def _g(p):
        return m8.loc[m8["parameter"] == p, "post_mean"].values[0]

    intercept = _g("ADAMTS4 (intercept)")
    a15       = _g("ADAMTS15_present")
    c6        = _g("C6ST1_present")
    ix        = _g("Interaction (A15 x C6ST1)")
    p_a4      = intercept
    p_c6      = intercept + c6
    p_a15     = intercept + a15
    p_combo   = intercept + a15 + c6 + ix
    p_add     = intercept + a15 + c6

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.axhline(0, color="#CCCCCC", lw=0.8, ls="--", zorder=0)
    ax.plot([0, 1], [p_a4, p_c6],
            color=TX_COLORS["ADAMTS4"], lw=2.5,
            marker=TX_MARKERS["ADAMTS4"], markersize=9,
            markeredgecolor="black", markeredgewidth=1.0, zorder=3)
    ax.plot([0, 1], [p_a15, p_combo],
            color=TX_COLORS["C6ST1_ADAMTS15"], lw=2.5,
            marker=TX_MARKERS["C6ST1_ADAMTS15"], markersize=9,
            markeredgecolor="black", markeredgewidth=1.0, zorder=3)
    ax.plot(1, p_add, marker="x", markersize=11, mew=2.5, color="black", zorder=5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["", ""])
    ax.set_xlim(-0.25, 1.55)
    ax.set_ylim(-0.88, 0.08)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("")
    plt.tight_layout()
    return fig

print("Generating fig_3_4_1 ...")
fig_3_4_1().savefig(f"{OUTPUT_DIR}/fig_3_4_1.png")
plt.close("all")

# =============================================================================
# SECTION 3.5
# =============================================================================

# Zone profile ADAMTS-4 vs ADAMTS-4-MD
def fig_3_5_1():
    zl = q("00b_enwrapment_by_zone", "zone_lmm_results.csv")[
        ["treatment", "zone", "beta", "ci_lo", "ci_hi", "p"]
    ].dropna(subset=["beta"])
    zl = zl[zl["treatment"].isin(["ADAMTS4", "ADAMTS4_MD"]) &
            zl["zone"].isin(["Core", "Penumbra", "Outside"])]

    zones = ["Core", "Penumbra", "Outside"]
    zx    = {z: i for i, z in enumerate(zones)}
    OFF   = 0.10

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.axhline(0, color="#AAAAAA", ls="--", lw=1, zorder=1)
    ax.axhspan(-0.08, 0.08, color="#F0F0F0", zorder=0)

    for z in zones:
        _filled(ax, zx[z], 0, "o", TX_COLORS["mScarlet"], False, s=80)

    for tx, sign in [("ADAMTS4", -1), ("ADAMTS4_MD", +1)]:
        c, m = TX_COLORS[tx], TX_MARKERS[tx]
        grp  = zl[zl["treatment"] == tx]
        xs, ys = [], []
        for z in zones:
            row = grp[grp["zone"] == z]
            if len(row) == 0:
                continue
            r = row.iloc[0]
            x = zx[z] + sign * OFF
            xs.append(x)
            ys.append(r["beta"])
        ax.plot(xs, ys, color=c, lw=1.6, alpha=0.7, zorder=2)
        for x, y, z in zip(xs, ys, zones):
            r = grp[grp["zone"] == z].iloc[0]
            _filled(ax, x, y, m, c, r["p"] < 0.05, s=160)

    ax.set_xticks(list(zx.values()))
    ax.set_xticklabels([""] * 3)
    ax.set_xlim(-0.55, 2.70)
    ax.set_ylim(-0.45, 0.30)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("")
    plt.tight_layout()
    return fig

print("Generating fig_3_5_1 ...")
fig_3_5_1().savefig(f"{OUTPUT_DIR}/fig_3_5_1.png")
plt.close("all")


# =============================================================================
# SECTION 3.7
# =============================================================================

# Fig Virus-positive vs virus-negative enwrapment (Bayesian)
def fig_3_7_1():
    bv = q("24_virus_cell_enwrapment", "partD_bayesian_results.csv")[
        ["treatment", "post_mean_diff", "ci_89_lo", "ci_89_hi",
         "p_viruspos_lower"]].dropna(subset=["post_mean_diff"])

    order = list(reversed(TX5))

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.axvline(0, color="#AAAAAA", ls="--", lw=1, zorder=1)
    ax.axvspan(-0.02, 0.02, color="#F5F5F5", zorder=0)

    ref_y = len(order)
    ax.scatter(0, ref_y, s=120, marker="o",
               c=TX_COLORS["mScarlet"], edgecolors="black",
               linewidths=1.1, zorder=4)

    for i, tx in enumerate(order):
        row = bv[bv["treatment"] == tx]
        if len(row) == 0:
            continue
        row  = row.iloc[0]
        c, m = TX_COLORS[tx], TX_MARKERS[tx]
        d    = row["post_mean_diff"]
        p    = row["p_viruspos_lower"]
        is_cred = (p > 0.95) or (p < 0.05)
        ax.errorbar(d, i,
                    xerr=[[d - row["ci_89_lo"]], [row["ci_89_hi"] - d]],
                    fmt="none", color=c,
                    capsize=4, capthick=1.6, linewidth=1.6, zorder=2)
        _filled(ax, d, i, m, c, is_cred, s=160)

    ax.set_yticks(list(range(len(order))) + [ref_y])
    ax.set_yticklabels([""] * (len(order) + 1))
    ax.set_ylim(-0.7, ref_y + 0.5)
    ax.set_xlim(-0.20, 0.16)
    ax.set_xlabel("")
    ax.set_title("")
    plt.tight_layout()
    return fig

print("Generating fig_3_7_1 ...")
fig_3_7_1().savefig(f"{OUTPUT_DIR}/fig_3_7_1.png")
plt.close("all")


# Fig Distance-bin enwrapment (line plot)
def fig_3_7_2():
    bd = q("25_cell_level_deeper", "partB_binned_enwrapment.csv")[
        ["treatment", "dist_bin", "frac_enwrapped", "n_cells"]].dropna()
    bin_order = ["0\u201330px", "30\u2013100px", "100\u2013300px",
                 "300\u2013600px", "600px+"]
    bx  = {b: i for i, b in enumerate(bin_order)}
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ms = bd[bd["treatment"] == "mScarlet"]
    xs = [bx[b] for b in bin_order if len(ms[ms["dist_bin"] == b]) > 0]
    ys = [ms[ms["dist_bin"] == b].iloc[0]["frac_enwrapped"]
          for b in bin_order if len(ms[ms["dist_bin"] == b]) > 0]
    ax.plot(xs, ys, color=TX_COLORS["mScarlet"], lw=1.8, ls="--", zorder=2)
    ax.scatter(xs, ys, s=52, marker=TX_MARKERS["mScarlet"],
               facecolor=TX_COLORS["mScarlet"], edgecolor="black",
               linewidth=0.8, zorder=3)
    for tx in TX5:
        c, m = TX_COLORS[tx], TX_MARKERS[tx]
        grp  = bd[bd["treatment"] == tx]
        xs, ys = [], []
        for b in bin_order:
            row = grp[grp["dist_bin"] == b]
            if len(row) > 0:
                xs.append(bx[b])
                ys.append(row.iloc[0]["frac_enwrapped"])
        if len(xs) >= 2:
            ax.plot(xs, ys, color=c, lw=1.8, zorder=3)
            ax.scatter(xs, ys, s=52, marker=m, facecolor=c,
                       edgecolor="black", linewidth=0.8, zorder=4)
    ax.set_xticks(range(len(bin_order)))
    ax.set_xticklabels([""] * len(bin_order))
    ax.set_ylim(0, 0.55)
    ax.set_xlim(-0.3, len(bin_order) - 0.7)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("")
    plt.tight_layout()
    return fig

print("Generating fig_3_7_2 ...")
fig_3_7_2().savefig(f"{OUTPUT_DIR}/fig_3_7_2.png")
plt.close("all")

# =============================================================================
# SECTION 3.8
# =============================================================================

# Virus-status enwrapment comparison (fill = virus status)
def fig_3_8():
    FOCUS    = ["mScarlet", "C6ST1", "ADAMTS15"]
    animal_p = ANIM[ANIM["treatment"].isin(FOCUS)]
    ms_mean  = ANIM[ANIM["treatment"] == "mScarlet"]["ratio"].mean()
    tc = q("24_virus_cell_enwrapment", "partA_triple_colocalization.csv")
    tc = tc[tc["treatment"].isin(["C6ST1", "ADAMTS15"])]
    tc["frac_vp"] = pd.to_numeric(tc["pct_viruspos_enwrap"], errors="coerce") / 100
    tc["frac_vn"] = pd.to_numeric(tc["pct_virusneg_enwrap"], errors="coerce") / 100
    paired  = tc.groupby(["treatment", "mouse_id"])[["frac_vp", "frac_vn"]].mean().reset_index()
    grp_x   = {"C6ST1": {"vn": 0.0, "vp": 0.6}, "ADAMTS15": {"vn": 1.5, "vp": 2.1}}
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(FIG_W, FIG_H),
                                    gridspec_kw={"width_ratios": [1.0, 1.3], "wspace": 0.38})
    rng = np.random.default_rng(42)
    axA.axhline(1.0, color="#DDDDDD", ls=":", lw=0.8, zorder=0)
    axA.axhline(ms_mean, color="#AAAAAA", ls="--", lw=1.0, zorder=1)
    for i, tx in enumerate(FOCUS):
        vals = animal_p[animal_p["treatment"] == tx]["ratio"].values
        c, m = TX_COLORS[tx], TX_MARKERS[tx]
        axA.scatter(i, np.mean(vals), s=200, marker=m, c=c,
                    edgecolors="black", linewidths=1.5, zorder=3)
        jit = rng.uniform(-0.12, 0.12, len(vals))
        axA.scatter(i + jit, vals, s=50, c=c, alpha=0.85,
                    marker=m, edgecolors="black", linewidths=0.5, zorder=4)
    axA.set_xticks(range(len(FOCUS)))
    axA.set_xticklabels([""] * len(FOCUS))
    axA.set_xlim(-0.55, len(FOCUS) - 0.45)
    axA.set_ylim(0.45, 1.22)
    axA.set_ylabel("")
    axA.set_xlabel("")
    axA.set_title("")
    axB.axhline(0, color="#DDDDDD", ls=":", lw=0.8, zorder=0)
    for tx in ["C6ST1", "ADAMTS15"]:
        c, m = TX_COLORS[tx], TX_MARKERS[tx]
        xvn, xvp = grp_x[tx]["vn"], grp_x[tx]["vp"]
        for _, row in paired[paired["treatment"] == tx].iterrows():
            axB.plot([xvn, xvp], [row["frac_vn"], row["frac_vp"]],
                     color=c, lw=1.0, alpha=0.55, zorder=2)
            axB.scatter(xvn, row["frac_vn"], s=60, c="white", marker=m,
                        edgecolors=c, linewidths=1.5, zorder=3)
            axB.scatter(xvp, row["frac_vp"], s=60, c=c, marker=m,
                        edgecolors="black", linewidths=0.7, zorder=3)
    axB.set_xticks([])
    axB.set_xlim(-0.15, 2.75)
    axB.set_ylim(-0.06, 0.72)
    axB.set_ylabel("")
    axB.set_xlabel("")
    axB.set_title("")
    plt.tight_layout()
    return fig

print("Generating fig_3_8 ...")
fig_3_8().savefig(f"{OUTPUT_DIR}/fig_3_8.png")
plt.close("all")

# =============================================================================
# SECTION 3.11
# =============================================================================

# Fig Viral load dose-response curves by treatment
def fig_viral_load_lines():
    sl = q("28_viral_load_heterogeneity", "analysis1_per_tx_slopes.csv")[
        ["treatment", "slope"]].dropna(subset=["slope"])

    mean_enw = (ANIM.groupby("treatment")["ratio"].mean() * 0.43).to_dict()

    def _logit(p):
        return np.log(np.clip(p, 1e-6, 1 - 1e-6) / (1 - np.clip(p, 1e-6, 1 - 1e-6)))

    def _expit(x):
        return 1 / (1 + np.exp(-x))

    x         = np.linspace(-3, 3, 300)
    x_ticks   = [-2, -1, 0, 1, 2]
    x_markers = np.array(x_ticks)

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.axvline(x=0, color="#000000", ls=":", lw=0.8, zorder=1)

    for tx in TX_ORDER:
        row = sl[sl["treatment"] == tx]
        if len(row) == 0:
            continue
        slope     = row.iloc[0]["slope"]
        color     = TX_COLORS[tx]
        mkr       = TX_MARKERS[tx]
        anchor    = mean_enw.get(tx, 0.35)
        intercept = _logit(anchor)

        x_plot = np.linspace(-3, 3, 300) if tx == "mScarlet" else x
        x_m    = np.array([-2, -1, 0, 1, 2]) if tx == "mScarlet" else x_markers

        y         = _expit(intercept + slope * x_plot)
        y_markers = _expit(intercept + slope * x_m)

        ax.plot(x_plot, y, color=color, linewidth=1.8,
                ls="--" if tx == "mScarlet" else "-", zorder=3)
        ax.scatter(x_m, y_markers, s=48, marker=mkr,
                   facecolor=color, edgecolor="black",
                   linewidth=0.8, zorder=4)

    ax.set_xlim(-3.2, 3.2)
    ax.set_ylim(0.08, 0.62)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([str(v) for v in x_ticks], fontsize=10)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("")
    plt.tight_layout()
    return fig

print("Generating fig_viral_load_lines ...")
fig_viral_load_lines().savefig(f"{OUTPUT_DIR}/fig_viral_load_lines.png")
plt.close("all")

# =============================================================================
# SUPPLEMENTARY — SECTION 3.2
# =============================================================================

# WFA cell density (hemisphere-level LMM)
def fig_3_2_s1():
    dens  = q("11_density_lmm", "primary_lmm_results.csv")
    wfa_d = dens[dens["staining"] == "WFA"][
        ["term", "beta", "p", "p_adj"]].dropna(subset=["beta"])
    wfa_d["treatment"] = wfa_d["term"].str.replace("treatment", "", regex=False)
    order = [t for t in TX5 if t in wfa_d["treatment"].values]

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.axvline(0, color="gray", ls="--", lw=1, zorder=1)

    ref_y = -0.5
    _filled(ax, 0, ref_y, "o", TX_COLORS["mScarlet"], False, s=120)

    y_pos, sig_dict = {}, {}
    for i, tx in enumerate(order):
        row = wfa_d[wfa_d["treatment"] == tx]
        if len(row) == 0:
            continue
        row    = row.iloc[0]
        c, m   = TX_COLORS[tx], TX_MARKERS[tx]
        is_sig = pd.notna(row["p_adj"]) and row["p_adj"] < 0.05
        _filled(ax, row["beta"], i, m, c, is_sig)
        y_pos[tx] = i
        if is_sig:
            sig_dict[tx] = stars(row["p_adj"])

    bracket_B_h(ax, y_pos, sig_dict, ref_y=ref_y, rail_x=0.22)

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([""] * len(order))
    ax.set_xlim(-0.95, 0.40)
    ax.set_ylim(ref_y - 0.3, len(order) - 0.4)
    ax.invert_yaxis()
    ax.set_xlabel("")
    ax.set_title("")
    plt.tight_layout()
    return fig

print("Generating fig_3_2_s1 ...")
fig_3_2_s1().savefig(f"{OUTPUT_DIR}/fig_3_2_s1.png")
plt.close("all")

# =============================================================================
# SUPPLEMENTARY — SECTION 3.7
# =============================================================================

# Viral load dose-response slopes
def fig_3_7_s1():
    sl = q("28_viral_load_heterogeneity", "analysis1_per_tx_slopes.csv")[
        ["treatment", "slope", "SE", "p", "OR_per_unit"]
    ].dropna(subset=["slope"])
    order = [t for t in TX_ORDER if t in sl["treatment"].values]

    fig, ax = plt.subplots(figsize=(FIG_W, 4.2))
    ax.axvline(0, color="#AAAAAA", ls="--", lw=1, zorder=1)

    y_pos, sig_dict = {}, {}
    for i, tx in enumerate(order):
        row    = sl[sl["treatment"] == tx].iloc[0]
        c, m   = TX_COLORS[tx], TX_MARKERS[tx]
        is_sig = (row["p"] < 0.05) and (tx != "mScarlet")
        ax.errorbar(row["slope"], i, xerr=1.96 * row["SE"],
                    fmt="none", color=c,
                    capsize=4, capthick=1.5, linewidth=1.5, zorder=2)
        _filled(ax, row["slope"], i, m, c,
                is_sig if tx != "mScarlet" else False, s=140)
        y_pos[tx] = i
        if is_sig:
            sig_dict[tx] = stars(row["p"])

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([""] * len(order))
    ax.set_xlim(-0.85, 0.22)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.invert_yaxis()
    ax.set_xlabel("")
    ax.set_title("")
    plt.tight_layout()
    return fig

print("Generating fig_3_7_s1 ...")
fig_3_7_s1().savefig(f"{OUTPUT_DIR}/fig_3_7_s1.png")
plt.close("all")

# =============================================================================
# APPENDIX D — PRELIMINARY ANALYSES
# =============================================================================

# Hemisphere WFA and PV fluorescence (preliminary)
s3 = pd.read_csv(os.path.join(PRE_LMM_DIR, "S3_hemisphere.csv"))
s3["p_adj"]     = pd.to_numeric(s3["p_adj"],     errors="coerce")
s3["mean_diff"] = pd.to_numeric(s3["mean_diff"], errors="coerce")

metric_offsets = {"diffFluo": 0.15, "avgPxIntensity": -0.15}
plot_order_D1  = [t for t in treatment_order if t != "mScarlet"]


def fig_D1():
    fig, axes = plt.subplots(1, 2, figsize=(FIG_WIDTH, FIG_HEIGHT),
                              sharey=True, gridspec_kw={"wspace": 0.08})
    for ax, staining in zip(axes, ["WFA", "PV"]):
        ax.axvline(x=0, color="gray", linestyle="--", linewidth=1, zorder=1)
        sub = s3[s3["staining"] == staining]
        for i, treatment in enumerate(reversed(plot_order_D1)):
            color = treatment_colors[treatment]
            mkr   = treatment_markers[treatment]
            for metric, y_off in metric_offsets.items():
                row = sub[
                    (sub["treatment"] == treatment) &
                    (sub["metric"]    == metric)
                ]
                if len(row) == 0:
                    continue
                row    = row.iloc[0]
                is_sig = row["p_adj"] < 0.05
                ax.scatter(row["mean_diff"], i + y_off,
                           s=90, marker=mkr,
                           c=color if is_sig else "white",
                           edgecolors=color, linewidths=1.5, zorder=3)
        ax.set_yticks(range(len(plot_order_D1)))
        ax.set_yticklabels([""] * len(plot_order_D1))
        ax.set_ylim(-0.6, len(plot_order_D1) - 0.4)
        ax.set_xlabel("")
        ax.set_title("")
    plt.tight_layout()
    return fig

print("Generating fig_D1 ...")
fig_D1().savefig(f"{OUTPUT_DIR}/D1_hemisphere_fluorescence.png",
                 dpi=300, bbox_inches="tight")
plt.close("all")


# Zone-stratified virus+/virus- fluorescence (preliminary)
s1 = pd.read_csv(os.path.join(PRE_LMM_DIR, "S1_virus_intensity.csv"))
s1["p_raw"]     = pd.to_numeric(s1["p_raw"],     errors="coerce")
s1["mean_diff"] = pd.to_numeric(s1["mean_diff"], errors="coerce")

zone_offsets  = {"Core": 0.15, "Penumbra": -0.15}
plot_order_D2 = [t for t in treatment_order if t != "mScarlet"]


def fig_D2():
    fig, axes = plt.subplots(1, 2, figsize=(FIG_WIDTH, FIG_HEIGHT),
                              sharey=True, gridspec_kw={"wspace": 0.08})
    for ax, stain in zip(axes, ["PNN", "PV"]):
        ax.axvline(x=0, color="gray", linestyle="--", linewidth=1, zorder=1)
        sub = s1[s1["stain"] == stain]
        for i, treatment in enumerate(reversed(plot_order_D2)):
            color = treatment_colors[treatment]
            mkr   = treatment_markers[treatment]
            for zone, y_off in zone_offsets.items():
                row = sub[
                    (sub["treatment"] == treatment) &
                    (sub["zone"]      == zone)
                ]
                if len(row) == 0:
                    continue
                row    = row.iloc[0]
                is_sig = row["p_raw"] < 0.05
                ax.scatter(row["mean_diff"], i + y_off,
                           s=90, marker=mkr,
                           c=color if is_sig else "white",
                           edgecolors=color, linewidths=1.5, zorder=3)
        ax.set_yticks(range(len(plot_order_D2)))
        ax.set_yticklabels([""] * len(plot_order_D2))
        ax.set_ylim(-0.6, len(plot_order_D2) - 0.4)
        ax.set_xlabel("")
        ax.set_title("")
    plt.tight_layout()
    return fig

print("Generating fig_D2 ...")
fig_D2().savefig(f"{OUTPUT_DIR}/D2_zone_fluorescence.png",
                 dpi=300, bbox_inches="tight")
plt.close("all")


# Power analysis: Cohen's d per comparison
s5 = pd.read_csv(os.path.join(PRE_LMM_DIR, "S5_power_analysis.csv"))
s5["cohen_d"] = pd.to_numeric(s5["cohen_d"], errors="coerce")

s5_plot = s5[s5["section"] == "S2: PV\u2013PNN enwrapment"].copy().reset_index(drop=True)
s5_plot["treatment"] = s5_plot["label"].apply(lambda x: x.split()[-1])
s5_plot["zone"]      = s5_plot["label"].apply(lambda x: x.split()[0])

zone_stem_colors = {
    "Core":     "#2171B5",
    "Penumbra": "#FDAE6B",
    "Outside":  "#74C476",
}


def _color(t):
    for key in ["C6ST1_ADAMTS15", "ADAMTS4_MD", "ADAMTS15", "C6ST1", "ADAMTS4", "mScarlet"]:
        if key.lower() == t.strip().lower():
            return treatment_colors[key]
    return "#888888"


def _marker(t):
    for key in ["C6ST1_ADAMTS15", "ADAMTS4_MD", "ADAMTS15", "C6ST1", "ADAMTS4", "mScarlet"]:
        if key.lower() == t.strip().lower():
            return treatment_markers[key]
    return "o"


def fig_D3():
    fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))
    ax.axvline(x=1.9, color="#AAAAAA", linestyle="--", linewidth=1.2, zorder=0)
    ax.axvspan(0, 1.9, color="#F8F8F8", zorder=0)

    for i, (_, row) in enumerate(s5_plot.iterrows()):
        color = _color(row["treatment"])
        mkr   = _marker(row["treatment"])
        d     = row["cohen_d"]
        zcol  = zone_stem_colors.get(row["zone"], "#888888")
        ax.plot([0, d], [i, i], color=zcol, linewidth=1.2, alpha=0.6, zorder=1)
        ax.scatter(d, i, s=80, marker=mkr,
                   c=color, edgecolors="black",
                   linewidths=0.6, zorder=3)

    ax.set_yticks(range(len(s5_plot)))
    ax.set_yticklabels([""] * len(s5_plot))
    ax.set_ylim(-0.7, len(s5_plot) - 0.3)
    ax.set_xlim(0, 2.5)
    ax.set_xlabel("")
    ax.set_title("")
    plt.tight_layout()
    return fig

print("Generating fig_D3 ...")
fig_D3().savefig(f"{OUTPUT_DIR}/D3_power_analysis.png",
                 dpi=300, bbox_inches="tight")
plt.close("all")

# =============================================================================
# SECTION 3.1 — REPRESENTATIVE MICROSCOPY FIGURES
# =============================================================================
# These figures render multi-channel TIFF microscopy data. Configure the paths
# below. If any required input is missing, the entire microscopy section is
# skipped and the data-plot figures above still run.

# -- MICROSCOPY CONFIGURATION --
from pathlib import Path
import json as _json

TIFF_ROOT       = Path("/path/to/originals_renamed")
CELLS_CSV       = Path("/path/to/cells_with_zones.csv")
VISUALIGN_DIR   = Path("/path/to/Visualign_json")
ATLAS_DIR       = Path("/path/to/atlas_pngs")
RAINBOW_PATH    = Path("/path/to/Rainbow_2017.json")
MICROSCOPY_OUT  = Path(f"{OUTPUT_DIR}/microscopy")

# -- Channel indices and display colours --
CH_WFA, CH_PV, CH_MS = 0, 1, 2

# -- Reference slice for AP-matching (visual continuity with Figure 2.2) --
REF_ANIMAL = "mScarlet_3"
REF_SLICE  = "s051"

# -- Crop geometry --
CROP_HALF_10X    = 350        # 700x700 px ≈ 455 µm at 10x (0.65 µm/px)
MIN_PV_IN_WINDOW = 25
PX_PER_UM_10X    = 1 / 0.65

# -- Display normalisation (uniform across all six tiles) --
NORM_LO_WFA, NORM_HI_WFA, GAMMA_WFA = 1.0, 98.0, 0.85
NORM_LO_PV,  NORM_HI_PV,  GAMMA_PV  = 1.0, 98.5, 0.90
NORM_LO_MS,  NORM_HI_MS,  GAMMA_MS  = 1.0, 99.0, 0.80

# When True, WFA bounds are taken from the Control crop and applied as absolute
# values to all panels. This preserves real intensity differences between
# treatments instead of auto-stretching each tile to its own range.
USE_REFERENCE_NORM = True

# When True, blend mScarlet into the composite (cyan tint over transduced
# regions). When False, the composite is two-channel WFA + PV only.
SHOW_MSCARLET   = False
MS_GREEN_WEIGHT = 0.5
MS_BLUE_WEIGHT  = 0.5

# -- Anatomical filters --
TARGET_BRAIN_AREAS = [
    "Primary visual area, layer 1",
    "Primary visual area, layer 2/3",
    "Primary visual area, layer 4",
    "Primary visual area, layer 5",
    "Primary visual area, layer 6a",
    "Primary visual area, layer 6b",
]
PRIORITY_LAYERS = [
    "Primary visual area, layer 2/3",
    "Primary visual area, layer 4",
]
IPSI_HEMISPHERE = "left"

# -- Atlas overlay settings --
SHOW_ATLAS_OVERLAY   = True
SHOW_VISUAL_FILL     = True
ATLAS_EDGE_COLOR     = (1.0, 0.85, 0.20)
ATLAS_EDGE_ALPHA     = 0.65
ATLAS_EDGE_THICKNESS = 1
V1_FILL_COLOR        = (1.0, 0.85, 0.20)
V1_FILL_ALPHA        = 0.35
CROP_BOX_LW          = 1.5
CROP_BOX_COLOR       = "#00FFFF"
OVERVIEW_TARGET_W_PX = 1500
OVERVIEW_SCALEBAR_UM = 1000
OVERVIEW_ZOOM        = True
ZOOM_FRAC            = 0.30
ZOOM_SCALEBAR_UM     = 500

# -- Output and figure sizing --
MICROSCOPY_DPI    = 300
COMPOSITE_W_IN    = 12
SCALEBAR_UM       = 200
TILE_LABEL_FS     = 22
TILE_PCT_FS       = 20
COMP_LABEL_FS     = 14
COMP_PCT_FS       = 12

# -- Group selection --
# 'animal' = mouse_id, 'slice' = section ID (sNNN), 'ic_pct_label' = caption.
GROUPS = [
    dict(treatment="mScarlet",       label="Control",
         animal="mScarlet_3",  slice="s073", manual_crop_10x=None,
         ic_pct_label="ref."),
    dict(treatment="ADAMTS15",       label="ADAMTS-15",
         animal=None,          slice="s019", manual_crop_10x=None,
         ic_pct_label="−29.2%"),
    dict(treatment="C6ST1_ADAMTS15", label="C6ST-1 / ADAMTS-15",
         animal=None,          slice="s037", manual_crop_10x=None,
         ic_pct_label="−51.0%"),
    dict(treatment="ADAMTS4",        label="ADAMTS-4",
         animal=None,          slice="s035", manual_crop_10x=None,
         ic_pct_label="–19.1%"),
    dict(treatment="ADAMTS4_MD",     label="ADAMTS-4 + MD",
         animal=None,          slice="s045", manual_crop_10x=None,
         ic_pct_label="n.s."),
    dict(treatment="C6ST1",          label="C6ST-1",
         animal="C6ST1_2",     slice="s019", manual_crop_10x=None,
         ic_pct_label="n.s."),
]

# -- Skip the section entirely if microscopy inputs are absent --
_microscopy_inputs = [TIFF_ROOT, CELLS_CSV, VISUALIGN_DIR, ATLAS_DIR, RAINBOW_PATH]
_microscopy_available = all(p.exists() for p in _microscopy_inputs)

if not _microscopy_available:
    print("\n" + "=" * 60)
    print("Microscopy section skipped — required inputs not found:")
    for p in _microscopy_inputs:
        marker = "OK" if p.exists() else "MISSING"
        print(f"  [{marker}] {p}")
    print("=" * 60)
else:
    try:
        import tifffile
        from scipy.ndimage import laplace
        import matplotlib.patheffects as pe
        from PIL import Image
    except ImportError as _exc:
        print(f"\nMicroscopy section skipped — missing package: {_exc}")
        _microscopy_available = False

if _microscopy_available:
    MICROSCOPY_OUT.mkdir(parents=True, exist_ok=True)
    print("\n" + "=" * 60)
    print("Generating microscopy figures (§3.1)")
    print("=" * 60)

    # -- Helpers --
    def norm_gamma(img, lo_pct, hi_pct, gamma, *, abs_bounds=None):
        """Percentile-stretch to [0,1] then apply gamma.
        If abs_bounds=(lo_val, hi_val) is given, those raw-intensity values
        are used instead of computing percentiles from img."""
        if abs_bounds is not None:
            lo, hi = abs_bounds
        else:
            lo = np.percentile(img, lo_pct)
            hi = np.percentile(img, hi_pct)
        if hi == lo:
            return np.zeros_like(img, dtype=float)
        return np.clip((img.astype(float) - lo) / (hi - lo), 0, 1) ** gamma

    def compute_channel_bounds(img, lo_pct, hi_pct):
        return float(np.percentile(img, lo_pct)), float(np.percentile(img, hi_pct))

    def outline(color="black", lw=4):
        return [pe.withStroke(linewidth=lw, foreground=color)]

    def load_tiff_three_channel(path):
        raw = tifffile.imread(str(path))
        if raw.ndim == 3 and raw.shape[0] <= 5:
            return raw[CH_WFA], raw[CH_PV], raw[CH_MS]
        if raw.ndim == 3 and raw.shape[2] <= 5:
            return raw[:, :, CH_WFA], raw[:, :, CH_PV], raw[:, :, CH_MS]
        raise ValueError(f"Unexpected TIFF shape: {raw.shape}")

    def normalise_three_channel(c1, c2, c3, *, ref_bounds=None):
        """Apply the uniform per-channel display settings.
        If ref_bounds is provided, only the WFA channel uses absolute bounds;
        PV and mScarlet remain per-crop normalised so cells stay visible."""
        bw = ref_bounds.get("wfa") if ref_bounds else None
        c1n = norm_gamma(c1, NORM_LO_WFA, NORM_HI_WFA, GAMMA_WFA, abs_bounds=bw)
        c2n = norm_gamma(c2, NORM_LO_PV,  NORM_HI_PV,  GAMMA_PV)
        c3n = norm_gamma(c3, NORM_LO_MS,  NORM_HI_MS,  GAMMA_MS)
        return c1n, c2n, c3n

    def make_rgb(c1n, c2n, c3n):
        """RGB composite: WFA → G; PV → R+B (magenta); mScarlet → G+B (cyan)."""
        H, W = c1n.shape
        rgb = np.zeros((H, W, 3), dtype=float)
        if SHOW_MSCARLET:
            rgb[:, :, 0] = np.clip(c2n,                              0, 1)
            rgb[:, :, 1] = np.clip(c1n + MS_GREEN_WEIGHT * c3n,      0, 1)
            rgb[:, :, 2] = np.clip(c2n + MS_BLUE_WEIGHT  * c3n,      0, 1)
        else:
            rgb[:, :, 0] = np.clip(c2n, 0, 1)
            rgb[:, :, 1] = np.clip(c1n, 0, 1)
            rgb[:, :, 2] = np.clip(c2n, 0, 1)
        return rgb

    def focus_score(channel, cx, cy, half, H, W):
        x0 = max(0, cx - half); x1 = min(W, cx + half)
        y0 = max(0, cy - half); y1 = min(H, cy + half)
        patch = channel[y0:y1, x0:x1].astype(float)
        return float(laplace(patch).var()) if patch.size > 0 else 0.0

    def add_scalebar(ax, w, h, scale_um, px_per_um, color="white",
                     margin_y=0.08, lw=5, fs=14):
        bar_px = px_per_um * scale_um
        x1 = w * 0.96; x0 = x1 - bar_px
        y  = h - h * margin_y
        ax.plot([x0, x1], [y, y], color=color, lw=lw,
                solid_capstyle="butt", transform=ax.transData)
        ax.text((x0 + x1) / 2, y - h * 0.025, f"{scale_um} µm",
                color=color, ha="center", va="bottom",
                fontsize=fs, fontweight="bold",
                transform=ax.transData, path_effects=outline(lw=3))

    def mouse_id_to_folder(mouse_id):
        parts = mouse_id.rsplit("_", 1)
        if len(parts) == 2 and parts[1].isdigit():
            return parts[0]
        return mouse_id

    def tiff_path_for(mouse_id, slice_id):
        folder = mouse_id_to_folder(mouse_id)
        return TIFF_ROOT / folder / mouse_id / f"{mouse_id}_{slice_id}.tiff"

    def atlas_path_for(mouse_id, slice_id):
        return ATLAS_DIR / f"{mouse_id}_{slice_id}_nl.png"

    def load_atlas_rgb(path):
        return np.asarray(Image.open(path).convert("RGB"))

    def resize_atlas_nearest(atlas_rgb, target_w, target_h):
        return np.asarray(
            Image.fromarray(atlas_rgb).resize((target_w, target_h), Image.NEAREST)
        )

    def detect_region_boundaries(atlas_rgb, thickness=1):
        labels = (atlas_rgb[..., 0].astype(np.int32) << 16) | \
                 (atlas_rgb[..., 1].astype(np.int32) << 8) | \
                  atlas_rgb[..., 2].astype(np.int32)
        edges = np.zeros_like(labels, dtype=bool)
        edges[:-1, :] |= labels[:-1, :] != labels[1:, :]
        edges[1:,  :] |= labels[1:,  :] != labels[:-1, :]
        edges[:, :-1] |= labels[:, :-1] != labels[:, 1:]
        edges[:, 1:]  |= labels[:, 1:]  != labels[:, :-1]
        if thickness > 1:
            from scipy.ndimage import binary_dilation
            edges = binary_dilation(edges, iterations=thickness - 1)
        return edges

    def composite_edges(rgb, edges, edge_color, alpha):
        out = rgb.copy()
        for c in range(3):
            out[..., c] = np.where(
                edges,
                alpha * edge_color[c] + (1 - alpha) * out[..., c],
                out[..., c],
            )
        return out

    def make_region_mask(atlas_rgb, packed_colors):
        r = atlas_rgb[..., 0].astype(np.uint32)
        g = atlas_rgb[..., 1].astype(np.uint32)
        b = atlas_rgb[..., 2].astype(np.uint32)
        return np.isin((r << 16) | (g << 8) | b, packed_colors)

    def block_downsample(arr, factor):
        if factor <= 1:
            return arr
        H, W = arr.shape
        H2 = (H // factor) * factor
        W2 = (W // factor) * factor
        return arr[:H2, :W2].reshape(H2 // factor, factor,
                                     W2 // factor, factor).mean(axis=(1, 3))

    # -- Load cells_with_zones.csv for crop centring --
    print("\nLoading cells_with_zones.csv ...")
    cells = pd.read_csv(
        CELLS_CSV,
        usecols=["mouse_id", "slice_id", "cell_type",
                 "x_hires", "y_hires", "zone", "brain_area", "hemisphere"],
        low_memory=False,
    )
    print(f"  {len(cells)} cells")

    # -- Auto-fill any animal=None by per-group median I/C ratio --
    print("\nResolving animals ...")
    _animal_data = ANIM.groupby(["treatment", "mouse_id"])["ratio"].mean().reset_index()
    for g in GROUPS:
        if g["animal"] is None:
            sub = _animal_data[_animal_data["treatment"] == g["treatment"]].sort_values("ratio")
            if len(sub):
                median_idx = len(sub) // 2 if len(sub) % 2 == 1 else (len(sub) - 1) // 2
                g["animal"] = sub.iloc[median_idx]["mouse_id"]
                print(f"  {g['label']:25s} → {g['animal']}  (median)")
            else:
                print(f"  {g['label']:25s} → NO ANIMALS in {g['treatment']}")
        else:
            print(f"  {g['label']:25s} → {g['animal']}  (set)")

    # -- Find crop centre helper (uses loaded `cells` global) --
    def find_crop_centre(c2_10x, mouse_id, slice_id, H, W):
        sl = cells[(cells["mouse_id"] == mouse_id) &
                   (cells["slice_id"] == slice_id)]
        pv = sl[(sl["cell_type"] == "PV") &
                (sl["hemisphere"] == IPSI_HEMISPHERE) &
                (sl["brain_area"].isin(TARGET_BRAIN_AREAS)) &
                (sl["zone"].isin(["Core", "Penumbra"]))].copy()

        if len(pv) < MIN_PV_IN_WINDOW:
            pv = sl[(sl["cell_type"] == "PV") &
                    (sl["hemisphere"] == IPSI_HEMISPHERE) &
                    (sl["brain_area"].isin(TARGET_BRAIN_AREAS))].copy()

        pv_priority = pv[pv["brain_area"].isin(PRIORITY_LAYERS)]
        if len(pv_priority) >= MIN_PV_IN_WINDOW:
            pv = pv_priority

        margin = CROP_HALF_10X + 20
        pv = pv[(pv["x_hires"] > margin) & (pv["x_hires"] < W - margin) &
                (pv["y_hires"] > margin) & (pv["y_hires"] < H - margin)]
        if len(pv) == 0:
            return W // 2, H // 2

        candidates = []
        for _, r in pv.iterrows():
            cx, cy = int(r["x_hires"]), int(r["y_hires"])
            in_win = ((pv["x_hires"] >= cx - CROP_HALF_10X) &
                      (pv["x_hires"] <  cx + CROP_HALF_10X) &
                      (pv["y_hires"] >= cy - CROP_HALF_10X) &
                      (pv["y_hires"] <  cy + CROP_HALF_10X))
            n_pv = int(in_win.sum())
            if n_pv < MIN_PV_IN_WINDOW:
                continue
            f = focus_score(c2_10x, cx, cy, CROP_HALF_10X, H, W)
            candidates.append((cx, cy, n_pv, f))

        if not candidates:
            return int(pv["x_hires"].median()), int(pv["y_hires"].median())

        cdf = pd.DataFrame(candidates, columns=["cx", "cy", "n_pv", "focus"])
        cdf["score"] = ((cdf["n_pv"]  - cdf["n_pv"].mean())  / (cdf["n_pv"].std()  + 1e-6)
                      + (cdf["focus"] - cdf["focus"].mean()) / (cdf["focus"].std() + 1e-6))
        best = cdf.sort_values("score", ascending=False).iloc[0]
        return int(best["cx"]), int(best["cy"])

    # -- Compute reference bounds from the Control crop --
    def compute_reference_bounds():
        g0 = GROUPS[0]
        c1r, c2r, c3r = load_tiff_three_channel(tiff_path_for(g0["animal"], g0["slice"]))
        Hr, Wr = c1r.shape
        cxr, cyr = (g0["manual_crop_10x"] if g0["manual_crop_10x"]
                    else find_crop_centre(c2r, g0["animal"], g0["slice"], Hr, Wr))
        cxr = max(CROP_HALF_10X, min(Wr - CROP_HALF_10X, cxr))
        cyr = max(CROP_HALF_10X, min(Hr - CROP_HALF_10X, cyr))
        x0r, x1r = cxr - CROP_HALF_10X, cxr + CROP_HALF_10X
        y0r, y1r = cyr - CROP_HALF_10X, cyr + CROP_HALF_10X
        return {
            "wfa": compute_channel_bounds(c1r[y0r:y1r, x0r:x1r], NORM_LO_WFA, NORM_HI_WFA),
            "pv":  compute_channel_bounds(c2r[y0r:y1r, x0r:x1r], NORM_LO_PV,  NORM_HI_PV),
            "ms":  compute_channel_bounds(c3r[y0r:y1r, x0r:x1r], NORM_LO_MS,  NORM_HI_MS),
        }

    REF_BOUNDS = compute_reference_bounds() if USE_REFERENCE_NORM else None
    if REF_BOUNDS:
        print(f"\nReference bounds (from {GROUPS[0]['label']}):")
        for ch, (lo, hi) in REF_BOUNDS.items():
            print(f"  {ch:>3s}: [{lo:.1f}, {hi:.1f}]")

    # -- Build V1 colour set from Rainbow_2017 (shared by overview figures) --
    with open(RAINBOW_PATH) as _f:
        _rb = _json.load(_f)
    V1_COLORS = np.array(
        [(e["red"] << 16) | (e["green"] << 8) | e["blue"]
         for e in _rb if "Primary visual area" in e["name"]],
        dtype=np.uint32,
    )
    print(f"\nV1 colour set: {len(V1_COLORS)} entries")

    # -- Standardised display size for whole-section overviews --
    _widths, _heights = [], []
    for _g in GROUPS:
        with tifffile.TiffFile(str(tiff_path_for(_g["animal"], _g["slice"]))) as _tif:
            _s = sorted(_tif.series[0].shape)
            _heights.append(_s[-2])
            _widths.append(_s[-1])
    OVERVIEW_FACTOR = max(1, max(_widths) // OVERVIEW_TARGET_W_PX)
    OVERVIEW_STD_W  = max(_widths)  // OVERVIEW_FACTOR
    OVERVIEW_STD_H  = max(_heights) // OVERVIEW_FACTOR
    print(f"Overview display: {OVERVIEW_STD_W} x {OVERVIEW_STD_H} px "
          f"(factor {OVERVIEW_FACTOR}x)")

    # ─────────────────────────────────────────────────────────────────────
    # Six-group V1 crop composite
    # ─────────────────────────────────────────────────────────────────────
    def render_one_tile(g, ax, *, label_fs, pct_fs, scalebar=False, ref_bounds=None):
        tiff_path = tiff_path_for(g["animal"], g["slice"])
        c1, c2, c3 = load_tiff_three_channel(tiff_path)
        H, W = c1.shape

        cx, cy = (g["manual_crop_10x"] if g["manual_crop_10x"]
                  else find_crop_centre(c2, g["animal"], g["slice"], H, W))
        cx = max(CROP_HALF_10X, min(W - CROP_HALF_10X, cx))
        cy = max(CROP_HALF_10X, min(H - CROP_HALF_10X, cy))
        x0, x1 = cx - CROP_HALF_10X, cx + CROP_HALF_10X
        y0, y1 = cy - CROP_HALF_10X, cy + CROP_HALF_10X

        c1n, c2n, c3n = normalise_three_channel(c1[y0:y1, x0:x1],
                                                 c2[y0:y1, x0:x1],
                                                 c3[y0:y1, x0:x1],
                                                 ref_bounds=ref_bounds)
        rgb = make_rgb(c1n, c2n, c3n)
        crop_h, crop_w = c1n.shape

        ax.set_facecolor("black")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)
        ax.imshow(rgb, interpolation="lanczos", aspect="equal")
        ax.set_xlim(0, crop_w); ax.set_ylim(crop_h, 0)

        ax.text(0.03, 0.96, g["label"], transform=ax.transAxes,
                fontsize=label_fs, fontweight="bold", color="white",
                va="top", ha="left", path_effects=outline(lw=2))
        if g.get("ic_pct_label"):
            ax.text(0.03, 0.04, g["ic_pct_label"], transform=ax.transAxes,
                    fontsize=pct_fs, fontweight="bold", color="white",
                    va="bottom", ha="left", path_effects=outline(lw=2))
        if scalebar:
            add_scalebar(ax, crop_w, crop_h, scale_um=SCALEBAR_UM,
                         px_per_um=PX_PER_UM_10X, margin_y=0.08, lw=5, fs=14)
        return rgb, crop_h, crop_w

    def fig_3_1_composite():
        rendered = []
        for g in GROUPS:
            fig_t, ax_t = plt.subplots(figsize=(6, 6), dpi=MICROSCOPY_DPI)
            fig_t.patch.set_facecolor("black")
            rgb, crop_h, crop_w = render_one_tile(
                g, ax_t, label_fs=TILE_LABEL_FS, pct_fs=TILE_PCT_FS,
                scalebar=False, ref_bounds=REF_BOUNDS,
            )
            plt.close(fig_t)
            rendered.append((g, rgb, crop_h, crop_w))

        fig, axes = plt.subplots(2, 3,
                                  figsize=(COMPOSITE_W_IN, COMPOSITE_W_IN * 2 / 3),
                                  dpi=MICROSCOPY_DPI)
        fig.patch.set_facecolor("black")
        plt.subplots_adjust(wspace=0.02, hspace=0.02,
                            left=0.005, right=0.995, top=0.995, bottom=0.005)

        for i, (ax, (g, rgb, crop_h, crop_w)) in enumerate(zip(axes.flat, rendered)):
            ax.set_facecolor("black")
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values(): s.set_visible(False)
            ax.imshow(rgb, interpolation="lanczos", aspect="equal")
            ax.set_xlim(0, crop_w); ax.set_ylim(crop_h, 0)
            ax.text(0.50, 0.96, g["label"], transform=ax.transAxes,
                    fontsize=COMP_LABEL_FS, fontweight="bold", color="white",
                    va="top", ha="center", path_effects=outline(lw=2))
            if g.get("ic_pct_label"):
                ax.text(0.03, 0.04, g["ic_pct_label"], transform=ax.transAxes,
                        fontsize=COMP_PCT_FS, fontweight="bold", color="white",
                        va="bottom", ha="left", path_effects=outline(lw=2))
            if i == 5:
                add_scalebar(ax, crop_w, crop_h, scale_um=SCALEBAR_UM,
                             px_per_um=PX_PER_UM_10X, margin_y=0.08, lw=5, fs=14)
        return fig, rendered

    print("\nGenerating fig_3_1 (six-group composite) ...")
    _fig, _rendered_tiles = fig_3_1_composite()
    _fig.savefig(f"{MICROSCOPY_OUT}/fig_3_1_composite.png",
                 dpi=MICROSCOPY_DPI, bbox_inches="tight",
                 facecolor="black", pad_inches=0.05)
    plt.close(_fig)

    # ─────────────────────────────────────────────────────────────────────
    # mScarlet-only greyscale panel
    # ─────────────────────────────────────────────────────────────────────
    def fig_3_1_mscarlet_supp():
        fig, axes = plt.subplots(2, 3,
                                  figsize=(COMPOSITE_W_IN, COMPOSITE_W_IN * 2 / 3),
                                  dpi=MICROSCOPY_DPI)
        fig.patch.set_facecolor("black")
        plt.subplots_adjust(wspace=0.02, hspace=0.02,
                            left=0.005, right=0.995, top=0.995, bottom=0.005)

        for i, (ax, (g, _rgb, crop_h, crop_w)) in enumerate(zip(axes.flat, _rendered_tiles)):
            c1, c2, c3 = load_tiff_three_channel(tiff_path_for(g["animal"], g["slice"]))
            H, W = c1.shape
            cx, cy = (g["manual_crop_10x"] if g["manual_crop_10x"]
                      else find_crop_centre(c2, g["animal"], g["slice"], H, W))
            cx = max(CROP_HALF_10X, min(W - CROP_HALF_10X, cx))
            cy = max(CROP_HALF_10X, min(H - CROP_HALF_10X, cy))
            x0, x1 = cx - CROP_HALF_10X, cx + CROP_HALF_10X
            y0, y1 = cy - CROP_HALF_10X, cy + CROP_HALF_10X
            c3n = norm_gamma(c3[y0:y1, x0:x1], NORM_LO_MS, NORM_HI_MS, GAMMA_MS)

            ax.set_facecolor("black")
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values(): s.set_visible(False)
            ax.imshow(c3n, cmap="gray", vmin=0, vmax=1,
                      interpolation="lanczos", aspect="equal")
            ax.set_xlim(0, crop_w); ax.set_ylim(crop_h, 0)
            ax.text(0.50, 0.96, g["label"], transform=ax.transAxes,
                    fontsize=COMP_LABEL_FS, fontweight="bold", color="white",
                    va="top", ha="center", path_effects=outline(lw=2))
            if i == 5:
                add_scalebar(ax, crop_w, crop_h, scale_um=SCALEBAR_UM,
                             px_per_um=PX_PER_UM_10X, margin_y=0.08, lw=5, fs=14)
        return fig

    print("Generating fig_3_1 mScarlet supplement ...")
    fig_3_1_mscarlet_supp().savefig(
        f"{MICROSCOPY_OUT}/fig_3_1_mscarlet_supp.png",
        dpi=MICROSCOPY_DPI, bbox_inches="tight",
        facecolor="black", pad_inches=0.05,
    )
    plt.close("all")

    # ─────────────────────────────────────────────────────────────────────
    # Whole-section overview composite (V1 fill, atlas, crop box)
    # ─────────────────────────────────────────────────────────────────────
    def render_overview_panel(ax, g, ref_bounds, *, label_pos="centre",
                              show_scalebar=False):
        """Render a single whole-section overview panel onto `ax`."""
        c1, c2, c3 = load_tiff_three_channel(tiff_path_for(g["animal"], g["slice"]))
        H, W = c1.shape

        cx, cy = (g["manual_crop_10x"] if g["manual_crop_10x"]
                  else find_crop_centre(c2, g["animal"], g["slice"], H, W))
        cx = max(CROP_HALF_10X, min(W - CROP_HALF_10X, cx))
        cy = max(CROP_HALF_10X, min(H - CROP_HALF_10X, cy))

        c1 = block_downsample(c1, OVERVIEW_FACTOR)
        c2 = block_downsample(c2, OVERVIEW_FACTOR)
        c3 = block_downsample(c3, OVERVIEW_FACTOR)
        Hds, Wds = c1.shape

        c1n, c2n, c3n = normalise_three_channel(c1, c2, c3, ref_bounds=ref_bounds)
        rgb = make_rgb(c1n, c2n, c3n)

        if SHOW_VISUAL_FILL or SHOW_ATLAS_OVERLAY:
            atlas_path = atlas_path_for(g["animal"], g["slice"])
            if atlas_path.exists():
                atlas_ds = resize_atlas_nearest(load_atlas_rgb(atlas_path), Wds, Hds)
                if SHOW_VISUAL_FILL:
                    v1_mask = make_region_mask(atlas_ds, V1_COLORS)
                    for c in range(3):
                        rgb[..., c] = np.where(
                            v1_mask,
                            V1_FILL_ALPHA * V1_FILL_COLOR[c]
                            + (1 - V1_FILL_ALPHA) * rgb[..., c],
                            rgb[..., c],
                        )
                if SHOW_ATLAS_OVERLAY:
                    rgb = composite_edges(
                        rgb,
                        detect_region_boundaries(atlas_ds,
                                                  thickness=ATLAS_EDGE_THICKNESS),
                        ATLAS_EDGE_COLOR, ATLAS_EDGE_ALPHA,
                    )

        ax.set_facecolor("black")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)
        ax.imshow(rgb, interpolation="lanczos", aspect="equal")

        if OVERVIEW_ZOOM:
            _zw    = int(OVERVIEW_STD_W * ZOOM_FRAC / 2)
            _zh    = int(OVERVIEW_STD_H * ZOOM_FRAC / 2)
            _cx_ds = cx / OVERVIEW_FACTOR
            _cy_ds = cy / OVERVIEW_FACTOR
            ax.set_xlim(_cx_ds - _zw, _cx_ds + _zw)
            ax.set_ylim(_cy_ds + _zh, _cy_ds - _zh)
        else:
            dx = (OVERVIEW_STD_W - Wds) // 2
            dy = (OVERVIEW_STD_H - Hds) // 2
            ax.set_xlim(-dx, OVERVIEW_STD_W - dx)
            ax.set_ylim(OVERVIEW_STD_H - dy, -dy)

        box_x0   = (cx - CROP_HALF_10X) / OVERVIEW_FACTOR
        box_y0   = (cy - CROP_HALF_10X) / OVERVIEW_FACTOR
        box_side = (2 * CROP_HALF_10X)  / OVERVIEW_FACTOR
        rect = plt.Rectangle((box_x0, box_y0), box_side, box_side,
                             fill=False, edgecolor=CROP_BOX_COLOR,
                             linewidth=CROP_BOX_LW)
        rect.set_path_effects(outline(color="black", lw=CROP_BOX_LW + 2))
        ax.add_patch(rect)

        ha = "center" if label_pos == "centre" else "left"
        xt = 0.50    if label_pos == "centre" else 0.03
        ax.text(xt, 0.96, g["label"], transform=ax.transAxes,
                fontsize=COMP_LABEL_FS, fontweight="bold", color="white",
                va="top", ha=ha, path_effects=outline(lw=2))
        if g.get("ic_pct_label"):
            ax.text(0.001, 0.04, g["ic_pct_label"], transform=ax.transAxes,
                    fontsize=COMP_PCT_FS, fontweight="bold", color="white",
                    va="bottom", ha="left", path_effects=outline(lw=2))
        if show_scalebar:
            _bar_um   = ZOOM_SCALEBAR_UM if OVERVIEW_ZOOM else OVERVIEW_SCALEBAR_UM
            _bar_span = OVERVIEW_STD_W * ZOOM_FRAC if OVERVIEW_ZOOM else OVERVIEW_STD_W
            _bfrac    = _bar_um * PX_PER_UM_10X / OVERVIEW_FACTOR / _bar_span
            _xr, _yb = 0.97, 0.03
            _xl = _xr - _bfrac
            ax.plot([_xl, _xr], [_yb, _yb], transform=ax.transAxes,
                    color="white", linewidth=5, solid_capstyle="butt",
                    path_effects=outline(color="black", lw=7))
            ax.text((_xl + _xr) / 2, _yb + 0.012, f"{_bar_um:g} µm",
                    transform=ax.transAxes,
                    color="white", fontsize=14, ha="center", va="bottom",
                    path_effects=outline(color="black", lw=2))

    def fig_D5_overview():
        _cell_h = (COMPOSITE_W_IN / 3) * (OVERVIEW_STD_H / OVERVIEW_STD_W)
        fig, axes = plt.subplots(2, 3, figsize=(COMPOSITE_W_IN, _cell_h * 2),
                                  dpi=MICROSCOPY_DPI)
        fig.patch.set_facecolor("black")
        plt.subplots_adjust(wspace=0.02, hspace=0.02,
                            left=0.005, right=0.995, top=0.995, bottom=0.005)
        for i, (ax, g) in enumerate(zip(axes.flat, GROUPS)):
            render_overview_panel(ax, g, REF_BOUNDS, label_pos="centre",
                                   show_scalebar=(i == 5))
        return fig

    print("Generating fig_D5 (whole-section overview) ...")
    fig_D5_overview().savefig(f"{MICROSCOPY_OUT}/fig_D5_overview.png",
                              dpi=MICROSCOPY_DPI, bbox_inches="tight",
                              facecolor="black", pad_inches=0.05)
    plt.close("all")

    # ─────────────────────────────────────────────────────────────────────
    # Narrative paired figures (overview row + V1 crop row)
    # ─────────────────────────────────────────────────────────────────────
    _label_to_g = {g["label"]: g for g in GROUPS}

    def render_paired_figure(group_labels, out_path):
        n      = len(group_labels)
        groups = [_label_to_g[lbl] for lbl in group_labels]

        _cell_w = COMPOSITE_W_IN / n
        _hr_ov  = OVERVIEW_STD_H / OVERVIEW_STD_W
        _fig_h  = _cell_w * (_hr_ov + 1.0)

        fig = plt.figure(figsize=(COMPOSITE_W_IN, _fig_h), dpi=MICROSCOPY_DPI)
        fig.patch.set_facecolor("black")
        gs = fig.add_gridspec(
            2, n, height_ratios=[_hr_ov, 1.0],
            wspace=0.02, hspace=0.04,
            left=0.005, right=0.995, top=0.995, bottom=0.005,
        )

        for col, g in enumerate(groups):
            ax_ov = fig.add_subplot(gs[0, col])
            ax_cr = fig.add_subplot(gs[1, col])

            render_overview_panel(
                ax_ov, g, REF_BOUNDS,
                label_pos="centre",
                show_scalebar=(col == n - 1),
            )

            # V1 crop row
            c1, c2, c3 = load_tiff_three_channel(tiff_path_for(g["animal"], g["slice"]))
            H, W = c1.shape
            cx, cy = (g["manual_crop_10x"] if g["manual_crop_10x"]
                      else find_crop_centre(c2, g["animal"], g["slice"], H, W))
            cx = max(CROP_HALF_10X, min(W - CROP_HALF_10X, cx))
            cy = max(CROP_HALF_10X, min(H - CROP_HALF_10X, cy))
            x0, x1 = cx - CROP_HALF_10X, cx + CROP_HALF_10X
            y0, y1 = cy - CROP_HALF_10X, cy + CROP_HALF_10X

            c1n, c2n, c3n = normalise_three_channel(
                c1[y0:y1, x0:x1], c2[y0:y1, x0:x1], c3[y0:y1, x0:x1],
                ref_bounds=REF_BOUNDS,
            )
            rgb_cr = make_rgb(c1n, c2n, c3n)
            crop_h, crop_w = c1n.shape

            ax_cr.set_facecolor("black")
            ax_cr.set_xticks([]); ax_cr.set_yticks([])
            for s in ax_cr.spines.values(): s.set_visible(False)
            ax_cr.imshow(rgb_cr, interpolation="lanczos", aspect="equal")
            ax_cr.set_xlim(0, crop_w); ax_cr.set_ylim(crop_h, 0)

            if col == n - 1:
                _bfrac_cr      = SCALEBAR_UM * PX_PER_UM_10X / crop_w
                _xr_cr, _yb_cr = 0.97, 0.03
                _xl_cr         = _xr_cr - _bfrac_cr
                ax_cr.plot([_xl_cr, _xr_cr], [_yb_cr, _yb_cr],
                           transform=ax_cr.transAxes,
                           color="white", linewidth=5, solid_capstyle="butt",
                           path_effects=outline(color="black", lw=7))
                ax_cr.text((_xl_cr + _xr_cr) / 2, _yb_cr + 0.022,
                           f"{SCALEBAR_UM:g} µm",
                           transform=ax_cr.transAxes,
                           color="white", fontsize=14, ha="center", va="bottom",
                           path_effects=outline(color="black", lw=2))

        fig.savefig(str(out_path), dpi=MICROSCOPY_DPI, bbox_inches="tight",
                    facecolor="black", pad_inches=0.05)
        return fig

    PAIRED_FIGURES = [
        (["Control", "ADAMTS-15", "C6ST-1 / ADAMTS-15"], "fig_3_1_primary_findings.png"),
        (["Control", "C6ST-1",    "C6ST-1 / ADAMTS-15"], "fig_3_4_c6st1_interaction.png"),
        (["Control", "ADAMTS-4",  "ADAMTS-4 + MD"],      "fig_3_5_adamts4_md.png"),
    ]

    for _labels, _fname in PAIRED_FIGURES:
        print(f"Generating {_fname} ...")
        _f = render_paired_figure(_labels, MICROSCOPY_OUT / _fname)
        plt.close(_f)

    print("\nMicroscopy figures complete.")

# =============================================================================
# Summary
# =============================================================================
print("\n" + "=" * 60)
print("All figures complete.")
print("=" * 60)
print(f"\nOutput directory: {OUTPUT_DIR}")
if _microscopy_available:
    print(f"Microscopy outputs: {MICROSCOPY_OUT}")
print("\nRepresentative microscopy figures (§3.1, §3.4, §3.5, §D.5):")
print("  fig_3_1_composite.png             — six-group V1 crop composite")
print("  fig_3_1_mscarlet_supp.png         — mScarlet-only greyscale (supp)")
print("  fig_D5_overview.png               — whole-section overview composite")
print("  fig_3_1_primary_findings.png      — overview + crop pairs (§3.1)")
print("  fig_3_4_c6st1_interaction.png     — overview + crop pairs (§3.4)")
print("  fig_3_5_adamts4_md.png            — overview + crop pairs (§3.5)")
print("\nMain figures:")
print("  fig_3_2_1.png         — Primary LMM forest plot")
print("  fig_3_2_2.png         — Animal-level enwrapment ratios")
print("  fig_3_2_3.png         — Layer-stratified forest plots")
print("  fig_3_3_1.png         — Distance-decay curves")
print("  fig_3_4_1.png         — Factorial interaction")
print("  fig_3_5_1.png         — Zone profile ADAMTS-4 vs ADAMTS-4-MD")
print("  fig_3_6_1.png         — Zone-stratified WFA and PV density")
print("  fig_3_7_1.png         — Virus+/virus- enwrapment (Bayesian)")
print("  fig_3_7_2.png         — Distance-bin enwrapment")
print("  fig_3_8.png           — Virus-status enwrapment comparison")
print("  fig_viral_load_lines.png — Viral load dose-response curves")
print("\nSupplementary figures:")
print("  fig_3_2_s1.png        — WFA cell density (hemisphere LMM)")
print("  fig_3_3_s1.png        — Off-target specificity (Bayesian)")
print("  fig_3_7_s1.png        — Viral load dose-response slopes")
print("\nAppendix D figures:")
print("  D1_hemisphere_fluorescence.png — Hemisphere WFA/PV (preliminary)")
print("  D2_zone_fluorescence.png       — Zone virus+/virus- fluor. (preliminary)")
print("  D3_power_analysis.png          — Cohen's d power analysis")