#!/usr/bin/env python3
"""
fig_virus_proximal_40x.py
==========================
Virus-proximal / virus-distal classification figure using high-resolution
(40×) warped TIFFs and injection marker (C3) cell detections.

Classification:
    A PV cell is "virus-proximal" if its nearest detected injection-marker
    cell (C3) lies within COLOC_PX_10X pixels in 10× space.  Cells beyond
    that threshold are "virus-distal".  The threshold and crop dimensions
    are scaled from 10× to 40× using WARP_SCALE.

Outputs (saved to OUT_DIR):
    fullslide_proximal_<animal>_<slice>.png  — 10× overview with crop box
    crop_proximal_<animal>_<slice>_40x.png  — 40× crop with PV overlay
    schematic_proximal.png                  — two-panel classification schematic

Crop selection:
    Each candidate window is scored by min(n_proximal, n_distal) × focus,
    maximising representation of both PV populations in the same frame.
    Override auto-selection by setting MANUAL_CX_40X / MANUAL_CY_40X.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import tifffile
from pathlib import Path
from scipy.spatial import cKDTree
from scipy.ndimage import laplace
import pandas as pd
from skimage.transform import resize as sk_resize

# ==============================================================================
# CONFIGURATION
# ==============================================================================

ANIMAL    = "ANIMAL_ID"
SLICE_ID  = "s000"
TREATMENT = "GROUP_NAME"

# Override auto crop selection (set after an initial run, or leave as None)
MANUAL_CX_40X = None
MANUAL_CY_40X = None

# ── Paths ─────────────────────────────────────────────────────────────────────
TIFF_ROOT = Path("/path/to/originals")   # Root containing <group>/<animal>/ subfolders
TIFF_10X  = TIFF_ROOT / TREATMENT / ANIMAL / f"{ANIMAL}_{SLICE_ID}.tiff"
TIFF_40X  = Path(f"/path/to/40x/{ANIMAL}_{SLICE_ID}_40x_warped.tiff")

# C3 injection-marker detection CSV (10× pixel coordinates)
C3_CSV    = Path(f"/path/to/analysis_results/Counts_C3-{ANIMAL}_{SLICE_ID}.csv")

CELLS_CSV = Path("/path/to/cells_with_zones.csv")
OUT_DIR   = Path("/path/to/output/proximal")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Channel indices (0-based) ──────────────────────────────────────────────────
CH_PV = 1   # Interneuron marker
CH_MS = 2   # Injection marker (mScarlet or equivalent)

# ── Imaging parameters ────────────────────────────────────────────────────────
WARP_SCALE    = 4
UM_PER_PX_10X = 0.65
UM_PER_PX_40X = 0.1625
PX_PER_UM_10X = 1 / UM_PER_PX_10X
PX_PER_UM_40X = 1 / UM_PER_PX_40X

# Classification threshold and distance bins (10× pixels; scaled to 40× below)
COLOC_PX_10X  = 30
DIST_BINS_10X = [30, 100, 300, 600]

COLOC_PX_40X  = COLOC_PX_10X * WARP_SCALE
DIST_BINS_40X = [b * WARP_SCALE for b in DIST_BINS_10X]

CROP_HALF_10X = 500
CROP_HALF_40X = CROP_HALF_10X * WARP_SCALE

# ── Crop selection ─────────────────────────────────────────────────────────────
MIN_PV_IN_WINDOW = 8   # Minimum PV cells required in a candidate window

# ── Output ─────────────────────────────────────────────────────────────────────
DPI   = 300
MAX_W = 2000

# ── Brain areas used for spatial filtering ────────────────────────────────────
TARGET_BRAIN_AREAS = [
    "Primary visual area, layer 1",
    "Primary visual area, layer 2/3",
    "Primary visual area, layer 4",
    "Primary visual area, layer 5",
    "Primary visual area, layer 6a",
    "Primary visual area, layer 6b",
    "Anterolateral visual area",
    "Posteromedial visual area",
    "Laterointermediate area",
    "Anteromedial visual area",
]

# ── Colours ───────────────────────────────────────────────────────────────────
COL_PV   = "#FF44FF"   # Interneuron marker channel (magenta)
COL_PROX = "#FFE800"   # Virus-proximal PV cells (yellow)
COL_DIST = "#00FFFF"   # Virus-distal PV cells (cyan)
COL_MS   = "#74C476"   # Injection-marker+ cells (green)

# ── Scale bar ─────────────────────────────────────────────────────────────────
SCALEBAR_LW = 10
SCALEBAR_FS = 32

# ==============================================================================
# HELPERS
# ==============================================================================

def norm(img, lo_pct=0.5, hi_pct=99.8):
    """Percentile-stretch a 2-D array to [0, 1]."""
    lo = np.percentile(img, lo_pct)
    hi = np.percentile(img, hi_pct)
    if hi == lo:
        return np.zeros_like(img, dtype=float)
    return np.clip((img.astype(float) - lo) / (hi - lo), 0, 1)


def outline(color="black", lw=2):
    """Path effect: filled stroke behind text for contrast."""
    return [pe.withStroke(linewidth=lw, foreground=color)]


def add_scalebar(ax, img_w, img_h, scale_um, px_per_um,
                 color="white", margin_y=0.06):
    """Draw a scale bar with label using the global SCALEBAR_LW / _FS settings."""
    bar_px = px_per_um * scale_um
    x0 = img_w * 0.04
    x1 = x0 + bar_px
    y  = img_h - img_h * margin_y
    ax.plot([x0, x1], [y, y], color=color, lw=SCALEBAR_LW,
            solid_capstyle="butt", transform=ax.transData)
    ax.text((x0 + x1) / 2, y - img_h * 0.025, f"{scale_um} µm",
            color=color, ha="center", va="bottom",
            fontsize=SCALEBAR_FS, fontweight="bold",
            transform=ax.transData, path_effects=outline())


def focus_score(ch, cx, cy, half, H, W):
    """Laplacian variance of a square patch — higher means sharper focus."""
    x0 = max(0, cx - half); x1 = min(W, cx + half)
    y0 = max(0, cy - half); y1 = min(H, cy + half)
    patch = ch[y0:y1, x0:x1].astype(float)
    return float(laplace(patch).var()) if patch.size > 0 else 0.0


def make_rgb_pv_ms(c_pv, c_ms):
    """RGB composite: PV → magenta, injection marker → orange-red shift."""
    H, W = c_pv.shape
    rgb = np.zeros((H, W, 3), dtype=float)
    rgb[:, :, 0] = np.clip(norm(c_pv) + norm(c_ms), 0, 1)
    rgb[:, :, 1] = np.clip(norm(c_pv) * 0.27,        0, 1)
    rgb[:, :, 2] = np.clip(norm(c_pv),                0, 1)
    return rgb


# ==============================================================================
# LOAD DATA
# ==============================================================================

print("Loading 10× TIFF for full-slide render ...")
raw_10x = tifffile.imread(str(TIFF_10X))
if raw_10x.ndim == 3 and raw_10x.shape[0] <= 5:
    c_pv_10x = raw_10x[CH_PV]; c_ms_10x = raw_10x[CH_MS]
else:
    c_pv_10x = raw_10x[:, :, CH_PV]; c_ms_10x = raw_10x[:, :, CH_MS]
H10, W10 = c_pv_10x.shape
print(f"  10×: {H10}×{W10}")

print("Loading 40× TIFF for crop ...")
raw_40x = tifffile.imread(str(TIFF_40X))
if raw_40x.ndim == 3 and raw_40x.shape[0] <= 5:
    c_pv_40x = raw_40x[CH_PV]; c_ms_40x = raw_40x[CH_MS]
else:
    c_pv_40x = raw_40x[:, :, CH_PV]; c_ms_40x = raw_40x[:, :, CH_MS]
H40, W40 = c_pv_40x.shape
print(f"  40×: {H40}×{W40}")

print("Loading injection-marker detections ...")
c3_df = pd.read_csv(str(C3_CSV)).rename(columns={"Global_X": "x", "Global_Y": "y"})
print(f"  {len(c3_df)} detections (10× pixel space)")
c3_df["x_40x"] = c3_df["x"] * WARP_SCALE
c3_df["y_40x"] = c3_df["y"] * WARP_SCALE

print("Loading 10× PV cells for classification ...")
cols = ["mouse_id", "slice_id", "cell_type", "x_hires", "y_hires",
        "zone", "brain_area"]
cells = pd.read_csv(str(CELLS_CSV), usecols=cols, low_memory=False)
pv_10x = cells[
    (cells["mouse_id"] == ANIMAL) &
    (cells["slice_id"] == SLICE_ID) &
    (cells["cell_type"] == "PV")
].copy()
print(f"  {len(pv_10x)} PV cells in 10× space")

# Classify proximal / distal by nearest-neighbour distance to C3 detections
tree_c3       = cKDTree(c3_df[["x", "y"]].values)
dists, _      = tree_c3.query(pv_10x[["x_hires", "y_hires"]].values, k=1)
pv_10x["dist_to_c3"] = dists
pv_10x["proximal"]   = dists <= COLOC_PX_10X

# Scale PV coordinates to 40× warped space
pv_10x["x_40x"] = pv_10x["x_hires"] * WARP_SCALE
pv_10x["y_40x"] = pv_10x["y_hires"] * WARP_SCALE

# Restrict to target brain areas
pv_visp = pv_10x[pv_10x["brain_area"].isin(TARGET_BRAIN_AREAS)].copy()
if len(pv_visp) == 0:
    print("  [WARN] No target-area cells — using all PV cells")
    pv_visp = pv_10x.copy()

print(f"  Target-area PV: {len(pv_visp)}  "
      f"proximal={pv_visp['proximal'].sum()}  "
      f"distal={(~pv_visp['proximal']).sum()}")

# ==============================================================================
# CROP CENTRE SELECTION
# ==============================================================================

if MANUAL_CX_40X is not None and MANUAL_CY_40X is not None:
    cx_40x, cy_40x = MANUAL_CX_40X, MANUAL_CY_40X
    print(f"\nManual crop centre: ({cx_40x}, {cy_40x})")
else:
    print("\nAuto-selecting crop centre ...")
    margin = CROP_HALF_40X + 20
    candidates = pv_visp[
        pv_visp["proximal"] &
        (pv_visp["x_40x"] > margin) & (pv_visp["x_40x"] < W40 - margin) &
        (pv_visp["y_40x"] > margin) & (pv_visp["y_40x"] < H40 - margin)
    ]

    if len(candidates) == 0:
        raise RuntimeError("No valid proximal PV candidates in target brain areas.")

    best_cx, best_cy, best_score = 0, 0, -1

    for _, row in candidates.iterrows():
        ccx, ccy = int(row["x_40x"]), int(row["y_40x"])
        in_win = (
            (pv_visp["x_40x"] >= ccx - CROP_HALF_40X) &
            (pv_visp["x_40x"] <  ccx + CROP_HALF_40X) &
            (pv_visp["y_40x"] >= ccy - CROP_HALF_40X) &
            (pv_visp["y_40x"] <  ccy + CROP_HALF_40X)
        )
        n_prox = int(pv_visp.loc[in_win, "proximal"].sum())
        n_dist = int((~pv_visp.loc[in_win, "proximal"]).sum())
        if (n_prox + n_dist) < MIN_PV_IN_WINDOW:
            continue
        # Score: maximise min(proximal, distal) so both populations appear
        score = min(n_prox, n_dist) * focus_score(
            c_pv_40x, ccx, ccy, CROP_HALF_40X, H40, W40
        )
        if score > best_score:
            best_score = score
            best_cx, best_cy = ccx, ccy

    cx_40x, cy_40x = best_cx, best_cy
    cx_10x = cx_40x // WARP_SCALE
    cy_10x = cy_40x // WARP_SCALE
    print(f"  Selected 40× crop centre: ({cx_40x}, {cy_40x})  "
          f"→ 10×: ({cx_10x}, {cy_10x})")

# Clamp to image bounds
cx_40x = max(CROP_HALF_40X, min(W40 - CROP_HALF_40X, cx_40x))
cy_40x = max(CROP_HALF_40X, min(H40 - CROP_HALF_40X, cy_40x))
cx_10x = cx_40x // WARP_SCALE
cy_10x = cy_40x // WARP_SCALE

# ==============================================================================
# PNG 1 — FULL SLIDE (10×)
# ==============================================================================

print("\nRendering full slide ...")
rgb_full = make_rgb_pv_ms(c_pv_10x, c_ms_10x)

scale  = MAX_W / W10
fit_w, fit_h = MAX_W, int(H10 * scale)
canvas = sk_resize(rgb_full, (fit_h, fit_w, 3), order=1,
                   preserve_range=True).astype(float)

fig_w = 12
fig   = plt.figure(figsize=(fig_w, round(fig_w * fit_h / fit_w, 2)), dpi=DPI)
ax    = fig.add_axes([0, 0, 1, 1])
fig.patch.set_facecolor("black")
ax.set_facecolor("black"); ax.set_axis_off()

ax.imshow(np.clip(canvas, 0, 1), interpolation="lanczos", aspect="auto")
ax.set_xlim(0, fit_w); ax.set_ylim(fit_h, 0)

rx0 = int((cx_10x - CROP_HALF_10X) * scale)
rx1 = int((cx_10x + CROP_HALF_10X) * scale)
ry0 = int((cy_10x - CROP_HALF_10X) * scale)
ry1 = int((cy_10x + CROP_HALF_10X) * scale)
ax.add_patch(mpatches.Rectangle(
    (rx0, ry0), rx1 - rx0, ry1 - ry0,
    fill=False, edgecolor="white", linewidth=1.2, zorder=5
))
add_scalebar(ax, fit_w, fit_h, scale_um=1000,
             px_per_um=PX_PER_UM_10X * scale)
ax.text(0.97, 0.97, "Injection marker", color=COL_MS, ha="right", va="top",
        fontsize=11, fontweight="bold", transform=ax.transAxes,
        path_effects=outline())
ax.text(0.97, 0.91, "PV", color=COL_PV, ha="right", va="top",
        fontsize=11, fontweight="bold", transform=ax.transAxes,
        path_effects=outline())
ax.text(0.03, 0.03, f"{ANIMAL} · {SLICE_ID}",
        color="white", ha="left", va="bottom", fontsize=9,
        transform=ax.transAxes, path_effects=outline())

out_full = OUT_DIR / f"fullslide_proximal_{ANIMAL}_{SLICE_ID}.png"
fig.savefig(str(out_full), dpi=DPI, bbox_inches=None, pad_inches=0,
            facecolor=fig.get_facecolor())
plt.close(fig)
print(f"  Saved: {out_full.name}")

# ==============================================================================
# PNG 2 — CROP (40×)
# ==============================================================================

print("Rendering 40× crop ...")

x0_40 = max(0, cx_40x - CROP_HALF_40X)
x1_40 = min(W40, cx_40x + CROP_HALF_40X)
y0_40 = max(0, cy_40x - CROP_HALF_40X)
y1_40 = min(H40, cy_40x + CROP_HALF_40X)

crop_pv  = c_pv_40x[y0_40:y1_40, x0_40:x1_40]
crop_ms  = c_ms_40x[y0_40:y1_40, x0_40:x1_40]
rgb_crop = make_rgb_pv_ms(crop_pv, crop_ms)
crop_h, crop_w = crop_pv.shape

# Filter cells and markers to the crop window
pv_in = pv_visp[
    (pv_visp["x_40x"] >= x0_40) & (pv_visp["x_40x"] < x1_40) &
    (pv_visp["y_40x"] >= y0_40) & (pv_visp["y_40x"] < y1_40)
].copy()
pv_in["lx"] = pv_in["x_40x"] - x0_40
pv_in["ly"] = pv_in["y_40x"] - y0_40

c3_in = c3_df[
    (c3_df["x_40x"] >= x0_40) & (c3_df["x_40x"] < x1_40) &
    (c3_df["y_40x"] >= y0_40) & (c3_df["y_40x"] < y1_40)
].copy()
c3_in["lx"] = c3_in["x_40x"] - x0_40
c3_in["ly"] = c3_in["y_40x"] - y0_40

sel_lx = cx_40x - x0_40
sel_ly = cy_40x - y0_40

fig = plt.figure(figsize=(12, 12), dpi=DPI)
ax  = fig.add_axes([0, 0, 1, 1])
fig.patch.set_facecolor("black")
ax.set_facecolor("black"); ax.set_axis_off()

ax.imshow(np.clip(rgb_crop, 0, 1), interpolation="lanczos", aspect="auto")
ax.set_xlim(0, crop_w); ax.set_ylim(crop_h, 0)

# Concentric distance rings
ring_styles  = ["--", "-.", ":",  (0, (3, 1, 1, 1))]
ring_alphas  = [1.0,  0.95, 0.85, 0.75]
ring_colours = ["#FFFFFF", "#CCCCCC", "#AAAAAA", "#888888"]
for i, r in enumerate(DIST_BINS_40X):
    ax.add_patch(mpatches.Circle(
        (sel_lx, sel_ly), r,
        fill=False, edgecolor=ring_colours[i],
        linewidth=2.5, linestyle=ring_styles[i],
        alpha=ring_alphas[i], zorder=2
    ))

# Injection-marker detections
for _, row in c3_in.iterrows():
    ax.plot(row["lx"], row["ly"], "s",
            color=COL_MS, markersize=6,
            markeredgecolor="black", markeredgewidth=0.5,
            zorder=3, alpha=0.85)

# PV cells coloured by proximity classification
for _, row in pv_in.iterrows():
    col = COL_PROX if row["proximal"] else COL_DIST
    ax.plot(row["lx"], row["ly"], "o",
            color=col, markersize=7,
            markeredgecolor="white", markeredgewidth=0.7, zorder=5)

add_scalebar(ax, crop_w, crop_h, scale_um=100,
             px_per_um=PX_PER_UM_40X, margin_y=0.13)
ax.text(0.03, 0.03, f"{ANIMAL} · {SLICE_ID} · 40×",
        color="white", ha="left", va="bottom", fontsize=8,
        transform=ax.transAxes, path_effects=outline())

out_crop = OUT_DIR / f"crop_proximal_{ANIMAL}_{SLICE_ID}_40x.png"
fig.savefig(str(out_crop), dpi=DPI, bbox_inches=None, pad_inches=0,
            facecolor=fig.get_facecolor())
plt.close(fig)
print(f"  Saved: {out_crop.name}")

# ==============================================================================
# SCHEMATIC — two-panel proximity classification diagram
# ==============================================================================

def render_schematic(out_path):
    """
    Two-panel schematic illustrating the virus-proximity classification.

    Left panel  — binary threshold: a proximal PV cell (yellow) with one
                  injection-marker cell inside the 30 px threshold circle and
                  several distant ones outside.
    Right panel — continuous distance binning: a distal PV cell (cyan) at
                  centre, surrounded by concentric annular zones with scattered
                  injection-marker dots placed by collision-free sampling.
    """
    np.random.seed(123)

    fig, axes = plt.subplots(1, 2, figsize=(12, 6), dpi=300)
    fig.patch.set_facecolor("#000000")
    plt.subplots_adjust(wspace=0.05, left=0.02, right=0.98,
                        top=0.98, bottom=0.02)

    PV_X, PV_Y = 150, 150
    PV_R       = 20
    THRESH_R   = 60   # Visual radius representing COLOC_PX_10X

    for ax in axes:
        ax.set_facecolor("#000000")
        ax.set_xlim(0, 300); ax.set_ylim(300, 0)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_aspect("equal")
        for spine in ax.spines.values():
            spine.set_visible(False)

    # ── Left: binary threshold ─────────────────────────────────────────────────
    ax = axes[0]

    ax.add_patch(mpatches.Circle(
        (PV_X, PV_Y), THRESH_R,
        fill=False, edgecolor="white", linewidth=2.5,
        linestyle="--", alpha=0.9, zorder=3
    ))
    ax.add_patch(mpatches.Circle(
        (PV_X, PV_Y), PV_R,
        facecolor=COL_PROX, edgecolor="white", linewidth=1.5, zorder=6
    ))

    # Injection-marker cell inside the threshold
    ms_x, ms_y = PV_X + 42, PV_Y - 18
    ax.add_patch(mpatches.Circle(
        (ms_x, ms_y), 8,
        facecolor=COL_MS, edgecolor="white", linewidth=1.0, zorder=5
    ))

    # Distance annotation arrow
    ax.annotate("", xy=(ms_x, ms_y), xytext=(PV_X, PV_Y),
                arrowprops=dict(arrowstyle="<->", color="white", lw=1.5,
                                shrinkA=0, shrinkB=8),
                zorder=12)
    ax.text(PV_X - 18, (PV_Y + ms_y) / 2 - 6,
            f"≤ {COLOC_PX_10X} px", color="white", fontsize=9,
            ha="right", va="bottom",
            path_effects=outline(lw=3), zorder=13)

    # Distant injection-marker cells outside the threshold
    for ang, dist in [(0.5, 95), (2.1, 110), (4.0, 100), (5.2, 85), (1.2, 130)]:
        ox, oy = PV_X + dist * np.cos(ang), PV_Y + dist * np.sin(ang)
        if 0 < ox < 300 and 0 < oy < 300:
            ax.add_patch(mpatches.Circle(
                (ox, oy), 6,
                facecolor=COL_MS, edgecolor="white", linewidth=0.8,
                alpha=0.45, zorder=4
            ))

    # ── Right: continuous distance binning ────────────────────────────────────
    ax = axes[1]

    ax.add_patch(mpatches.Circle(
        (PV_X, PV_Y), PV_R,
        facecolor=COL_DIST, edgecolor="white", linewidth=1.5, zorder=6
    ))

    bin_radii  = [60, 90, 125, 146]   # Visual ring radii (schematic units)
    bin_styles = ["--", "-.", ":", (0, (3, 1, 1, 1))]
    bin_cols   = ["#FFFFFF", "#AAAAAA", "#888888", "#666666"]

    for i, r in enumerate(bin_radii):
        ax.add_patch(mpatches.Circle(
            (PV_X, PV_Y), r,
            fill=False, edgecolor=bin_cols[i],
            linewidth=1.2, linestyle=bin_styles[i], alpha=0.7, zorder=3
        ))

    # Scatter injection-marker dots per annular bin using collision avoidance
    placed = [(PV_X, PV_Y, PV_R + 5)]   # Seed with the soma as an obstacle

    for i in range(len(bin_radii)):
        r_in  = 0 if i == 0 else bin_radii[i - 1]
        r_out = bin_radii[i]
        n_dots = [7, 11, 15, 12][i]

        placed_in_bin = 0
        attempts      = 0
        while placed_in_bin < n_dots and attempts < 1500:
            attempts += 1
            ang  = np.random.uniform(0, 2 * np.pi)
            dist = np.random.uniform(r_in + 7, r_out - 7)
            mx   = PV_X + dist * np.cos(ang)
            my   = PV_Y + dist * np.sin(ang)
            dot_r = np.random.uniform(3.5, 5.0)

            if any(np.hypot(mx - px, my - py) < dot_r + pr + 4
                   for px, py, pr in placed):
                continue
            if not (8 < mx < 292 and 8 < my < 292):
                continue

            ax.add_patch(mpatches.Circle(
                (mx, my), dot_r,
                facecolor=COL_MS, edgecolor="white", linewidth=0.6,
                alpha=np.random.uniform(0.5, 0.85), zorder=4
            ))
            placed.append((mx, my, dot_r))
            placed_in_bin += 1

    fig.savefig(str(out_path), facecolor=fig.get_facecolor(),
                bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    print(f"  Saved: {out_path.name}")


render_schematic(OUT_DIR / "schematic_proximal.png")

print("\nAll done.")
print(f"To fix the crop centre for future runs, add to the configuration:")
print(f"  MANUAL_CX_40X = {cx_40x}")
print(f"  MANUAL_CY_40X = {cy_40x}")
