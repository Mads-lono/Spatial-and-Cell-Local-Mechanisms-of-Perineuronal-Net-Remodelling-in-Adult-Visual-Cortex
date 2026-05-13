#!/usr/bin/env python3
"""
fig_enwrapment_comparison_40x.py
=================================
Generates enwrapment comparison figures using high-resolution (40×) warped
TIFFs and enwrapment CSVs produced by a CPN detection pipeline.

Outputs per run (saved to OUT_DIR):
    fullslide_<animal>_<slice>.png   — 10× overview with crop-box rectangle
    crop_<animal>_<slice>_40x.png    — 40× crop with PV marker overlay
    schematic_enwrapment.png         — two-panel enwrapment concept schematic

Crop selection:
    Candidates are scored by focus quality (Laplacian variance of the PV
    channel) and nearest-PNN distance.  Set A1/A2_SKIP_TOP_N > 0 to skip the
    top-ranked window(s) and explore alternatives.

Spatial filtering:
    Cells are restricted to the ipsilateral hemisphere and specified brain
    areas using spatial bounds derived from a 10× cells CSV, scaled to 40×
    coordinate space by WARP_SCALE.  Disable by setting USE_SPATIAL_FILTER =
    False to use all detections.

40× images are assumed to be warped into 10× coordinate space at 40×
resolution; the scale relationship is simply warped_xy = 10x_xy × WARP_SCALE.
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

# ── Animal 1 (control) ────────────────────────────────────────────────────────
A1_ANIMAL    = "CONTROL_ANIMAL_ID"
A1_SLICE     = "s000"
A1_TREATMENT = "CONTROL_GROUP"
A1_USE_40X   = True   # Set False to fall back to 10× crop

A1_TIFF_40X  = Path("/path/to/40x/CONTROL_ANIMAL_ID_s000_40x_warped.tiff")
A1_ENW_CSV   = Path("/path/to/40x/CONTROL_ANIMAL_ID_s000_enwrapment.csv")

# Skip top N crop candidates (0 = best-scoring window; increase to explore)
A1_SKIP_TOP_N = 0

# ── Animal 2 (treated) ────────────────────────────────────────────────────────
A2_ANIMAL    = "TREATED_ANIMAL_ID"
A2_SLICE     = "s000"
A2_TREATMENT = "TREATED_GROUP"
A2_USE_40X   = True

A2_TIFF_40X  = Path("/path/to/40x/TREATED_ANIMAL_ID_s000_40x_warped.tiff")
A2_ENW_CSV   = Path("/path/to/40x/TREATED_ANIMAL_ID_s000_enwrapment.csv")

A2_SKIP_TOP_N = 0

# ── Paths ─────────────────────────────────────────────────────────────────────
TIFF_ROOT = Path("/path/to/originals")   # Root containing <group>/<animal>/ subfolders
CELLS_CSV = Path("/path/to/cells_with_zones.csv")
OUT_DIR   = Path("/path/to/output/enwrapment")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Coordinate relationship between 10× and warped 40× space ──────────────────
WARP_SCALE = 4   # warped_xy = 10x_xy × WARP_SCALE

# ── Channel indices (0-based) ──────────────────────────────────────────────────
CH_WFA = 0   # WFA / PNN
CH_PV  = 1   # Interneuron marker

# ── Enwrapment detection parameters ───────────────────────────────────────────
COLOC_PX_10X  = 30    # Enwrapment threshold in 10× pixels
CROP_HALF_10X = 150   # Crop half-width in 10× pixels

UM_PER_PX_40X = 0.1625
# Threshold and crop scaled to the same physical size as their 10× equivalents
COLOC_PX_40X  = int(COLOC_PX_10X  * (0.65 / UM_PER_PX_40X))
CROP_HALF_40X = int(CROP_HALF_10X * (0.65 / UM_PER_PX_40X))

PX_PER_UM_10X = 1 / 0.65
PX_PER_UM_40X = 1 / UM_PER_PX_40X

# ── Spatial filter settings ────────────────────────────────────────────────────
USE_SPATIAL_FILTER      = True
FOCUS_PERCENTILE_THRESHOLD = 50   # Retain top N% of candidates by focus score
MIN_PV_IN_WINDOW           = 5    # Minimum PV cells required in a candidate window

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

# ── Output ────────────────────────────────────────────────────────────────────
DPI   = 300
MAX_W = 2000   # Maximum render width for full-slide images

# ── Colours ───────────────────────────────────────────────────────────────────
COL_WFA           = "#00FF00"   # WFA / PNN channel
COL_PV            = "#FF44FF"   # Interneuron marker channel
COL_ENWRAPPED     = "#00FFFF"   # Enwrapped PV cells (cyan)
COL_NOT_ENWRAPPED = "#FFFFFF"   # Non-enwrapped PV cells (white)

MRK_ENWRAPPED     = "^"   # Upward triangle — enwrapped
MRK_NOT_ENWRAPPED = "o"   # Circle — non-enwrapped
SZ_ENWRAPPED      = 9
SZ_NOT_ENWRAPPED  = 7

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
    """Path effect: filled stroke behind text or lines for contrast."""
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


def load_tiff_channels(path):
    """Load WFA (C1) and PV (C2) channels from a multi-channel TIFF."""
    raw = tifffile.imread(str(path))
    if raw.ndim == 3 and raw.shape[0] <= 5:
        return raw[CH_WFA], raw[CH_PV]
    if raw.ndim == 3 and raw.shape[2] <= 5:
        return raw[:, :, CH_WFA], raw[:, :, CH_PV]
    raise ValueError(f"Unexpected TIFF shape: {raw.shape}")


def map_10x_to_40x(x, y):
    """Convert 10× pixel coordinates to warped 40× space."""
    return int(x * WARP_SCALE), int(y * WARP_SCALE)


def make_rgb(c1n, c2n):
    """Build a two-channel RGB composite: C1 → green, C2 → magenta."""
    H, W = c1n.shape
    rgb = np.zeros((H, W, 3), dtype=float)
    rgb[:, :, 0] = np.clip(c2n, 0, 1)   # R: C2
    rgb[:, :, 1] = np.clip(c1n, 0, 1)   # G: C1
    rgb[:, :, 2] = np.clip(c2n, 0, 1)   # B: C2
    return rgb


def focus_score_40x(c2_40x, cx, cy, crop_half, H, W):
    """
    Estimate focus quality at a given crop centre using Laplacian variance.
    Higher variance indicates sharper edges (better focus).
    """
    x0 = max(0, cx - crop_half); x1 = min(W, cx + crop_half)
    y0 = max(0, cy - crop_half); y1 = min(H, cy + crop_half)
    patch = c2_40x[y0:y1, x0:x1].astype(float)
    return float(laplace(patch).var()) if patch.size > 0 else 0.0


def load_cells_40x(animal, slice_id, enw_csv_path):
    """
    Load a 40× CPN enwrapment CSV and apply spatial filters.

    Spatial bounds are derived from the 10× cells CSV and scaled by
    WARP_SCALE to map into 40× warped coordinate space.  A 10% padding
    is added on each side of the target-area bounding box.

    Returns (pv_df, pnn_df) with columns x, y, enwrapped (and others
    from the source CSV).
    """
    df = pd.read_csv(str(enw_csv_path))
    pv_all  = df[df["cell_type"] == "PV"].copy()
    pnn_all = df[df["cell_type"] == "PNN"].copy()
    print(f"  40× CSV: {len(pv_all)} PV, {len(pnn_all)} PNN")

    if not USE_SPATIAL_FILTER:
        return pv_all, pnn_all

    cols = ["mouse_id", "slice_id", "cell_type",
            "x_hires", "y_hires", "zone", "brain_area", "hemisphere"]
    try:
        cells_10x = pd.read_csv(str(CELLS_CSV), usecols=cols, low_memory=False)
    except Exception as e:
        print(f"  [WARN] Could not load cells CSV: {e} — spatial filter disabled")
        return pv_all, pnn_all

    sl = cells_10x[
        (cells_10x["mouse_id"] == animal) &
        (cells_10x["slice_id"] == slice_id)
    ]

    if len(sl) == 0:
        print(f"  [WARN] No 10× cells for {animal} {slice_id} — spatial filter disabled")
        return pv_all, pnn_all

    ipsi = sl[sl["hemisphere"] == "ipsilateral"]
    if len(ipsi) == 0:
        print("  [WARN] No ipsilateral cells in 10× data — using all hemispheres")
        ipsi = sl

    visp = ipsi[ipsi["brain_area"].isin(TARGET_BRAIN_AREAS)]
    if len(visp) == 0:
        print("  [WARN] No target-area cells in 10× data — using full ipsilateral bounds")
        visp = ipsi

    pad_x = (visp["x_hires"].max() - visp["x_hires"].min()) * 0.10
    pad_y = (visp["y_hires"].max() - visp["y_hires"].min()) * 0.10
    x_lo  = (float(visp["x_hires"].min()) - pad_x) * WARP_SCALE
    x_hi  = (float(visp["x_hires"].max()) + pad_x) * WARP_SCALE
    y_lo  = (float(visp["y_hires"].min()) - pad_y) * WARP_SCALE
    y_hi  = (float(visp["y_hires"].max()) + pad_y) * WARP_SCALE

    print(f"  40× filter bounds (±10% pad): "
          f"x=[{x_lo:.0f},{x_hi:.0f}] y=[{y_lo:.0f},{y_hi:.0f}]")

    def _spatial_filter(df):
        return df[
            (df["x"] >= x_lo) & (df["x"] <= x_hi) &
            (df["y"] >= y_lo) & (df["y"] <= y_hi)
        ].copy()

    pv_filt  = _spatial_filter(pv_all)
    pnn_filt = _spatial_filter(pnn_all)
    print(f"  After filter: {len(pv_filt)} PV "
          f"({len(pv_filt)/max(len(pv_all),1)*100:.0f}%), "
          f"{len(pnn_filt)} PNN "
          f"({len(pnn_filt)/max(len(pnn_all),1)*100:.0f}%)")

    enw = pv_filt["enwrapped"].sum()
    print(f"  Filtered enwrapment: {enw}/{len(pv_filt)} = "
          f"{enw/max(len(pv_filt),1)*100:.1f}%")

    return pv_filt, pnn_filt


def load_cells(animal, slice_id):
    """Load PV and PNN cells from the 10× cells CSV for a given animal/slice."""
    cols  = ["mouse_id", "slice_id", "cell_type",
             "x_hires", "y_hires", "zone", "brain_area"]
    cells = pd.read_csv(str(CELLS_CSV), usecols=cols, low_memory=False)
    df    = cells[
        (cells["mouse_id"] == animal) &
        (cells["slice_id"] == slice_id)
    ].copy()
    if len(df) == 0:
        raise RuntimeError(f"No cells found for {animal}, slice {slice_id}.")
    return df


# ==============================================================================
# RENDER FUNCTIONS
# ==============================================================================

def render_fullslide(c1, c2, cx_10x, cy_10x, crop_half_px,
                     px_per_um, label, slice_id, out_path):
    """
    Full-slide 10× overview with a rectangle marking the high-res crop region.
    The image is downsampled to MAX_W for output.
    """
    H, W   = c1.shape
    rgb    = make_rgb(norm(c1), norm(c2))
    scale  = MAX_W / W
    dW, dH = MAX_W, int(H * scale)
    rgb_ds = sk_resize(rgb, (dH, dW, 3), order=1,
                       preserve_range=True).astype(float)

    rx0 = int((cx_10x - crop_half_px) * scale)
    rx1 = int((cx_10x + crop_half_px) * scale)
    ry0 = int((cy_10x - crop_half_px) * scale)
    ry1 = int((cy_10x + crop_half_px) * scale)

    fig, ax = plt.subplots(figsize=(10, 10 * dH / dW), dpi=DPI)
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.imshow(np.clip(rgb_ds, 0, 1), interpolation="lanczos", aspect="equal")
    ax.set_xlim(0, dW); ax.set_ylim(dH, 0)
    ax.add_patch(mpatches.Rectangle(
        (rx0, ry0), rx1 - rx0, ry1 - ry0,
        fill=False, edgecolor="white", linewidth=1.2, zorder=5
    ))
    add_scalebar(ax, dW, dH, scale_um=1000, px_per_um=px_per_um * scale)
    ax.text(0.97, 0.97, "WFA / PNN", color=COL_WFA, ha="right", va="top",
            fontsize=11, fontweight="bold", transform=ax.transAxes,
            path_effects=outline())
    ax.text(0.97, 0.91, "PV", color=COL_PV, ha="right", va="top",
            fontsize=11, fontweight="bold", transform=ax.transAxes,
            path_effects=outline())
    ax.text(0.03, 0.03, f"{label} · {slice_id}",
            color="white", ha="left", va="bottom", fontsize=9,
            transform=ax.transAxes, path_effects=outline())

    fig.savefig(str(out_path), dpi=DPI, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  Saved: {out_path.name}")


def render_crop(c1, c2, cx, cy, pv_cells_img, pnn_cells_img,
                crop_half, coloc_px, px_per_um,
                label, slice_id, use_40x, out_path):
    """
    High-resolution crop panel with PV marker overlay.

    cx, cy         — crop centre in the image's own pixel space
    pv_cells_img   — DataFrame with columns 'cx', 'cy' in image pixel space
    pnn_cells_img  — same
    coloc_px       — enwrapment threshold in image pixels
    """
    H, W = c1.shape
    x0 = cx - crop_half; x1 = cx + crop_half
    y0 = cy - crop_half; y1 = cy + crop_half

    c1_crop = norm(c1[y0:y1, x0:x1])
    c2_crop = norm(c2[y0:y1, x0:x1])
    rgb     = make_rgb(c1_crop, c2_crop)
    crop_h, crop_w = c1_crop.shape

    # Filter cells to the crop window
    pnn_in = pnn_cells_img[
        (pnn_cells_img["cx"] >= x0) & (pnn_cells_img["cx"] < x1) &
        (pnn_cells_img["cy"] >= y0) & (pnn_cells_img["cy"] < y1)
    ].copy()
    pnn_in["lx"] = pnn_in["cx"] - x0
    pnn_in["ly"] = pnn_in["cy"] - y0

    pv_in = pv_cells_img[
        (pv_cells_img["cx"] >= x0) & (pv_cells_img["cx"] < x1) &
        (pv_cells_img["cy"] >= y0) & (pv_cells_img["cy"] < y1)
    ].copy()
    pv_in["lx"] = pv_in["cx"] - x0
    pv_in["ly"] = pv_in["cy"] - y0

    # Re-compute enwrapment within the crop window via nearest-neighbour distance
    if len(pnn_in) > 0:
        tree              = cKDTree(pnn_in[["lx", "ly"]].values)
        dists, _          = tree.query(pv_in[["lx", "ly"]].values, k=1)
        pv_in             = pv_in.copy()
        pv_in["enwrapped"] = dists <= coloc_px
    else:
        pv_in             = pv_in.copy()
        pv_in["enwrapped"] = False

    fig, ax = plt.subplots(figsize=(6, 6), dpi=DPI)
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.imshow(rgb, interpolation="lanczos", aspect="equal")
    ax.set_xlim(0, crop_w); ax.set_ylim(crop_h, 0)

    # Enwrapment threshold circle centred on the selected PV cell
    ax.add_patch(mpatches.Circle(
        (cx - x0, cy - y0), coloc_px,
        fill=False, edgecolor="white", linewidth=1.5,
        linestyle="--", alpha=0.85
    ))

    for _, row in pv_in.iterrows():
        if row["enwrapped"]:
            ax.plot(row["lx"], row["ly"], MRK_ENWRAPPED,
                    color=COL_ENWRAPPED, markersize=SZ_ENWRAPPED,
                    markeredgecolor="black", markeredgewidth=0.8, zorder=5)
        else:
            ax.plot(row["lx"], row["ly"], MRK_NOT_ENWRAPPED,
                    color=COL_NOT_ENWRAPPED, markersize=SZ_NOT_ENWRAPPED,
                    markeredgecolor="black", markeredgewidth=0.6, zorder=4)

    res_label = "40×" if use_40x else "10×"
    add_scalebar(ax, crop_w, crop_h, scale_um=50, px_per_um=px_per_um,
                 margin_y=0.13)
    ax.text(0.03, 0.03, f"{label} · {slice_id} · {res_label}",
            color="white", ha="left", va="bottom", fontsize=8,
            transform=ax.transAxes, path_effects=outline())

    fig.savefig(str(out_path), dpi=DPI, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  Saved: {out_path.name}")


# ── Schematic helpers ──────────────────────────────────────────────────────────

def _draw_pnn_polygons(ax, x, y, inner_r, outer_r,
                       n_polys=3000, alpha_range=(0.4, 0.8),
                       beta_params=(1.2, 2.0)):
    """
    Scatter small random polygons in an annular region to represent
    extracellular matrix texture (e.g. WFA-stained PNN coat).
    """
    for _ in range(n_polys):
        ang  = np.random.uniform(0, 2 * np.pi)
        dist = inner_r + (np.random.beta(*beta_params) * (outer_r - inner_r))
        px, py = x + dist * np.cos(ang), y + dist * np.sin(ang)

        n_v   = np.random.randint(3, 7)
        base  = np.random.uniform(0.15, 0.5)
        v_ang = np.linspace(0, 2 * np.pi, n_v, endpoint=False)
        v_r   = base * (0.6 + np.random.rand(n_v) * 0.7)
        verts = np.stack([px + v_r * np.cos(v_ang),
                          py + v_r * np.sin(v_ang)], axis=1)
        ax.add_patch(mpatches.Polygon(
            verts, facecolor=COL_WFA, edgecolor="none",
            alpha=np.random.uniform(*alpha_range), zorder=2
        ))


def _draw_tapered_process(ax, x, y, soma_r, thresh_r, angle_rad, n_polys=600):
    """
    Draw a tapered cellular process from the soma edge toward the threshold
    radius, widening at the base and narrowing distally.
    """
    for _ in range(n_polys):
        progress = np.random.uniform(0, 1.1)
        dist     = soma_r + progress * (thresh_r - soma_r + 5)
        sigma    = 0.15 * (1.0 - progress * 0.8)
        ang      = angle_rad + np.random.normal(0, sigma)
        px, py   = x + dist * np.cos(ang), y + dist * np.sin(ang)

        n_v   = np.random.randint(3, 6)
        base  = np.random.uniform(0.1, 0.4)
        v_ang = np.linspace(0, 2 * np.pi, n_v, endpoint=False)
        v_r   = base * (0.7 + np.random.rand(n_v) * 0.7)
        verts = np.stack([px + v_r * np.cos(v_ang),
                          py + v_r * np.sin(v_ang)], axis=1)
        ax.add_patch(mpatches.Polygon(
            verts, facecolor=COL_WFA, edgecolor="none",
            alpha=np.random.uniform(0.4, 0.8), zorder=2
        ))


def render_schematic(out_path):
    """
    Two-panel schematic illustrating enwrapment classification.

    Left panel  — enwrapped:  PNN coat encircles the soma with two
                               tapered processes reaching the threshold.
    Right panel — not enwrapped: diffuse WFA clusters outside the threshold.
    """
    np.random.seed(88)

    fig, axes = plt.subplots(1, 2, figsize=(12, 8), dpi=300)
    fig.patch.set_facecolor("#000000")
    plt.subplots_adjust(wspace=0.02, left=0, right=1, bottom=0, top=1)

    PV_X, PV_Y = 100, 100
    PV_R       = 30
    THRESH_R   = 60

    for i, ax in enumerate(axes):
        ax.set_facecolor("#000000")
        ax.set_xlim(30, 170); ax.set_ylim(170, 30)
        ax.axis("off"); ax.set_aspect("equal")

        is_enwrapped = (i == 0)

        # Threshold ring
        ax.add_patch(mpatches.Circle(
            (PV_X, PV_Y), THRESH_R,
            fill=False, edgecolor="#888888", linewidth=2.5,
            linestyle=(0, (5, 5)), alpha=0.9, zorder=10
        ))

        if is_enwrapped:
            # Dense PNN coat around the soma + two tapered processes
            _draw_pnn_polygons(ax, PV_X, PV_Y, PV_R - 4, PV_R + 12, n_polys=5500)
            _draw_tapered_process(ax, PV_X, PV_Y, PV_R, THRESH_R, angle_rad=0.7)
            _draw_tapered_process(ax, PV_X, PV_Y, PV_R, THRESH_R, angle_rad=3.6)
        else:
            # Diffuse, loosely clustered WFA signal outside the threshold
            for _ in range(12):
                dist_c = np.random.uniform(THRESH_R + 10, 65)
                ang_c  = np.random.uniform(0, 2 * np.pi)
                cx_c   = PV_X + dist_c * np.cos(ang_c)
                cy_c   = PV_Y + dist_c * np.sin(ang_c)
                _draw_pnn_polygons(ax, cx_c, cy_c, 0, 12,
                                   n_polys=150, alpha_range=(0.2, 0.5))

        # PV soma
        ax.add_patch(mpatches.Circle(
            (PV_X, PV_Y), PV_R,
            facecolor=COL_PV, edgecolor="white", linewidth=1.5, zorder=15
        ))

    # Threshold-distance annotation on the enwrapped panel
    axes[0].annotate(
        "", xy=(PV_X + THRESH_R, PV_Y), xytext=(PV_X, PV_Y),
        arrowprops=dict(arrowstyle="<->", color="white", lw=2.5,
                        shrinkA=0, shrinkB=0,
                        path_effects=outline(lw=6)), zorder=25
    )
    axes[0].text(
        PV_X + THRESH_R / 2, PV_Y - 3,
        f"{COLOC_PX_10X} px",
        color="white", ha="center", va="bottom",
        fontsize=14, fontweight="bold",
        path_effects=outline(lw=6), zorder=30
    )

    fig.savefig(str(out_path), facecolor=fig.get_facecolor(),
                bbox_inches=None, pad_inches=0)
    plt.close(fig)
    print(f"  Saved: {out_path.name}")


# ==============================================================================
# ANIMAL 1 (CONTROL)
# ==============================================================================

print(f"\n{'='*60}")
print(f"Animal 1: {A1_ANIMAL}  slice {A1_SLICE}  (40× crop: {A1_USE_40X})")
print(f"{'='*60}")

# Load 10× TIFF for full-slide render and cell-selection geometry
tiff_10x_a1 = TIFF_ROOT / A1_TREATMENT / A1_ANIMAL / f"{A1_ANIMAL}_{A1_SLICE}.tiff"
c1_10x_a1, c2_10x_a1 = load_tiff_channels(tiff_10x_a1)
H1_10x, W1_10x = c1_10x_a1.shape

# Load 40× enwrapment detections with spatial filter
print("  Loading 40× enwrapment CSV ...")
pv_a1_40x, pnn_a1_40x = load_cells_40x(A1_ANIMAL, A1_SLICE, A1_ENW_CSV)
pv_a1_40x  = pv_a1_40x.rename(columns={"x": "cx", "y": "cy"})
pnn_a1_40x = pnn_a1_40x.rename(columns={"x": "cx", "y": "cy"})

# Load 10× cells to determine the crop-box location on the full-slide image
cells_a1    = load_cells(A1_ANIMAL, A1_SLICE)
pv_core_a1  = cells_a1[(cells_a1["cell_type"] == "PV")  & (cells_a1["zone"] == "Core")].reset_index(drop=True)
pnn_core_a1 = cells_a1[(cells_a1["cell_type"] == "PNN") & (cells_a1["zone"] == "Core")].reset_index(drop=True)
tree1    = cKDTree(pnn_core_a1[["x_hires", "y_hires"]].values)
dists1, _ = tree1.query(pv_core_a1[["x_hires", "y_hires"]].values, k=1)
pv_core_a1 = pv_core_a1.copy()
pv_core_a1["enwrapped"] = dists1 <= COLOC_PX_10X
best1_10x  = pv_core_a1[pv_core_a1["enwrapped"]].sort_values(
    dists1[pv_core_a1["enwrapped"].values].tolist()
    if False else "enwrapped"  # placeholder — sort by index
).iloc[0]
cx1_10x = int(pv_core_a1[pv_core_a1["enwrapped"]].iloc[0]["x_hires"])
cy1_10x = int(pv_core_a1[pv_core_a1["enwrapped"]].iloc[0]["y_hires"])

if A1_USE_40X:
    print("  Loading 40× TIFF ...")
    c1_40x_a1, c2_40x_a1 = load_tiff_channels(A1_TIFF_40X)
    H1_40x, W1_40x = c1_40x_a1.shape
    print(f"  40× shape: {H1_40x}×{W1_40x}")

    pv_enw_40x = pv_a1_40x[pv_a1_40x["enwrapped"]].copy()

    if len(pv_enw_40x) == 0:
        print("  [WARN] No enwrapped PV cells in 40× CSV — falling back to 10× coords")
        cx1_40x, cy1_40x = map_10x_to_40x(cx1_10x, cy1_10x)
    else:
        margin = CROP_HALF_40X + 20
        pv_enw_40x = pv_enw_40x[
            (pv_enw_40x["cx"] > margin) & (pv_enw_40x["cx"] < W1_40x - margin) &
            (pv_enw_40x["cy"] > margin) & (pv_enw_40x["cy"] < H1_40x - margin)
        ]

        print(f"  Scoring {len(pv_enw_40x)} enwrapped PV candidates by focus ...")
        scores = [focus_score_40x(c2_40x_a1, int(r["cx"]), int(r["cy"]),
                                   CROP_HALF_40X, H1_40x, W1_40x)
                  for _, r in pv_enw_40x.iterrows()]
        pv_enw_40x = pv_enw_40x.copy()
        pv_enw_40x["focus"] = scores
        thr = np.percentile(scores, 100 - FOCUS_PERCENTILE_THRESHOLD)
        focused = pv_enw_40x[pv_enw_40x["focus"] >= thr]
        if len(focused) == 0:
            focused = pv_enw_40x

        best1   = focused.sort_values("dist_to_nearest_pnn").iloc[A1_SKIP_TOP_N]
        cx1_40x = int(best1["cx"])
        cy1_40x = int(best1["cy"])
        cx1_10x = cx1_40x // WARP_SCALE
        cy1_10x = cy1_40x // WARP_SCALE
        print(f"  Selected 40× crop centre: ({cx1_40x}, {cy1_40x})  "
              f"focus={best1['focus']:.1f}  "
              f"dist_pnn={best1['dist_to_nearest_pnn']:.1f} px")

    cx1_40x = max(CROP_HALF_40X, min(W1_40x - CROP_HALF_40X, cx1_40x))
    cy1_40x = max(CROP_HALF_40X, min(H1_40x - CROP_HALF_40X, cy1_40x))

render_fullslide(
    c1_10x_a1, c2_10x_a1, cx1_10x, cy1_10x, CROP_HALF_10X,
    PX_PER_UM_10X, "Control", A1_SLICE,
    OUT_DIR / f"fullslide_{A1_ANIMAL}_{A1_SLICE}.png"
)

if A1_USE_40X:
    render_crop(
        c1_40x_a1, c2_40x_a1, cx1_40x, cy1_40x,
        pv_a1_40x, pnn_a1_40x,
        CROP_HALF_40X, COLOC_PX_40X, PX_PER_UM_40X,
        "Control", A1_SLICE, True,
        OUT_DIR / f"crop_{A1_ANIMAL}_{A1_SLICE}_40x.png"
    )
else:
    pv_a1_10x  = cells_a1[cells_a1["cell_type"] == "PV"].copy()
    pnn_a1_10x = cells_a1[cells_a1["cell_type"] == "PNN"].copy()
    pv_a1_10x["cx"]  = pv_a1_10x["x_hires"]
    pv_a1_10x["cy"]  = pv_a1_10x["y_hires"]
    pnn_a1_10x["cx"] = pnn_a1_10x["x_hires"]
    pnn_a1_10x["cy"] = pnn_a1_10x["y_hires"]
    render_crop(
        c1_10x_a1, c2_10x_a1, cx1_10x, cy1_10x,
        pv_a1_10x, pnn_a1_10x,
        CROP_HALF_10X, COLOC_PX_10X, PX_PER_UM_10X,
        "Control", A1_SLICE, False,
        OUT_DIR / f"crop_{A1_ANIMAL}_{A1_SLICE}_10x.png"
    )

# ==============================================================================
# ANIMAL 2 (TREATED)
# ==============================================================================

print(f"\n{'='*60}")
print(f"Animal 2: {A2_ANIMAL}  slice {A2_SLICE}  (40× crop: {A2_USE_40X})")
print(f"{'='*60}")

tiff_10x_a2 = TIFF_ROOT / A2_TREATMENT / A2_ANIMAL / f"{A2_ANIMAL}_{A2_SLICE}.tiff"
c1_10x_a2, c2_10x_a2 = load_tiff_channels(tiff_10x_a2)
H2_10x, W2_10x = c1_10x_a2.shape

print("  Loading 40× enwrapment CSV ...")
pv_a2_40x, pnn_a2_40x = load_cells_40x(A2_ANIMAL, A2_SLICE, A2_ENW_CSV)
pv_a2_40x  = pv_a2_40x.rename(columns={"x": "cx", "y": "cy"})
pnn_a2_40x = pnn_a2_40x.rename(columns={"x": "cx", "y": "cy"})

cells_a2 = load_cells(A2_ANIMAL, A2_SLICE)

if A2_USE_40X:
    print("  Loading 40× TIFF ...")
    c1_40x_a2, c2_40x_a2 = load_tiff_channels(A2_TIFF_40X)
    H2_40x, W2_40x = c1_40x_a2.shape
    print(f"  40× shape: {H2_40x}×{W2_40x}")

    # Select the crop window with the lowest enwrapment ratio among focused regions.
    # This surfaces areas where PNN degradation is clearly visible.
    pv_all_40x = pv_a2_40x.copy()
    margin = CROP_HALF_40X + 20
    pv_all_40x = pv_all_40x[
        (pv_all_40x["cx"] > margin) & (pv_all_40x["cx"] < W2_40x - margin) &
        (pv_all_40x["cy"] > margin) & (pv_all_40x["cy"] < H2_40x - margin)
    ]

    pv_enw_cands = pv_all_40x[pv_all_40x["enwrapped"]]
    print(f"  Scoring {len(pv_enw_cands)} enwrapped PV candidates ...")

    scores = []
    for _, row in pv_enw_cands.iterrows():
        ccx, ccy = int(row["cx"]), int(row["cy"])
        in_win = (
            (pv_all_40x["cx"] >= ccx - CROP_HALF_40X) &
            (pv_all_40x["cx"] <  ccx + CROP_HALF_40X) &
            (pv_all_40x["cy"] >= ccy - CROP_HALF_40X) &
            (pv_all_40x["cy"] <  ccy + CROP_HALF_40X)
        )
        n_total     = int(in_win.sum())
        n_enwrapped = int(pv_all_40x.loc[in_win, "enwrapped"].sum())
        fs = focus_score_40x(c2_40x_a2, ccx, ccy, CROP_HALF_40X, H2_40x, W2_40x)
        scores.append({"cx": ccx, "cy": ccy,
                        "n_total": n_total, "n_enwrapped": n_enwrapped,
                        "enw_ratio": n_enwrapped / max(n_total, 1),
                        "focus": fs})

    scores_df = pd.DataFrame(scores)
    valid     = scores_df[scores_df["n_total"] >= MIN_PV_IN_WINDOW]
    if len(valid) == 0:
        valid = scores_df

    focus_thr = np.percentile(valid["focus"], 100 - FOCUS_PERCENTILE_THRESHOLD)
    focused   = valid[valid["focus"] >= focus_thr]
    if len(focused) == 0:
        focused = valid

    best2   = focused.sort_values(["enw_ratio", "n_total"],
                                   ascending=[True, False]).iloc[A2_SKIP_TOP_N]
    cx2_40x = int(best2["cx"])
    cy2_40x = int(best2["cy"])
    cx2_10x = cx2_40x // WARP_SCALE
    cy2_10x = cy2_40x // WARP_SCALE
    print(f"  Selected 40× crop centre: ({cx2_40x}, {cy2_40x})  "
          f"enwrapment={best2['n_enwrapped']:.0f}/{best2['n_total']:.0f} "
          f"({best2['enw_ratio']*100:.0f}%)  focus={best2['focus']:.1f}")

    cx2_40x = max(CROP_HALF_40X, min(W2_40x - CROP_HALF_40X, cx2_40x))
    cy2_40x = max(CROP_HALF_40X, min(H2_40x - CROP_HALF_40X, cy2_40x))

render_fullslide(
    c1_10x_a2, c2_10x_a2, cx2_10x, cy2_10x, CROP_HALF_10X,
    PX_PER_UM_10X, "Treated", A2_SLICE,
    OUT_DIR / f"fullslide_{A2_ANIMAL}_{A2_SLICE}.png"
)

if A2_USE_40X:
    render_crop(
        c1_40x_a2, c2_40x_a2, cx2_40x, cy2_40x,
        pv_a2_40x, pnn_a2_40x,
        CROP_HALF_40X, COLOC_PX_40X, PX_PER_UM_40X,
        "Treated", A2_SLICE, True,
        OUT_DIR / f"crop_{A2_ANIMAL}_{A2_SLICE}_40x.png"
    )
else:
    pv_a2_10x  = cells_a2[cells_a2["cell_type"] == "PV"].copy()
    pnn_a2_10x = cells_a2[cells_a2["cell_type"] == "PNN"].copy()
    pv_a2_10x["cx"]  = pv_a2_10x["x_hires"]
    pv_a2_10x["cy"]  = pv_a2_10x["y_hires"]
    pnn_a2_10x["cx"] = pnn_a2_10x["x_hires"]
    pnn_a2_10x["cy"] = pnn_a2_10x["y_hires"]
    render_crop(
        c1_10x_a2, c2_10x_a2, cx2_10x, cy2_10x,
        pv_a2_10x, pnn_a2_10x,
        CROP_HALF_10X, COLOC_PX_10X, PX_PER_UM_10X,
        "Treated", A2_SLICE, False,
        OUT_DIR / f"crop_{A2_ANIMAL}_{A2_SLICE}_10x.png"
    )

# ==============================================================================
# SCHEMATIC
# ==============================================================================

render_schematic(OUT_DIR / "schematic_enwrapment.png")

print("\nAll done.")
