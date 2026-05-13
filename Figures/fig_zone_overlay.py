#!/usr/bin/env python3
"""
fig_zone_overlay.py
===================
Generates a five-panel figure illustrating injection zone assignment
for a single animal across one or more tissue sections.

Panel layout (assembled externally, e.g. in PowerPoint or Inkscape):
    A  WFA / PNN channel (C1)
    B  Interneuron marker channel (C2)
    C  Injection marker channel (C3)
    D  RGB composite  (C1 · C2 · C3)
    E  Composite + injection-zone overlay

One PNG is produced per panel per slice and saved into a per-slice
subfolder under OUT_DIR.  Set SLICES_TO_RENDER to render multiple
sections from the same animal in a single run.

Zone masks must be binary TIFFs in MASK_BASE named:
    MaskLoose_C3-<ANIMAL>_<slice_id>.tif   (penumbra outer boundary)
    MaskStrict_C3-<ANIMAL>_<slice_id>.tif  (core inner boundary)
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import tifffile
from pathlib import Path
from skimage.transform import resize as sk_resize

# ==============================================================================
# CONFIGURATION
# ==============================================================================

ANIMAL    = "ANIMAL_ID"        # Folder / filename prefix for the animal
TREATMENT = "GROUP_NAME"       # Parent folder (treatment group)

SLICES_TO_RENDER = [
    "s000",                    # Add slice IDs to render
]

# ── Paths ──────────────────────────────────────────────────────────────────────
TIFF_BASE = Path("/path/to/originals") / TREATMENT / ANIMAL
MASK_BASE = Path("/path/to/analysis_results")
OUT_DIR   = Path("/path/to/output/injection_overlay")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Channel indices (0-based within the multi-channel TIFF) ───────────────────
CH_WFA = 0   # C1 — WFA / PNN
CH_PV  = 1   # C2 — Interneuron marker (e.g. PV)
CH_MS  = 2   # C3 — Injection marker (e.g. mScarlet)

# ── Scale bar ──────────────────────────────────────────────────────────────────
PX_PER_UM   = 1 / 0.65   # Pixels per µm — adjust for your objective/camera
SCALEBAR_UM = 1000        # Physical scale bar length in µm
SCALEBAR_LW = 10          # Line width
SCALEBAR_FS = 32          # Font size

# ── Output ─────────────────────────────────────────────────────────────────────
DPI   = 300
MAX_W = 2000   # Downsample width for rendering; larger images are rescaled

# ── Colours ────────────────────────────────────────────────────────────────────
COL_WFA      = "#00FF00"             # WFA / PNN channel
COL_PV       = "#FF44FF"             # Interneuron marker channel
COL_MS       = "#00FFFF"             # Injection marker channel
COL_CORE     = (0.00, 1.00, 1.00)   # Core zone tint (RGB, normalised)
COL_PENUMBRA = (1.00, 1.00, 0.00)   # Penumbra zone tint (RGB, normalised)
ZONE_ALPHA   = 0.40                  # Zone overlay opacity

# ==============================================================================
# COLOURMAPS
# ==============================================================================

cmap_wfa = LinearSegmentedColormap.from_list("wfa", ["black", COL_WFA])
cmap_pv  = LinearSegmentedColormap.from_list("pv",  ["black", COL_PV])
cmap_ms  = LinearSegmentedColormap.from_list("ms",  ["black", COL_MS])

# ==============================================================================
# HELPERS
# ==============================================================================

def norm(img, lo_pct=0.5, hi_pct=99.5):
    """Percentile-stretch a 2-D array to [0, 1]."""
    lo = np.percentile(img, lo_pct)
    hi = np.percentile(img, hi_pct)
    if hi == lo:
        return np.zeros_like(img, dtype=float)
    return np.clip((img.astype(float) - lo) / (hi - lo), 0, 1)


# ==============================================================================
# PER-SLICE RENDER
# ==============================================================================

for slice_id in SLICES_TO_RENDER:

    tiff_path   = TIFF_BASE / f"{ANIMAL}_{slice_id}.tiff"
    loose_path  = MASK_BASE / f"MaskLoose_C3-{ANIMAL}_{slice_id}.tif"
    strict_path = MASK_BASE / f"MaskStrict_C3-{ANIMAL}_{slice_id}.tif"

    missing = [p for p in (tiff_path, loose_path, strict_path) if not p.exists()]
    if missing:
        print(f"[SKIP {slice_id}] Missing files:")
        for p in missing:
            print(f"  {p}")
        continue

    print(f"\nRendering {slice_id} ...")

    # ── Load TIFF ──────────────────────────────────────────────────────────────
    raw = tifffile.imread(str(tiff_path))
    print(f"  TIFF shape: {raw.shape}, dtype: {raw.dtype}")

    if raw.ndim == 3 and raw.shape[0] <= 5:
        c1 = raw[CH_WFA]; c2 = raw[CH_PV]; c3 = raw[CH_MS]
    elif raw.ndim == 3 and raw.shape[2] <= 5:
        c1 = raw[:, :, CH_WFA]; c2 = raw[:, :, CH_PV]; c3 = raw[:, :, CH_MS]
    else:
        print(f"  [SKIP] Unexpected TIFF shape: {raw.shape}")
        continue

    H, W = c1.shape

    # ── Load zone masks ────────────────────────────────────────────────────────
    loose  = tifffile.imread(str(loose_path)).astype(bool)
    strict = tifffile.imread(str(strict_path)).astype(bool)

    if loose.shape != c1.shape:
        print(f"  Resizing masks from {loose.shape} to {c1.shape}")
        loose  = sk_resize(loose,  c1.shape, order=0, preserve_range=True).astype(bool)
        strict = sk_resize(strict, c1.shape, order=0, preserve_range=True).astype(bool)

    core     = strict
    penumbra = loose & ~strict

    # ── Normalise channels ────────────────────────────────────────────────────
    c1n = norm(c1); c2n = norm(c2); c3n = norm(c3)

    # ── RGB composite ─────────────────────────────────────────────────────────
    # C1  → green  (G only)
    # C2  → magenta (R + B, no G)
    # C3  → shifts transduced pixels toward cyan via G + B
    rgb = np.zeros((H, W, 3), dtype=float)
    rgb[:, :, 0] = np.clip(c2n,       0, 1)   # R: C2
    rgb[:, :, 1] = np.clip(c1n + c3n, 0, 1)   # G: C1 + C3
    rgb[:, :, 2] = np.clip(c2n + c3n, 0, 1)   # B: C2 + C3

    # ── Zone overlay ──────────────────────────────────────────────────────────
    rgb_zones = rgb.copy() * 0.78   # Dim base image slightly under tint

    for mask, colour in [(core, COL_CORE), (penumbra, COL_PENUMBRA)]:
        for ch, val in enumerate(colour):
            rgb_zones[mask, ch] = np.clip(
                rgb_zones[mask, ch] * (1 - ZONE_ALPHA) + val * ZONE_ALPHA, 0, 1
            )
    rgb_zones[~loose] = rgb[~loose]   # Restore undimmed composite outside zones

    # ── Panel definitions ──────────────────────────────────────────────────────
    panels = [
        ("A", "C1_WFA_PNN",       c1n,       cmap_wfa),
        ("B", "C2_PV",            c2n,       cmap_pv),
        ("C", "C3_injection",     c3n,       cmap_ms),
        ("D", "Composite",        rgb,       None),
        ("E", "Composite_zones",  rgb_zones, None),
    ]

    slice_dir = OUT_DIR / slice_id
    slice_dir.mkdir(parents=True, exist_ok=True)

    for label, name, img_data, cmap in panels:

        # Downsample if the source exceeds MAX_W
        if W > MAX_W:
            scale  = MAX_W / W
            disp_h = int(H * scale)
            if img_data.ndim == 2:
                disp = sk_resize(img_data, (disp_h, MAX_W),
                                 order=1, preserve_range=True).astype(float)
            else:
                disp = sk_resize(img_data, (disp_h, MAX_W, img_data.shape[2]),
                                 order=1, preserve_range=True).astype(float)
            disp_w = MAX_W
        else:
            disp   = img_data
            disp_h, disp_w = H, W

        fig, ax = plt.subplots(figsize=(10, 10 * disp_h / disp_w), dpi=DPI)
        fig.patch.set_facecolor("black")
        ax.set_facecolor("black")
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        plt.tight_layout(pad=0)

        if cmap is not None:
            ax.imshow(disp, cmap=cmap, interpolation="lanczos", aspect="equal")
        else:
            ax.imshow(np.clip(disp, 0, 1), interpolation="lanczos", aspect="equal")

        ax.set_xlim(0, disp_w)
        ax.set_ylim(disp_h, 0)

        # ── Scale bar ─────────────────────────────────────────────────────────
        # Bar length is corrected for any downsampling applied above.
        bar_px = (PX_PER_UM * SCALEBAR_UM) * (disp_w / W)
        sb_x0  = disp_w * 0.05
        sb_x1  = sb_x0 + bar_px
        sb_y   = disp_h - disp_h * 0.05
        ax.plot([sb_x0, sb_x1], [sb_y, sb_y], color="white",
                lw=SCALEBAR_LW, solid_capstyle="butt", transform=ax.transData)
        ax.text((sb_x0 + sb_x1) / 2, sb_y - disp_h * 0.030,
                f"{SCALEBAR_UM} µm",
                color="white", ha="center", va="bottom",
                fontsize=SCALEBAR_FS, fontweight="bold",
                transform=ax.transData)

        # ── Panel label ───────────────────────────────────────────────────────
        ax.text(0.015, 0.97, label,
                transform=ax.transAxes, fontsize=18,
                fontweight="bold", color="white", va="top", ha="left")

        out_png = slice_dir / f"{label}_{name}_{slice_id}.png"
        fig.savefig(str(out_png), dpi=DPI, bbox_inches="tight",
                    facecolor=fig.get_facecolor())
        plt.close(fig)
        print(f"  Saved: {out_png.name}")

print("\nAll done.")
