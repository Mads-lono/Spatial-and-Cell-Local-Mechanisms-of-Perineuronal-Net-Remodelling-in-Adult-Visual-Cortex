#!/usr/bin/env python3
"""
04_generate_injection_masks.py
==============================
Generates injection zone masks for each C3 image using Cellpose (CPSAM model).

For each C3 TIFF produced by 03_extract_c3.py:
  1. Two threshold-based binary masks are created from the C3 fluorescence:
       MaskLoose  — permissive boundary (σ=6, relative threshold 15%)
       MaskStrict — tight injection core (σ=10, relative threshold 60%)
     These define the three downstream injection zones:
       Core     = MaskStrict (intersection with tissue mask)
       Penumbra = MaskLoose & ~MaskStrict
       Outside  = ~MaskLoose (within tissue)

  2. The image is cropped to the loose mask bounding box (+ padding) and
     passed to Cellpose (CPSAM pretrained model, GPU) for cell segmentation.

  3. Cells whose centroid falls within MaskLoose are retained. A watershed
     correction step splits merged detections using distance-transform peaks.

  4. Per-file outputs saved to OUTPUT_FOLDER:
       MaskLoose_{base}.tif    — binary mask (255 = inside, 0 = outside)
       MaskStrict_{base}.tif   — binary mask
       Counts_{base}.csv       — cell centroids and basic measurements
       SegLabels_{base}.tif    — labelled segmentation (uint16)

Files with all four outputs already present are skipped.

INPUT:  /path/to/C3_images/C3-*.tif*       (from 03_extract_c3.py)
OUTPUT: /path/to/analysis_results/          (masks, counts, segmentation)

Usage:
    conda activate cpgpu   # or whichever environment has Cellpose
    python 04_generate_injection_masks.py

Requirements:
    pip install cellpose tifffile scikit-image scipy pandas pillow numpy
"""

import os
import glob
import re
import numpy as np
import pandas as pd
import tifffile
from skimage.transform import rescale, resize
from skimage.filters import gaussian
from skimage.measure import label, regionprops
from skimage.morphology import binary_erosion, disk
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from scipy import ndimage as ndi
from cellpose import models
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

# =============================================================================
# CONFIGURATION
# =============================================================================
# Flat folder of C3 images produced by 03_extract_c3.py
INPUT_FOLDER  = "/path/to/C3_images"

# Destination for masks, counts CSV, and segmentation label images
OUTPUT_FOLDER = "/path/to/analysis_results"

# ── Mask settings ─────────────────────────────────────────────────────────────
# Loose mask — defines the Penumbra outer boundary and crops the image for
# Cellpose inference. Lower threshold = larger region.
LOOSE_SIGMA      = 6
LOOSE_REL_THRESH = 0.15   # 15% above background

# Strict mask — defines the Core zone. Higher threshold = tighter region.
STRICT_SIGMA      = 10
STRICT_REL_THRESH = 0.60  # 60% above background

# ── Cellpose settings ─────────────────────────────────────────────────────────
# Pretrained model to use. "cpsam" is the default CPSAM checkpoint.
CELLPOSE_MODEL   = "cpsam"
CELLPOSE_GPU     = True
FLOW_THRESHOLD   = 0.4
CELLPROB_THRESHOLD = 0.0

# ── Watershed correction ──────────────────────────────────────────────────────
ENABLE_WATERSHED  = True
WATERSHED_MIN_DIST = 5    # minimum distance between watershed peaks (pixels)

# ── Cropping padding (pixels added around the loose mask bounding box) ────────
PADDING = 100

# =============================================================================
# FILENAME PARSING
# =============================================================================

def parse_filename(filename: str):
    """
    Parse treatment, animal number, and slice ID from a C3 filename.

    Expected format: C3-{treatment}_{animal}_{slice}.tiff
    Example:         C3-ADAMTS15_3_s007.tiff
    Returns (treatment, animal, slice_id) or ("Unknown", "Unknown", "Unknown").
    """
    pattern = r"^C3-(.+)_([^_]+)_([sS]\d+)\.tiff?"
    match   = re.search(pattern, filename)
    if match:
        return match.group(1), match.group(2), match.group(3)
    return "Unknown", "Unknown", "Unknown"

# =============================================================================
# MASK GENERATION
# =============================================================================

def create_injection_mask(
    img: np.ndarray,
    original_shape: tuple,
    sigma: float = 5,
    rel_threshold: float = 0.2,
    erosion: int = 0,
) -> np.ndarray:
    """
    Create a binary injection zone mask from a C3 fluorescence image.

    The mask is computed at a downscaled resolution (2% of original) for speed,
    then upsampled back to full resolution with nearest-neighbour interpolation.

    Process:
      1. Downsample to 2% of original size.
      2. Gaussian blur (sigma).
      3. Compute background estimate as the 20th percentile of pixel values.
      4. Threshold = background + rel_threshold * (max - background).
         Triangle threshold is also computed; the higher of the two is used.
      5. Keep the largest connected component only.
      6. Optional binary erosion.
      7. Upsample to original_shape.

    Returns a binary uint8 mask at original_shape resolution.
    """
    scale = 0.02
    img_small   = rescale(img, scale, anti_aliasing=True, preserve_range=True)
    img_blurred = gaussian(img_small, sigma=sigma, preserve_range=True)

    background_val = np.percentile(img_blurred, 20)
    max_val        = np.max(img_blurred)

    if (max_val - background_val) < 10:
        return np.zeros(original_shape, dtype=np.uint8)

    thresh_rel = background_val + rel_threshold * (max_val - background_val)

    # Triangle threshold as a floor (avoids thresholding noise)
    from skimage.filters import threshold_triangle
    try:
        thresh_tri = threshold_triangle(img_blurred)
    except Exception:
        thresh_tri = 0

    final_thresh = max(thresh_tri, thresh_rel)
    binary_mask  = img_blurred > final_thresh

    # Keep largest connected component
    labeled  = label(binary_mask)
    regions  = regionprops(labeled)
    if not regions:
        return np.zeros(original_shape, dtype=np.uint8)

    largest   = max(regions, key=lambda r: r.area)
    mask_small = labeled == largest.label

    if erosion > 0:
        mask_small = binary_erosion(mask_small, footprint=disk(erosion))

    mask_full = resize(
        mask_small, original_shape,
        order=0, preserve_range=True, anti_aliasing=False
    )
    return (mask_full > 0).astype(np.uint8)

# =============================================================================
# WATERSHED CORRECTION
# =============================================================================

def apply_watershed_correction(labels: np.ndarray, min_dist: int = 3) -> np.ndarray:
    """
    Split merged Cellpose detections using distance-transform watershed.

    Each cell in the label image is expected to be a separate connected
    component. Distance transform peaks become watershed seeds.
    A rescue step ensures small cells that would otherwise lose their peak
    are still represented.
    """
    if not np.any(labels):
        return labels

    binary_mask = labels > 0
    distance    = ndi.distance_transform_edt(binary_mask)

    peaks_mask = np.zeros(distance.shape, dtype=bool)
    coords     = peak_local_max(distance, min_distance=min_dist, labels=binary_mask)
    if len(coords) > 0:
        peaks_mask[tuple(coords.T)] = True

    # Rescue any cells that have no peak assigned
    present_ids = np.unique(labels[peaks_mask])
    all_ids     = np.unique(labels)
    all_ids     = all_ids[all_ids != 0]
    missing_ids = np.setdiff1d(all_ids, present_ids)

    if len(missing_ids) > 0:
        objs = ndi.find_objects(labels)
        for mid in missing_ids:
            sl = objs[mid - 1]
            if sl is None:
                continue
            local_dist  = distance[sl]
            local_label = labels[sl]
            cell_dist   = np.where(local_label == mid, local_dist, 0)
            if np.max(cell_dist) > 0:
                max_y, max_x = np.unravel_index(np.argmax(cell_dist), cell_dist.shape)
                peaks_mask[sl[0].start + max_y, sl[1].start + max_x] = True

    markers, _ = ndi.label(peaks_mask)
    new_labels  = watershed(-distance, markers, mask=binary_mask)
    return new_labels.astype(labels.dtype)

# =============================================================================
# MAIN PROCESSING LOOP
# =============================================================================

def main() -> None:
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    print("Loading Cellpose model ...")
    model = models.CellposeModel(
        gpu=CELLPOSE_GPU, pretrained_model=CELLPOSE_MODEL
    )

    files = sorted(glob.glob(os.path.join(INPUT_FOLDER, "C3-*.tif*")))
    print(f"Found {len(files)} C3 image(s) to process.\n")

    for file_idx, filepath in enumerate(files):
        filename  = os.path.basename(filepath)
        treatment, animal, slice_id = parse_filename(filename)
        base_name = os.path.splitext(filename)[0]

        print(f"\n--- [{file_idx + 1}/{len(files)}] {filename} ---")

        path_loose  = os.path.join(OUTPUT_FOLDER, f"MaskLoose_{base_name}.tif")
        path_strict = os.path.join(OUTPUT_FOLDER, f"MaskStrict_{base_name}.tif")
        path_csv    = os.path.join(OUTPUT_FOLDER, f"Counts_{base_name}.csv")
        path_seg    = os.path.join(OUTPUT_FOLDER, f"SegLabels_{base_name}.tif")

        if all(os.path.exists(p) for p in [path_loose, path_strict, path_csv, path_seg]):
            print("  All outputs exist — skipping.")
            continue

        # ── Load image ────────────────────────────────────────────────────────
        try:
            img = tifffile.imread(filepath)
        except Exception as exc:
            print(f"  ERROR loading: {exc}")
            continue

        if img.ndim == 3:
            img = img[..., 0] if img.shape[-1] < 10 else img[0, :, :]

        # ── Generate masks ────────────────────────────────────────────────────
        print("  Generating strict mask ...")
        mask_strict = create_injection_mask(
            img, img.shape, sigma=STRICT_SIGMA, rel_threshold=STRICT_REL_THRESH
        )
        tifffile.imwrite(path_strict, mask_strict * 255, compression="zlib")

        print("  Generating loose mask ...")
        mask_loose = create_injection_mask(
            img, img.shape, sigma=LOOSE_SIGMA, rel_threshold=LOOSE_REL_THRESH
        )
        tifffile.imwrite(path_loose, mask_loose * 255, compression="zlib")

        # Guard: empty loose mask means no injection detected
        rows = np.any(mask_loose, axis=1)
        if not np.any(rows):
            print("  WARNING: empty loose mask — saving empty outputs.")
            pd.DataFrame(columns=["Filename", "Treatment", "Animal", "Slice"]).to_csv(
                path_csv, index=False
            )
            tifffile.imwrite(
                path_seg, np.zeros(img.shape, dtype=np.uint16), compression="zlib"
            )
            continue

        # ── Crop to loose mask bounding box ───────────────────────────────────
        cols = np.any(mask_loose, axis=0)
        y_min, y_max = np.where(rows)[0][[0, -1]]
        x_min, x_max = np.where(cols)[0][[0, -1]]
        y_min = max(0, y_min - PADDING)
        x_min = max(0, x_min - PADDING)
        y_max = min(mask_loose.shape[0], y_max + PADDING)
        x_max = min(mask_loose.shape[1], x_max + PADDING)

        img_crop        = img[y_min:y_max, x_min:x_max]
        mask_loose_crop = mask_loose[y_min:y_max, x_min:x_max]

        # ── Cellpose inference ────────────────────────────────────────────────
        print("  Running Cellpose (CPSAM) ...")
        masks_crop, _, _ = model.eval(
            img_crop,
            diameter=None,
            channels=[0, 0],
            normalize=True,
            flow_threshold=FLOW_THRESHOLD,
            cellprob_threshold=CELLPROB_THRESHOLD,
        )

        # ── Watershed correction ──────────────────────────────────────────────
        if ENABLE_WATERSHED:
            print("  Applying watershed correction ...")
            masks_crop = apply_watershed_correction(masks_crop, min_dist=WATERSHED_MIN_DIST)

        # ── Filter: keep cells whose centroid is inside loose mask ────────────
        print("  Filtering and saving ...")
        full_label_image = np.zeros(img.shape, dtype=np.uint16)
        valid_cells      = []

        for region in regionprops(masks_crop, intensity_image=img_crop):
            ly, lx = int(region.centroid[0]), int(region.centroid[1])

            if mask_loose_crop[ly, lx] == 0:
                continue

            gy = ly + y_min
            gx = lx + x_min

            valid_cells.append({
                "Filename":        filename,
                "Treatment":       treatment,
                "Animal":          animal,
                "Slice":           slice_id,
                "Cell_ID":         region.label,
                "Global_Y":        gy,
                "Global_X":        gx,
                "Area_px":         region.area,
                "Mean_Intensity":  region.intensity_mean,
            })

            # Reconstruct label image at full resolution
            gc_y = region.coords[:, 0] + y_min
            gc_x = region.coords[:, 1] + x_min
            valid = (gc_y < img.shape[0]) & (gc_x < img.shape[1])
            full_label_image[gc_y[valid], gc_x[valid]] = region.label

        if valid_cells:
            pd.DataFrame(valid_cells).to_csv(path_csv, index=False)
            print(f"  Saved {len(valid_cells)} cells.")
        else:
            print("  No cells found within loose injection area.")
            pd.DataFrame(
                columns=["Filename", "Treatment", "Animal", "Slice"]
            ).to_csv(path_csv, index=False)

        tifffile.imwrite(path_seg, full_label_image, compression="zlib")

    print("\n" + "=" * 60)
    print("Done.")


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    main()
