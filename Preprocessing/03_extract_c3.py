#!/usr/bin/env python3
"""
03_extract_c3.py
================
Extracts Channel 3 (index 2) from all multi-channel TIFFs found recursively
under the input directory and saves each as a flat (non-nested) C3-prefixed
TIFF in the output folder. The output folder is intentionally flat to simplify
the Cellpose batch step that follows (04_generate_injection_masks.py).

The script handles channel-first (C, H, W) and channel-last (H, W, C) layouts
with 3 or more channels. Already-extracted files (prefixed "C3-") are skipped,
as are any files already present in the output folder.

INPUT:  /path/to/originals/{group}/{animal}/*.tiff  (multi-channel TIFFs)
OUTPUT: /path/to/C3_images/C3-{animal}_{slice}.tiff  (flat, single-channel)

Usage:
    python 03_extract_c3.py

Requirements:
    pip install tifffile pillow
"""

import os
import tifffile
import numpy as np
from PIL import Image

# Suppress PIL decompression bomb warnings for large microscopy images
Image.MAX_IMAGE_PIXELS = None

# =============================================================================
# CONFIGURATION
# =============================================================================
# Root folder containing all treatment/animal subfolders with original TIFFs
INPUT_DIR = "/path/to/originals"

# Flat output folder for extracted C3 images
OUTPUT_DIR = "/path/to/C3_images"

# Subdirectory names to skip (e.g. thumbnail folders, output folder itself)
SKIP_DIRS = ["thum"]

# =============================================================================
# EXTRACTION
# =============================================================================

def extract_channel_3(img: np.ndarray) -> np.ndarray | None:
    """
    Return the third channel (index 2) from a multi-channel image array,
    or None if the array does not contain at least 3 channels.

    Supported layouts:
        (H, W, C) with C >= 3
        (C, H, W) with C >= 3
    """
    if img.ndim == 3:
        if img.shape[-1] >= 3:          # (H, W, C)
            return img[..., 2]
        elif img.shape[0] >= 3:         # (C, H, W)
            return img[2, :, :]
    return None


def run() -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print(f"Scanning {INPUT_DIR} for multi-channel TIFFs ...")
    processed = skipped_exists = skipped_noc3 = errors = 0

    for root, dirs, files in os.walk(INPUT_DIR):
        # Skip unwanted subdirectories
        if any(skip in root for skip in SKIP_DIRS):
            continue
        if OUTPUT_DIR in root:
            continue

        for filename in sorted(files):
            if not filename.lower().endswith((".tif", ".tiff")):
                continue
            if filename.startswith("C3-"):
                continue  # already an extracted file

            full_path  = os.path.join(root, filename)
            new_name   = f"C3-{filename}"
            save_path  = os.path.join(OUTPUT_DIR, new_name)

            if os.path.exists(save_path):
                skipped_exists += 1
                continue

            try:
                img = tifffile.imread(full_path)
                ch3 = extract_channel_3(img)

                if ch3 is None:
                    print(f"  [SKIP] Cannot identify C3 in {filename} "
                          f"(shape: {img.shape})")
                    skipped_noc3 += 1
                    continue

                tifffile.imwrite(save_path, ch3, compression="zlib")
                processed += 1
                print(f"  [{processed}] Extracted: {new_name}")

            except Exception as exc:
                print(f"  [ERROR] {filename}: {exc}")
                errors += 1

    print("\n" + "-" * 50)
    print(f"Extraction complete.")
    print(f"  Processed:           {processed}")
    print(f"  Skipped (exists):    {skipped_exists}")
    print(f"  Skipped (no C3):     {skipped_noc3}")
    print(f"  Errors:              {errors}")
    print(f"  Output:              {OUTPUT_DIR}")


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    run()
