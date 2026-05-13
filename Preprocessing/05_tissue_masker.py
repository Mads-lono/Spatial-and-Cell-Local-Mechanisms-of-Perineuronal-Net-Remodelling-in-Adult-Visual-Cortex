#!/usr/bin/env python3
"""
05_tissue_masker.py
===================
Converts ilastik "Simple Segmentation" PNG outputs to 1-bit binary tissue
masks. ilastik encodes the segmentation as pixel values: 1 = tissue (foreground),
2 = background. This script inverts that convention to a standard binary mask
(255 = tissue, 0 = background) and saves as 1-bit PNG.

The "_Simple Segmentation" suffix added by ilastik is stripped from output
filenames so that downstream scripts can match masks to slices by animal ID
and slice ID.

INPUT:  /path/to/ilastik_output/{animal}/thumbnails/*.png
OUTPUT: /path/to/DATA/{animal}/masks/*.png  (1-bit binary PNGs)

Usage:
    python 05_tissue_masker.py

Requirements:
    pip install pillow numpy
"""

import os
import numpy as np
from PIL import Image

# =============================================================================
# CONFIGURATION
# =============================================================================
# Folder containing ilastik Simple Segmentation PNG files
INPUT_FOLDER  = "/path/to/ilastik_output"

# Destination for 1-bit binary tissue masks
OUTPUT_FOLDER = "/path/to/DATA/masks"

# =============================================================================
# CONVERSION
# =============================================================================

def generate_1bitmask(input_folder: str, output_folder: str) -> None:
    """
    Convert all ilastik segmentation PNGs in input_folder to 1-bit binary masks.

    Pixel value 1 in the ilastik output = tissue (foreground).
    Pixel value 2 = background.
    Output is a 1-bit PNG where white (1) = tissue, black (0) = background.
    """
    os.makedirs(output_folder, exist_ok=True)
    processed = errors = 0

    for filename in sorted(os.listdir(input_folder)):
        if not filename.endswith(".png"):
            continue

        img_path = os.path.join(input_folder, filename)
        new_name = filename.replace("_Simple Segmentation", "")
        out_path = os.path.join(output_folder, new_name)

        try:
            img_array   = np.array(Image.open(img_path))
            binary_mask = np.where(img_array == 1, 1, 0).astype(np.uint8)
            binary_img  = Image.fromarray(binary_mask * 255).convert("1")
            binary_img.save(out_path)
            print(f"  ✓ {filename} → {new_name}")
            processed += 1
        except Exception as exc:
            print(f"  ✗ {filename}: {exc}")
            errors += 1

    print(f"\nDone. Processed: {processed}  Errors: {errors}")
    print(f"Output: {output_folder}")


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    generate_1bitmask(INPUT_FOLDER, OUTPUT_FOLDER)
