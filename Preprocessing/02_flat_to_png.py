#!/usr/bin/env python3
"""
02_flat_to_png.py
=================
Converts VisuAlign .flat binary atlas overlay files to RGB PNGs, preserving
the input folder structure in the output directory. The .flat format encodes
each pixel as an integer region index. This script uses the Rainbow 2017
colour lookup JSON to map region indices to RGB colours, producing the atlas
overlay PNGs consumed by the Diffu_Spot quantification script (08).

INPUT:  /path/to/visuAlign_flat/{animal}/.../*.flat
        /path/to/Rainbow_2017.json
OUTPUT: /path/to/flat_png_output/  (mirrored folder structure)

Usage:
    python 02_flat_to_png.py

Requirements:
    pip install numpy pillow
"""

import struct
import json
import numpy as np
from pathlib import Path
from PIL import Image

# =============================================================================
# CONFIGURATION
# =============================================================================
# Folder containing per-animal subfolders with .flat files from VisuAlign
INPUT_DIR = Path("/path/to/visuAlign_flat")

# Destination folder — mirrored structure with PNGs
OUTPUT_DIR = Path("/path/to/flat_png_output")

# Allen Brain Atlas colour lookup table (Rainbow 2017 JSON)
JSON_LOOKUP_PATH = Path("/path/to/Rainbow_2017.json")

# =============================================================================
# CORE DECODING
# =============================================================================

def load_palette(json_path: Path) -> np.ndarray:
    """
    Build an RGB palette array from the Rainbow 2017 JSON.

    The .flat file stores region indices as 8-bit (unit=1) or 16-bit (unit=2)
    integers. Index 0 is background (black). Every other index maps to an RGB
    colour in the JSON.

    Returns an (N, 3) uint8 array where palette[index] = [R, G, B].
    """
    with open(json_path, "r") as f:
        regions = json.load(f)

    # Find the maximum region index to size the palette
    max_index = max((int(r.get("index", 0)) for r in regions), default=0)
    palette = np.zeros((max_index + 1, 3), dtype=np.uint8)

    for region in regions:
        idx = int(region.get("index", 0))
        if idx > 0 and region.get("name") != "empty":
            palette[idx] = [
                int(region.get("red",   0)),
                int(region.get("green", 0)),
                int(region.get("blue",  0)),
            ]

    return palette


def decode_flat_to_rgb(flat_path: Path, palette: np.ndarray) -> np.ndarray:
    """
    Decode a single .flat file into an RGB NumPy array.

    .flat header format (big-endian):
        unit   : 1 byte  — 1 = uint8 indices, 2 = uint16 indices
        width  : 4 bytes — image width in pixels
        height : 4 bytes — image height in pixels
    Followed by width × height index values of size `unit` bytes each.

    Returns an (height, width, 3) uint8 RGB array.
    """
    with open(flat_path, "rb") as f:
        unit, width, height = struct.unpack(">BII", f.read(9))
        dtype_char  = {1: "B", 2: "H"}[unit]
        n_pixels    = width * height
        raw_bytes   = f.read(unit * n_pixels)
        region_ids  = struct.unpack(f">{n_pixels}{dtype_char}", raw_bytes)

    id_array = np.array(region_ids, dtype=np.uint16).reshape((height, width))

    # Clamp out-of-range IDs to avoid index errors
    id_array = np.clip(id_array, 0, len(palette) - 1)

    return palette[id_array]  # (H, W, 3)

# =============================================================================
# BATCH CONVERSION
# =============================================================================

def convert_all(
    input_dir: Path,
    output_dir: Path,
    json_path: Path,
) -> None:
    """
    Recursively find all .flat files under input_dir, decode each to RGB,
    and save as PNG in a mirrored structure under output_dir.
    """
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    if not json_path.exists():
        raise FileNotFoundError(f"JSON lookup not found: {json_path}")

    print(f"Loading palette from {json_path.name} ...")
    palette = load_palette(json_path)
    print(f"  Palette built — {len(palette)} region entries.\n")

    flat_files = sorted(input_dir.rglob("*.flat"))
    if not flat_files:
        print(f"No .flat files found under {input_dir}. Nothing to do.")
        return

    print(f"Found {len(flat_files)} .flat file(s).\n")
    print("=" * 60)

    processed = skipped = errors = 0

    for flat_path in flat_files:
        rel      = flat_path.relative_to(input_dir)
        out_path = (output_dir / rel).with_suffix(".png")

        if out_path.exists():
            skipped += 1
            continue

        out_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            rgb = decode_flat_to_rgb(flat_path, palette)
            Image.fromarray(rgb, mode="RGB").save(out_path)
            print(f"  ✓ {rel}")
            processed += 1
        except Exception as exc:
            print(f"  ✗ {rel} — {exc}")
            errors += 1

    print("\n" + "=" * 60)
    print(f"Done.  Converted: {processed}  Skipped: {skipped}  Errors: {errors}")
    print(f"Output: {output_dir.resolve()}")

# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    convert_all(INPUT_DIR, OUTPUT_DIR, JSON_LOOKUP_PATH)
