#!/usr/bin/env python3
"""
00_split_and_convert_8bit.py
============================
Reads multi-channel TIFFs from a backup location (READ-ONLY), splits each
file into separate single-channel TIFFs (-C1, -C2, -C3), and converts the
output files to 8-bit in place. The original source files are never modified.
Folder structure is mirrored exactly from input to output.

INPUT:  /path/to/originals/{group}/{animal}/*.tiff  (multi-channel, 16-bit)
OUTPUT: /path/to/split_8bit/{group}/{animal}/{stem}-C1.tif  etc.

Usage:
    python 00_split_and_convert_8bit.py

Requirements:
    pip install numpy tifffile
"""

import numpy as np
from pathlib import Path
from typing import List, Union
from tifffile import imread, imwrite

# =============================================================================
# CONFIGURATION
# =============================================================================
# Source directory — original multi-channel TIFFs (READ-ONLY, never modified)
INPUT_DIR  = Path("/path/to/originals")

# Output directory — split and 8-bit converted files written here
OUTPUT_DIR = Path("/path/to/split_8bit")

# Glob pattern to match the original files
PATTERN = "*.tiff"

# =============================================================================
# PATH VALIDATION
# =============================================================================

def validate_paths(input_dir: Path, output_dir: Path) -> None:
    """Refuse to run if output overlaps with the input backup location."""
    inp = input_dir.resolve()
    out = output_dir.resolve()

    if inp == out:
        raise ValueError(
            f"SAFETY ERROR: input and output are the same directory.\n"
            f"  Input:  {inp}\n  Output: {out}\n"
            "This would modify the backup. Aborting."
        )
    try:
        out.relative_to(inp)
        raise ValueError(
            f"SAFETY ERROR: output is inside the input directory.\n"
            f"  Input:  {inp}\n  Output: {out}\n"
            "This could modify the backup. Aborting."
        )
    except ValueError as e:
        if "SAFETY" in str(e):
            raise
        # Good — output is not nested inside input

    if not input_dir.exists():
        raise ValueError(f"Input directory does not exist: {input_dir}")

    print(f"Path validation passed.")
    print(f"  Reading from: {inp}")
    print(f"  Writing to:   {out}")

# =============================================================================
# CHANNEL SPLITTING
# =============================================================================

def split_tiff_channels(input_path: Path, output_dir: Path) -> List[Path]:
    """
    Split a multi-channel TIFF into separate single-channel files.

    Supports channel-first (C, H, W) and channel-last (H, W, C) layouts
    with 2 or 3 channels. Output files are named {stem}-C1.tif, -C2.tif, -C3.tif.

    Returns a list of paths to the saved channel files.
    """
    img = imread(input_path)

    if img.ndim < 3:
        raise ValueError(
            f"Expected at least 3 dimensions, got {img.ndim}: {input_path.name}"
        )

    # Determine channel axis
    if img.shape[0] in (2, 3, 4):
        channel_axis, n_channels = 0, img.shape[0]
    elif img.shape[-1] in (2, 3, 4):
        channel_axis, n_channels = -1, img.shape[-1]
    else:
        raise ValueError(
            f"Cannot find 2–4 channel axis in shape {img.shape}: {input_path.name}"
        )

    stem = input_path.stem
    output_paths = []

    for i in range(n_channels):
        ch = img[i, ...] if channel_axis == 0 else img[..., i]
        out_path = output_dir / f"{stem}-C{i + 1}.tif"
        imwrite(out_path, ch)
        output_paths.append(out_path)

    return output_paths

# =============================================================================
# 8-BIT CONVERSION
# =============================================================================

def convert_to_8bit(file_path: Path) -> bool:
    """
    Convert a 16-bit TIFF to 8-bit in place by scaling to the image maximum.

    Returns True if conversion was performed, False if the file was already 8-bit.
    """
    img = imread(file_path)

    if img.dtype == np.uint8:
        return False

    if img.dtype != np.uint16:
        raise ValueError(f"Unexpected dtype {img.dtype}: {file_path.name}")

    max_val = img.max()
    img_8bit = (img / max_val * 255).astype(np.uint8) if max_val > 0 else img.astype(np.uint8)
    imwrite(file_path, img_8bit)
    return True

# =============================================================================
# SINGLE-IMAGE PROCESSING
# =============================================================================

def process_single_image(input_path: Path, output_dir: Path) -> None:
    """Split channels then convert each output file to 8-bit."""
    split_paths = split_tiff_channels(input_path, output_dir)
    for sp in split_paths:
        converted = convert_to_8bit(sp)
        status = "16→8bit" if converted else "already 8-bit"
        print(f"      ✓ {sp.name} ({status})")

# =============================================================================
# BATCH PROCESSING
# =============================================================================

def find_all_tiff_files(input_dir: Path, pattern: str) -> List[Path]:
    """Recursively find all files matching pattern under input_dir."""
    return sorted(input_dir.rglob(pattern))


def preview_structure(input_dir: Path, pattern: str) -> None:
    """Print the folder structure and file counts before processing."""
    from collections import defaultdict
    files = find_all_tiff_files(input_dir, pattern)
    if not files:
        print(f"No files matching '{pattern}' found in {input_dir}")
        return

    folders: dict = defaultdict(list)
    for f in files:
        folders[f.parent.relative_to(input_dir)].append(f.name)

    print(f"Found {len(files)} file(s) in {len(folders)} folder(s):\n")
    for folder in sorted(folders):
        flist = folders[folder]
        print(f"  {folder}/  ({len(flist)} files)")
        for name in flist[:3]:
            print(f"    - {name}")
        if len(flist) > 3:
            print(f"    - ... and {len(flist) - 3} more")


def batch_process_recursive(
    input_dir: Union[str, Path],
    output_dir: Union[str, Path],
    pattern: str = "*.tiff",
) -> None:
    """
    Process all matching TIFFs recursively, preserving folder structure.

    For each file: split channels → convert each channel to 8-bit.
    Only the output directory is written to; the input is never touched.
    """
    input_dir  = Path(input_dir)
    output_dir = Path(output_dir)

    validate_paths(input_dir, output_dir)

    tiff_files = find_all_tiff_files(input_dir, pattern)
    if not tiff_files:
        print(f"No files matching '{pattern}' found. Nothing to do.")
        return

    subdirs = {f.parent.relative_to(input_dir) for f in tiff_files}
    print(f"\nFound {len(tiff_files)} file(s) across {len(subdirs)} folder(s).\n")
    print("=" * 60)

    current_subdir = None
    for i, tiff_file in enumerate(tiff_files, 1):
        rel = tiff_file.parent.relative_to(input_dir)
        if rel != current_subdir:
            current_subdir = rel
            print(f"\n  {rel}/")

        file_out_dir = output_dir / rel
        file_out_dir.mkdir(parents=True, exist_ok=True)

        print(f"  [{i}/{len(tiff_files)}] {tiff_file.name}")
        try:
            process_single_image(tiff_file, file_out_dir)
        except Exception as exc:
            print(f"      ✗ Error: {exc}")

    print("\n" + "=" * 60)
    print("Batch processing complete.")
    print(f"Output: {output_dir.resolve()}")


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    # Optional: preview structure before running
    preview_structure(INPUT_DIR, PATTERN)
    print()
    batch_process_recursive(INPUT_DIR, OUTPUT_DIR, PATTERN)
