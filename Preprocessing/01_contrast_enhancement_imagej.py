# =============================================================================
# 01_contrast_enhancement_imagej.py
# =============================================================================
# *** IMPORTANT: This is a Jython script. It must be run inside Fiji/ImageJ. ***
# *** It will NOT work with a standard Python interpreter.                    ***
#
# How to run:
#   1. Open Fiji
#   2. Plugins > Macros > Edit...  (or Script Editor, language = Python)
#   3. Paste this file and click Run
#   4. Select the input folder (split 8-bit TIFFs from script 00)
#   5. Select the output folder (will be created if it does not exist)
#
# What it does:
#   Walks all subdirectories of the input folder, finds -C1.tif and -C2.tif
#   files (skips C3), applies:
#       1. Subtract Background (rolling ball radius = 50 px)
#       2. Enhance Contrast (saturated = 0.15%, normalize)
#   and saves the result to a mirrored output directory structure.
#   Files that already exist in the output are skipped.
#   RAM is freed after every file to prevent memory exhaustion on large batches.
#
# Method: Background subtraction + contrast enhancement (Method 3 Opt5).
#   Parameters were empirically optimised to maximise CPN detection yield.
#
# INPUT:  /path/to/split_8bit/{group}/{animal}/*-C1.tif, *-C2.tif
# OUTPUT: /path/to/split_8bit_enhanced/{group}/{animal}/*-C1.tif, *-C2.tif
# =============================================================================

from ij import IJ
from ij.io import DirectoryChooser
import os
import sys

# =============================================================================
# CONFIGURATION
# =============================================================================
# Rolling ball radius for background subtraction (pixels)
ROLLING_RADIUS = 50

# Saturation percentage for contrast enhancement
SATURATED_PCT = 0.15

# Channel suffixes to process (C3 is deliberately excluded)
PROCESS_CHANNELS = ("-C1.tif", "-C2.tif", "-C1.tiff", "-C2.tiff")

# Free RAM after every N files (lower = more frequent GC, slower processing)
GC_EVERY_N = 5

# =============================================================================
# HELPERS
# =============================================================================

def is_target_channel(filename):
    """Return True if the file is a C1 or C2 TIFF that should be processed."""
    fname_lower = filename.lower()
    return any(fname_lower.endswith(suffix.lower()) for suffix in PROCESS_CHANNELS)


def clear_memory():
    """Force ImageJ garbage collection."""
    IJ.freeMemory()
    IJ.run("Collect Garbage")


def get_subdirs_with_targets(parent_dir):
    """Recursively find all subdirectories containing C1/C2 TIF files."""
    result = []
    for root, dirs, files in os.walk(parent_dir):
        if any(is_target_channel(f) for f in files):
            result.append(root)
    return result

# =============================================================================
# MAIN
# =============================================================================

print("=" * 60)
print("Contrast Enhancement")
print("Background subtraction (rolling=" + str(ROLLING_RADIUS) + ")")
print("Enhance Contrast (saturated=" + str(SATURATED_PCT) + "%, normalize)")
print("Processing C1 and C2 channels only (C3 skipped)")
print("=" * 60)
print("")

# ── Select directories ────────────────────────────────────────────────────────
print("Select PARENT INPUT directory ...")
input_chooser = DirectoryChooser("Choose Parent Input Directory")
parent_input  = input_chooser.getDirectory()

if not parent_input:
    print("No input directory selected. Exiting.")
    sys.exit()

print("Select PARENT OUTPUT directory ...")
output_chooser = DirectoryChooser("Choose Parent Output Directory")
parent_output  = output_chooser.getDirectory()

if not parent_output:
    print("No output directory selected. Exiting.")
    sys.exit()

# Safety: refuse if input == output
if os.path.abspath(parent_input) == os.path.abspath(parent_output):
    print("ERROR: Input and output are the same directory. Aborting.")
    sys.exit()

print("")
print("Input:  " + parent_input)
print("Output: " + parent_output)
print("")

# ── Discover directories ──────────────────────────────────────────────────────
print("Scanning directory structure ...")
subdirs = get_subdirs_with_targets(parent_input)

if not subdirs:
    print("No subdirectories with C1/C2 images found. Exiting.")
    sys.exit()

print("Found " + str(len(subdirs)) + " subdirectory/ies with target files.")
print("")

# ── Process ───────────────────────────────────────────────────────────────────
total_processed = 0
total_skipped   = 0
total_errors    = 0
file_counter    = 0

for subdir in subdirs:
    # Compute mirrored output path
    rel_path  = subdir[len(parent_input):]
    if rel_path.startswith(os.sep):
        rel_path = rel_path[1:]
    out_subdir = os.path.join(parent_output, rel_path)

    if not os.path.exists(out_subdir):
        os.makedirs(out_subdir)

    files = sorted(os.listdir(subdir))
    target_files = [f for f in files if is_target_channel(f)]

    if not target_files:
        continue

    print("-- " + rel_path + " (" + str(len(target_files)) + " files) --")

    for fname in target_files:
        in_path  = os.path.join(subdir, fname)
        out_path = os.path.join(out_subdir, fname)

        if os.path.exists(out_path):
            print("  SKIP: " + fname)
            total_skipped += 1
            continue

        print("  Processing: " + fname)
        try:
            imp = IJ.openImage(in_path)
            if imp is None:
                print("  ERROR: Could not open " + fname)
                total_errors += 1
                continue

            IJ.run(imp, "Subtract Background...",
                   "rolling=" + str(ROLLING_RADIUS))
            IJ.run(imp, "Enhance Contrast",
                   "saturated=" + str(SATURATED_PCT) + " normalize")
            IJ.saveAs(imp, "Tiff", out_path)
            imp.close()

            total_processed += 1
            file_counter    += 1
            print("  Saved: " + fname)

        except Exception as e:
            print("  ERROR processing " + fname + ": " + str(e))
            total_errors += 1

        # Periodic RAM clear
        if file_counter % GC_EVERY_N == 0:
            clear_memory()

# Final GC pass
clear_memory()

print("")
print("=" * 60)
print("Done.")
print("  Processed: " + str(total_processed))
print("  Skipped:   " + str(total_skipped))
print("  Errors:    " + str(total_errors))
print("  Output:    " + parent_output)
print("=" * 60)
