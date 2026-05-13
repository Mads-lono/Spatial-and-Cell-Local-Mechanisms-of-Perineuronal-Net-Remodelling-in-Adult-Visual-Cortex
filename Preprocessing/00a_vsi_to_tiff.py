# =============================================================================
# 00a_vsi_to_tiff.py
# =============================================================================
# *** IMPORTANT: This is a Jython script. It must be run inside Fiji/ImageJ. ***
# *** It will NOT work with a standard Python interpreter.                    ***
#
# How to run:
#   1. Open Fiji
#   2. Plugins > Macros > Edit... (or Script Editor, language = Python)
#   3. Paste this file and click Run
#   4. Fill in the parameter dialog that appears
#
# What it does:
#   Converts raw .vsi microscope files (Olympus CellSens format) to TIFF,
#   applying a systematic sXXX slice index with a configurable starting point
#   and a 2-step increment. This produces the canonical renamed TIFF files
#   that all downstream scripts (00b, 00, 01 ...) read from.
#
#   Files are sorted by the trailing number in their filename before indexing,
#   ensuring consistent ordering regardless of OS directory listing order.
#   A dry-run mode prints expected output paths without opening any files.
#
# Parameters (filled via Fiji dialog):
#   images_sequence_dir  — folder containing .vsi files for one animal/series
#   group                — animal identifier and series label, e.g. "ANIMAL_1_A1"
#                          the last "_XX" suffix is stripped for the output folder name
#   images_output_dir    — root output directory
#   dry_run              — if checked, prints paths without processing
#   index_starting_point — sXXX index of the first slice (default 1)
#
# OUTPUT: /path/to/output/{animal}/{animal}_A1_sXXX.tiff
#         Index increments by 2 per file (every-other-section convention).
#
# Requires: Bio-Formats plugin (bundled with Fiji)
# =============================================================================

#@ File(label="Directory of the image sequence (.vsi files)", style="directory") images_sequence_dir
#@ String(label="Animal and series (e.g. ANIMAL_1_A1)", required=true, value="ANIMAL_1_A1") group
#@ File(label="Output directory", style="directory") images_output_dir
#@ boolean(label="Dry-run?", dry_run=true) dry_run
#@ Integer(label="Index starting point (sXXX)?", value=1) index_starting_point

import os
import re
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from ij import IJ

images_sequence_dir = str(images_sequence_dir)
images_output_dir   = str(images_output_dir)

# Strip the last suffix from group to get the output folder name
# e.g. "ANIMAL_1_A1" -> output folder "ANIMAL_1"
output_folder_name = group.rsplit("_", 1)[0]
pathout = os.path.join(images_output_dir, output_folder_name)

if not os.path.exists(pathout):
    os.makedirs(pathout)

# Collect .vsi files
fnames = [
    os.path.join(images_sequence_dir, f)
    for f in os.listdir(images_sequence_dir)
    if f.endswith(".vsi")
]

if len(fnames) < 1:
    raise Exception("No .vsi files found in: %s" % images_sequence_dir)


def extract_last_number(fname):
    """Extract the trailing integer before .vsi for sorting."""
    match = re.search(r"(\d+)(?=\.vsi$)", fname)
    return int(match.group(0)) if match else float("inf")


fnames = sorted(fnames, key=extract_last_number)
print("Found %d .vsi file(s):" % len(fnames))
for f in fnames:
    print("  " + f)

# Pre-compute all output paths
save_paths = []
index = index_starting_point
for _ in fnames:
    savepath = os.path.join(pathout, "%s_s%03d.tiff" % (group, index))
    save_paths.append(savepath)
    index += 2   # every-other-section convention

# Dry run: print expected paths and exit
if dry_run:
    print("\nDry run enabled — the following files would be saved:")
    for sp in save_paths:
        print("  " + sp)
    print("All done (dry run).")

else:
    options = ImporterOptions()
    options.setAutoscale(True)

    counter = 0
    index   = index_starting_point

    for fname in fnames:
        try:
            options.setId(fname)
        except Exception as e:
            print("Error opening: %s\n  %s" % (fname, str(e)))
            continue

        imps = BF.openImagePlus(options)

        for imp in imps:
            savepath = save_paths[counter]
            counter += 1

            imp.setTitle(group + str(counter))

            if os.path.exists(savepath):
                print("SKIP (already exists): " + savepath)
                imp.close()
                continue

            IJ.saveAsTiff(imp, savepath)
            print("Saved: " + savepath)
            imp.close()

    print("\nAll done! %d file(s) processed." % counter)
