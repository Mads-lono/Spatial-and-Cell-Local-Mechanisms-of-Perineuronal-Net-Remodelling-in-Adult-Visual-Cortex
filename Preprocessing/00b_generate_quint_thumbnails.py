# =============================================================================
# 00b_generate_quint_thumbnails.py
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
#   Converts multi-channel TIFFs into contrast-enhanced RGB JPEG thumbnails
#   at 1500 px width for use in the QUINT atlas registration workflow
#   (QuickNII / VisuAlign). The 1500 px width matches the atlas overlay
#   resolution used throughout the downstream pipeline.
#
#   For each TIFF:
#     1. Split into individual channel images.
#     2. Apply histogram stretch (0.35% saturated) + equalisation per channel.
#     3. Merge channels back into a composite RGB image.
#     4. Rescale to 1500 × (proportional height) using bilinear interpolation.
#     5. Save as JPEG in {output_dir}/{animal}/.
#
#   A dry-run mode prints expected output paths without processing any files.
#
# Parameters (filled via Fiji dialog):
#   images_sequence_dir — folder containing the multi-channel TIFFs
#   animal              — animal identifier used to name the output subfolder
#   images_output_dir   — root output directory
#   dry_run             — if checked, prints paths without processing
#
# OUTPUT: /path/to/output/{animal}/{original_name}.jpg
#
# Note: These JPEG thumbnails are not used for quantitative analysis —
#       they are registration aids only. The full-resolution 16-bit TIFFs
#       remain the source for all downstream quantification.
# =============================================================================

#@ File(label="Directory of the image sequence (TIFFs)", style="directory") images_sequence_dir
#@ String(label="Animal identifier", required=true, value="ANIMAL_1") animal
#@ File(label="Output directory", style="directory") images_output_dir
#@ boolean(label="Dry-run?", dry_run=true) dry_run

import os
from ij import IJ
from ij.plugin import ContrastEnhancer
from ij.plugin import ChannelSplitter
from ij.plugin import RGBStackMerge
from ij.plugin import RGBStackConverter
from ij.plugin import Scaler

images_sequence_dir = str(images_sequence_dir)
images_output_dir   = str(images_output_dir)

# Collect all TIFF files in the input directory
images = sorted(os.listdir(images_sequence_dir))
fnames = [
    os.path.join(images_sequence_dir, img)
    for img in images
    if img.lower().endswith((".tif", ".tiff"))
]

if not fnames:
    print("No TIFF files found in: " + images_sequence_dir)
    raise Exception("No input files found.")

# Output subfolder
save_folder = os.path.join(images_output_dir, animal)

# Target thumbnail width in pixels — matches atlas overlay resolution
TARGET_WIDTH = 1500

# Contrast enhancement saturation fraction (percentage of pixels to clip)
CONTRAST_SATURATED = 0.35

# Dry run: print expected paths and exit
if dry_run:
    print("Dry run enabled — the following files would be saved:")
    for img in images:
        if not img.lower().endswith((".tif", ".tiff")):
            continue
        save_path = os.path.join(save_folder, img).replace(".tiff", ".jpg").replace(".tif", ".jpg")
        print("  " + os.path.join(images_sequence_dir, img) + "  -->  " + save_path)
    print("All done (dry run).")

else:
    # Instantiate plugins once outside the loop
    merge            = RGBStackMerge()
    splitter         = ChannelSplitter()
    rgb_converter    = RGBStackConverter()
    scaler           = Scaler()
    contrast_enhancer = ContrastEnhancer()

    if not os.path.exists(save_folder):
        os.makedirs(save_folder)

    for file_path in fnames:
        print("Opening: " + os.path.basename(file_path))
        imp = IJ.openImage(file_path)

        # Split channels
        channels = splitter.split(imp)

        # Contrast enhancement per channel
        for ch in channels:
            contrast_enhancer.stretchHistogram(ch, CONTRAST_SATURATED)
            contrast_enhancer.equalize(ch)

        # Merge back to composite, then convert to single RGB image
        imp2 = merge.mergeChannels(channels, False)
        rgb_converter.convertToRGB(imp2)

        # Rescale to target width maintaining aspect ratio
        original_width  = imp2.getWidth()
        original_height = imp2.getHeight()
        new_height = int(TARGET_WIDTH / float(original_width) * original_height)

        print("  Rescaling %d x %d -> %d x %d" % (
            original_width, original_height, TARGET_WIDTH, new_height))

        imp2 = scaler.resize(imp2, TARGET_WIDTH, new_height, 1, "bilinear")

        # Save as JPEG
        save_path = os.path.join(save_folder, os.path.basename(file_path))
        save_path = save_path.replace(".tiff", "").replace(".tif", "")
        IJ.saveAs(imp2, "jpg", save_path)
        print("  Saved: " + save_path + ".jpg")
        imp2.close()

    print("\nAll done!")
