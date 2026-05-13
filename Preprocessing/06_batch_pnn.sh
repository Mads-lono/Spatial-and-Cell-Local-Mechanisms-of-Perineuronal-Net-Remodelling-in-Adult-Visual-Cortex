#!/bin/bash
# =============================================================================
# 06_batch_pnn.sh
# =============================================================================
# Runs CPN cell detection on all split 8-bit enhanced TIFFs for one channel.
#
# Set MODE to "PNN" (channel C1) or "PV" (channel C2). The mode controls
# which detection and scoring models are used and which channel files are
# targeted. PNN and PV are mutually exclusive — run the script twice, once
# per mode, pointing INPUT_BASE and OUTPUT_BASE at the appropriate directories.
#
# For each animal subfolder the script:
#   1. Locates all -C{channel}.tif files
#   2. Skips the folder if a non-empty localizations CSV already exists
#   3. Runs predict.py (Faster R-CNN detection + rank learning rescoring)
#   4. On success, runs draw_predictions.py to render visual output
#
# Usage:
#   conda activate cpngpu
#   cd /path/to/CPN
#   bash 06_batch_pnn.sh
#
# To run detached (overnight):
#   nohup bash 06_batch_pnn.sh > pnn_run.log 2>&1 &
# =============================================================================

# =============================================================================
# CONFIGURATION
# =============================================================================
# Set to "PNN" to detect PNN-enwrapped cells (C1 channel)
# Set to "PV"  to detect PV interneurons         (C2 channel)
MODE="PNN"

INPUT_BASE="/path/to/split_8bit_enhanced"
OUTPUT_BASE="/path/to/CPN_output"
CPN_DIR="/path/to/CPN"

# Detection threshold (0.0 = keep all detections; filter to >= 0.5 in post-processing)
DETECTION_THRESHOLD="0.0"

# CUDA device index
CUDA_DEVICE="cuda:0"

# =============================================================================
# MODE SETUP  (do not edit below this line)
# =============================================================================
if [ "$MODE" == "PNN" ]; then
    CHANNEL="C1"
    MODEL_DETECT="pnn_v2_fasterrcnn_640"
    MODEL_SCORE="pnn_v2_scoring_rank_learning"
elif [ "$MODE" == "PV" ]; then
    CHANNEL="C2"
    MODEL_DETECT="pv_v2_fasterrcnn_640"
    MODEL_SCORE="pv_v2_scoring_rank_learning"
else
    echo "ERROR: MODE must be 'PNN' or 'PV' — got: $MODE"
    exit 1
fi

# Ensure Python output is not buffered (shows progress in real time)
export PYTHONUNBUFFERED=1

echo "========================================"
echo "CPN BATCH DETECTION"
echo "  Mode:    $MODE"
echo "  Channel: $CHANNEL"
echo "  Detect:  $MODEL_DETECT"
echo "  Score:   $MODEL_SCORE"
echo "  Input:   $INPUT_BASE"
echo "  Output:  $OUTPUT_BASE"
echo "  Started: $(date)"
echo "========================================"
echo ""

cd "$CPN_DIR" || { echo "ERROR: CPN_DIR not found: $CPN_DIR"; exit 1; }

# =============================================================================
# PROCESSING LOOP
# =============================================================================
for subdir in "$INPUT_BASE"/*/*/; do
    # Derive relative path (e.g. ADAMTS15/ADAMTS15_1)
    rel_path="${subdir#$INPUT_BASE/}"
    rel_path="${rel_path%/}"

    out_dir="$OUTPUT_BASE/$rel_path"
    mkdir -p "$out_dir"

    # Check that target channel files exist in this folder
    tif_files=("$subdir"*-${CHANNEL}.tif)
    if [ ! -e "${tif_files[0]}" ]; then
        echo "SKIP (no -${CHANNEL}.tif files): $rel_path"
        continue
    fi

    # Skip if already successfully processed
    out_csv="$out_dir/localizations_${CHANNEL}.csv"
    if [ -f "$out_csv" ] && [ -s "$out_csv" ]; then
        echo "SKIP (already done): $rel_path"
        continue
    fi

    echo "========================================"
    echo "Processing: $rel_path"
    echo "Output:     $out_dir"
    echo "Started:    $(date)"
    echo "========================================"

    # Run detection
    python -u predict.py "${MODEL_DETECT}/" \
        -r "${MODEL_SCORE}/" \
        -t "$DETECTION_THRESHOLD" \
        -d "$CUDA_DEVICE" \
        "$subdir"*-${CHANNEL}.tif \
        -o "$out_csv"

    if [ -f "$out_csv" ] && [ -s "$out_csv" ]; then
        echo "Detection OK — $(wc -l < "$out_csv") rows"

        # Render visual output
        python -u draw_predictions.py \
            --root "$subdir" \
            --output "$out_dir/" \
            --mode all \
            "$out_csv"
    else
        echo "ERROR: detection failed or produced empty CSV for $rel_path"
    fi

    echo "Completed: $rel_path at $(date)"
    echo ""
done

echo "========================================"
echo "ALL DONE — $(date)"
echo "========================================"
