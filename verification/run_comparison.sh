#!/usr/bin/env bash
set -euo pipefail
# Run both conversion pipelines on one MRXS slide and compare pixel-level output.
#
# Usage:
#   ./verification/run_comparison.sh SLIDE_NAME
#
# Example:
#   ./verification/run_comparison.sh 295182_PM_PatientCevi040534-Panel1_UT1.1
#
# Requires: bioformats2raw (v0.12.0), raw2ometiff (v0.10.0) in repo root.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"

SLIDE_NAME="${1:?Usage: $0 SLIDE_NAME}"
MRXS_PATH="$REPO_DIR/MRXS/${SLIDE_NAME}.mrxs"

REF_DIR="$SCRIPT_DIR/.ref"
RESULT_DIR="$SCRIPT_DIR/results"
mkdir -p "$REF_DIR" "$RESULT_DIR"

echo "=== Slide: $SLIDE_NAME ==="

# ---- 1. bioformats2raw + raw2ometiff (reference pipeline) ----
ZARR_DIR="/tmp/zarr_${SLIDE_NAME}"
rm -rf "$ZARR_DIR"
echo "[1/4] bioformats2raw -> Zarr ..."
"$REPO_DIR/bioformats2raw-0.12.0/bin/bioformats2raw" \
    "$MRXS_PATH" "$ZARR_DIR" \
    --no-minmax 2>&1

REF_TIFF="$REF_DIR/${SLIDE_NAME}.ref.ome.tif"
echo "[2/4] raw2ometiff -> reference TIFF ..."
"$REPO_DIR/raw2ometiff-0.10.0/bin/raw2ometiff" \
    "$ZARR_DIR" "$REF_TIFF" 2>&1

# ---- 2. Our converter (test pipeline) ----
TEST_TIFF="$RESULT_DIR/${SLIDE_NAME}.test.ome.tif"
echo "[3/4] mrxs2ometiff.py -> test TIFF ..."
python3 "$REPO_DIR/mrxs2ometiff.py" "$MRXS_PATH" -o "$TEST_TIFF" 2>&1

# ---- 3. Pixel-level comparison ----
echo "[4/4] Comparing ..."
python3 "$SCRIPT_DIR/compare_tiffs.py" "$REF_TIFF" "$TEST_TIFF" \
    --outdir "$RESULT_DIR"

# Cleanup
rm -rf "$ZARR_DIR"

echo "=== Done. Results in $RESULT_DIR ==="
