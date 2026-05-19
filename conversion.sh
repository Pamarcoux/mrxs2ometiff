#!/usr/bin/env bash
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate BsAbsPredict

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)" 2>/dev/null || pwd
mkdir -p TIFF
for file in MRXS/*.mrxs; do
    [ -e "$file" ] || continue
    base="$(basename "$file" .mrxs)"
    output="TIFF/${base}.ome.tif"
    echo "Processing: $file"
    python3 "$SCRIPT_DIR/mrxs2ometiff.py" "$file" -o "$output"
    echo "Done: $output"
done
