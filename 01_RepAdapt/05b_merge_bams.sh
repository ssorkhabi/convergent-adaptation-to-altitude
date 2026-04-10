#!/bin/bash

# 05b_merge_bams.sh SAMPLE THREADS IN_BAM OUT_BAM

set -euo pipefail

# variables
SAMPLE="$1"
THREADS="$2"
IN_BAM="$3"
OUT_BAM="$4"

OUTDIR="$(dirname "$OUT_BAM")"

# For now we just "merge" a single BAM (no technical replicates); this is compatible with future multi-BAM merging.
samtools merge -@ "${THREADS}" "${OUT_BAM}" "${IN_BAM}"
samtools index "${OUT_BAM}"

if [[ -s "${OUT_BAM}" && -s "${OUT_BAM}.bai" ]]; then
    rm -f "${IN_BAM}" "${IN_BAM}.bai"
fi