#!/bin/bash

# 05_change_RG.sh SAMPLE IN_BAM OUT_BAM

set -euo pipefail

# variables
SAMPLE="$1"
IN_BAM="$2"
OUT_BAM="$3"

OUTDIR="$(dirname "$OUT_BAM")"

picard AddOrReplaceReadGroups \
  I="${IN_BAM}" \
  O="${OUT_BAM}" \
  RGID="${SAMPLE}" \
  RGLB="${SAMPLE}_LB" \
  RGPL=ILLUMINA \
  RGPU=unit1 \
  RGSM="${SAMPLE}"

samtools index "${OUT_BAM}"

if [[ -s "${OUT_BAM}" && -s "${OUT_BAM}.bai" ]]; then
    rm -f "${IN_BAM}" "${IN_BAM}.bai"
fi