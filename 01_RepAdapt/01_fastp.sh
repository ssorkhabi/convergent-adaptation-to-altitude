#!/bin/bash

# 01_fastp.sh SAMPLE THREADS IN_R1 IN_R2 OUT_R1 OUT_R2

set -euo pipefail

SAMPLE="$1"
THREADS="$2"
IN_R1="$3"
IN_R2="$4"
OUT_R1="$5"
OUT_R2="$6"

OUTDIR="$(dirname "$OUT_R1")"
REPORTDIR="${OUTDIR}/01_reports"
mkdir -p "$REPORTDIR"

fastp \
  -w "${THREADS}" \
  -i "${IN_R1}" \
  -I "${IN_R2}" \
  -o "${OUT_R1}" \
  -O "${OUT_R2}" \
  -j "${REPORTDIR}/${SAMPLE}.json" \
  -h "${REPORTDIR}/${SAMPLE}.html"

echo "Fastp trimming completed for ${SAMPLE}"