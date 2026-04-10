#!/bin/bash

# 02_bwa_alignments.sh SAMPLE THREADS GENOME IN_R1 IN_R2 OUT_BAM

set -euo pipefail

SAMPLE="$1"
THREADS="$2"
GENOME="$3"
IN_R1="$4"
IN_R2="$5"
OUT_BAM="$6"

OUTDIR="$(dirname "$OUT_BAM")"
mkdir -p "$OUTDIR"

TMPDIR="${SLURM_TMPDIR:-${OUTDIR}/tmp_${SAMPLE}}"
mkdir -p "$TMPDIR"

bwa mem -t "${THREADS}" \
  -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
  "${GENOME}" \
  "${IN_R1}" \
  "${IN_R2}" \
| samtools sort \
    -@ "${THREADS}" \
    -T "${TMPDIR}/${SAMPLE}" \
    -o "${OUT_BAM}"

samtools index "${OUT_BAM}"

if [[ -s "${OUT_BAM}" ]]; then
    rm -f "${IN_R1}" "${IN_R2}"
    rm -rf "${TMPDIR}"
fi

echo "BWA alignment and indexing completed for ${SAMPLE}"