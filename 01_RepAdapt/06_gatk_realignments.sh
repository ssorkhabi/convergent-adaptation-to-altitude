#!/bin/bash

# 06_gatk_realignments.sh SAMPLE THREADS GENOME IN_BAM OUT_BAM

set -euo pipefail

SAMPLE="$1"
THREADS="$2"
GENOME="$3"
IN_BAM="$4"
OUT_BAM="$5"

# defining output directory and temporary directory
OUTDIR="$(dirname "$OUT_BAM")"
mkdir -p "$OUTDIR"

INTERVALS="${OUTDIR}/${SAMPLE}.intervals"

TMPDIR="${OUTDIR}/tmp_${SAMPLE}"
mkdir -p "$TMPDIR"

cleanup() {
    rm -rf "$TMPDIR"
}
trap cleanup EXIT

JAVA_MEM="${JAVA_MEM:-6g}"
export JAVA_OPTS="-Xmx${JAVA_MEM} -Djava.io.tmpdir=${TMPDIR}"

# Checking BAM index exists
if [[ ! -f "${IN_BAM}.bai" && ! -f "${IN_BAM%.bam}.bai" ]]; then
    samtools index -@ "$THREADS" "$IN_BAM"
fi

# RealignerTargetCreator
gatk3 $JAVA_OPTS \
    -T RealignerTargetCreator \
    -R "$GENOME" \
    -I "$IN_BAM" \
    -nt "$THREADS" \
    -o "$INTERVALS"

# IndelRealigner
gatk3 $JAVA_OPTS \
    -T IndelRealigner \
    -R "$GENOME" \
    -I "$IN_BAM" \
    -targetIntervals "$INTERVALS" \
    -o "$OUT_BAM"

samtools index -@ "$THREADS" "$OUT_BAM"

if [[ ! -s "$OUT_BAM" || ! -s "${OUT_BAM}.bai" ]]; then
    echo "ERROR: Output BAM or index missing" >&2
    exit 1
fi

echo "Indel realignment complete for $SAMPLE"

if [[ -s "${OUT_BAM}" && -s "${OUT_BAM}.bai" ]]; then
    rm -rf "${TMPDIR}"
    rm -f "${INTERVALS}"
    rm -f "${IN_BAM}" "${IN_BAM}.bai"
fi
