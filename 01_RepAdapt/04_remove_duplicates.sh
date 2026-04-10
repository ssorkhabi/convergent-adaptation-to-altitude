#!/bin/bash

# 04_remove_duplicates.sh SAMPLE THREADS MEM_MB IN_BAM OUT_BAM DUP_METRICS

set -euo pipefail

# variables
SAMPLE="$1"
THREADS="$2"
MEM_MB="$3"
IN_BAM="$4"
OUT_BAM="$5"
DUP_METRICS="$6"

# Temp directory
TMPDIR="${SLURM_TMPDIR:-$(dirname "$OUT_BAM")/tmp_${SAMPLE}}"
mkdir -p "$TMPDIR"

# Java heap
JAVA_HEAP_MB=$(( MEM_MB - 4000 ))
if (( JAVA_HEAP_MB < 2000 )); then
    JAVA_HEAP_MB=2000
fi

JAVA_OPTS="-Xmx${JAVA_HEAP_MB}m -Djava.io.tmpdir=${TMPDIR}"

# Locate Picard JAR from mamba env
if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "[ERROR] CONDA_PREFIX not set. Did you activate the environment?" >&2
    exit 1
fi

PICARD_JAR=$(ls "$CONDA_PREFIX"/share/picard-*/picard.jar 2>/dev/null | head -n 1)

if [[ -z "$PICARD_JAR" ]]; then
    echo "[ERROR] Could not find picard.jar in conda environment" >&2
    exit 1
fi

# MarkDuplicates (long story why we don't use the 'picard' command directly)
java ${JAVA_OPTS} -jar "$PICARD_JAR" MarkDuplicates \
  INPUT="$IN_BAM" \
  OUTPUT="$OUT_BAM" \
  METRICS_FILE="$DUP_METRICS" \
  REMOVE_DUPLICATES=true \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=SILENT \
  TMP_DIR="$TMPDIR"


# Index
samtools index "$OUT_BAM"

# Cleanup
if [[ -s "$OUT_BAM" && -s "${OUT_BAM}.bai" ]]; then
    rm -f "${IN_BAM}" "${IN_BAM}.bai"
    rm -rf "$TMPDIR"
fi
