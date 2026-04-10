#!/bin/bash

# 03_collect_metrics.sh SAMPLE GENOME IN_BAM OUT_ALIGNMENT OUT_INSERT OUT_INSERT_PDF OUT_WGS OUT_WGS_PDF

set -euo pipefail

# variables
SAMPLE="$1"
GENOME="$2"
IN_BAM="$3"
OUT_ALIGNMENT="$4"
OUT_INSERT="$5"
OUT_INSERT_PDF="$6"
OUT_WGS="$7"
OUT_WGS_PDF="$8"

OUTDIR="$(dirname "$OUT_ALIGNMENT")"

# setting a safe temp directory to avoid 'desk full' errors
TMPDIR="${SLURM_TMPDIR:-$(pwd)/tmp_${SAMPLE}}"
mkdir -p "$TMPDIR"

# Export temp dirs for EVERYTHING
export TMPDIR
export TEMP="$TMPDIR"
export TMP="$TMPDIR"

# Java temp
JAVA_OPTS="-Djava.io.tmpdir=$TMPDIR"
export JAVA_OPTS
export JAVA_TOOL_OPTIONS="-Dsamjdk.use_async_io_read_samtools=false \
                          -Dsamjdk.use_async_io_write_samtools=false \
                          -Dsamjdk.use_async_io_write_tribble=false"

picard $JAVA_OPTS CollectAlignmentSummaryMetrics \
  R="${GENOME}" \
  I="${IN_BAM}" \
  O="${OUT_ALIGNMENT}"

picard $JAVA_OPTS CollectInsertSizeMetrics \
  I="${IN_BAM}" \
  O="${OUT_INSERT}" \
  H="${OUT_INSERT_PDF}"

picard $JAVA_OPTS CollectWgsMetricsWithNonZeroCoverage \
  R="${GENOME}" \
  I="${IN_BAM}" \
  O="${OUT_WGS}" \
  CHART="${OUT_WGS_PDF}"