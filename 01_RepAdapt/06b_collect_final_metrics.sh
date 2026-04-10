#!/bin/bash

# 06b_collect_final_metrics.sh SAMPLE GENOME IN_BAM OUT_ALIGNMENT OUT_INSERT OUT_INSERT_PDF OUT_WGS OUT_WGS_PDF

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
mkdir -p "$OUTDIR"

picard CollectAlignmentSummaryMetrics \
  R="${GENOME}" \
  I="${IN_BAM}" \
  O="${OUT_ALIGNMENT}"

picard CollectInsertSizeMetrics \
  I="${IN_BAM}" \
  O="${OUT_INSERT}" \
  H="${OUT_INSERT_PDF}"

picard CollectWgsMetricsWithNonZeroCoverage \
  R="${GENOME}" \
  I="${IN_BAM}" \
  O="${OUT_WGS}" \
  CHART="${OUT_WGS_PDF}"
