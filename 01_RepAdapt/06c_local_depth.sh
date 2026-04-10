#!/bin/bash

# 06c_local_depth.sh SAMPLE BAM GENOME_BED WINDOWS_BED WINDOWS_LIST GENES_BED GENES_LIST OUT_GENES OUT_WINDOWS OUT_WG

# set -euo pipefail

SAMPLE="$1"
BAM="$2"
GENOME_BED="$3"
WINDOWS_BED="$4"
WINDOWS_LIST="$5"
GENES_BED="$6"
GENES_LIST="$7"
OUT_GENES="$8"
OUT_WINDOWS="$9"
OUT_WG="${10}"

SVDIR=$(dirname "$OUT_WG")

samtools depth -aa "$BAM" > "$SVDIR/${SAMPLE}.depth"

# gene depth
awk '{print $1"\t"($2-1)"\t"$2"\t"$3}' "$SVDIR/${SAMPLE}.depth" \
| bedtools map -a "$GENES_BED" -b stdin -c 4 -o mean -null 0 -g "$GENOME_BED" \
| awk -F"\t" '{print $1":"$2"-"$3"\t"$4}' | sort -k1,1 \
> "$SVDIR/${SAMPLE}-genes.tsv"

join -a 1 -e 0 -o '1.1 2.2' -t $'\t' "$GENES_LIST" "$SVDIR/${SAMPLE}-genes.tsv" \
> "$OUT_GENES"

# window depth
awk '{print $1"\t"($2-1)"\t"$2"\t"$3}' "$SVDIR/${SAMPLE}.depth" \
| bedtools map -a "$WINDOWS_BED" -b stdin -c 4 -o mean -null 0 -g "$GENOME_BED" \
| awk -F"\t" '{print $1":"$2"-"$3"\t"$4}' | sort -k1,1 \
> "$SVDIR/${SAMPLE}-windows.tsv"

join -a 1 -e 0 -o '1.1 2.2' -t $'\t' "$WINDOWS_LIST" "$SVDIR/${SAMPLE}-windows.tsv" \
> "$OUT_WINDOWS"

# whole-genome depth
awk '{sum+=$3; count++} END {if(count>0) print sum/count; else print "No data"}' \
"$SVDIR/${SAMPLE}.depth" > "$OUT_WG"

rm -f "$SVDIR/${SAMPLE}.depth"
