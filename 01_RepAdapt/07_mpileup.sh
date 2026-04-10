#!/bin/bash

# 07_mpileup.sh SCAFFOLD THREADS GENOME BAMMAP PLOIDY OUTVCF

set -euo pipefail

SCAFFOLD=$1
THREADS=$2
GENOME=$3
BAMMAP=$4
PLOIDY=$5
OUTVCF=$6

bcftools mpileup \
    --threads ${THREADS} \
    -Ou \
    -f "${GENOME}" \
    --bam-list "${BAMMAP}" \
    -q 5 \
    -r "${SCAFFOLD}" \
    -I \
    -a FMT/AD,FMT/DP \
| bcftools call \
    -S "${PLOIDY}" \
    -G - \
    -a GQ \
    -mv \
    -Oz \
    > "${OUTVCF}"

bcftools index -t "${OUTVCF}"