#!/bin/bash

# 08_scaffoldVCF_filtering.sh INVCF_GZ OUTVCF_GZ

set -euo pipefail

INVCF=$1
OUTVCF=$2

bcftools filter \
    -e 'AC=AN || MQ < 30' \
    -Oz \
    -o "${OUTVCF}" \
    "${INVCF}"

bcftools index -t "${OUTVCF}"