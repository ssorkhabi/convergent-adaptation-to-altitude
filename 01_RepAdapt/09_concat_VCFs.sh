#!/bin/bash

# 09_concat_VCFs.sh VCF1.vcf.gz VCF2.vcf.gz ... OUTPREFIX MODE

# MODE: parts | final | maf

set -euo pipefail

MODE="${@: -1}"

if [[ "${MODE}" == "parts" ]]; then
    VCFS=( "${@:1:$#-2}" )
else
    VCFS=()
fi

OUTPREFIX="${@: -2:1}"


echo "Mode: ${MODE}"
echo "Number of VCFs: ${#VCFS[@]}"

SUBSET_SIZE=25
TMPDIR="${OUTPREFIX}_tmp"
PARTSDIR="${TMPDIR}/parts"

mkdir -p "${PARTSDIR}"

# MODE: parts  (build chunked part VCFs)
if [[ "${MODE}" == "parts" ]]; then

    printf "%s\n" "${VCFS[@]}" > "${TMPDIR}/vcfs.list"
    split -l "${SUBSET_SIZE}" "${TMPDIR}/vcfs.list" "${TMPDIR}/subset_"

    i=0
    for subset in "${TMPDIR}"/subset_*; do
        i=$((i+1))
        part="${PARTSDIR}/part_${i}.vcf.gz"

        if [[ -f "${part}" && -f "${part}.tbi" ]]; then
            echo "part_${i} exists, skipping"
            continue
        fi

        echo "Building part_${i}"
        bcftools concat \
          -a \
         -Oz \
         -o "${part}" \
         $(cat "${subset}")

        bcftools index -t "${part}"
    done
fi


# MODE: final  (combine existing parts)
if [[ "${MODE}" == "final" ]]; then

    PARTS=( "${PARTSDIR}"/part_*.vcf.gz )

    if [[ "${#PARTS[@]}" -eq 0 ]]; then
        echo "ERROR: no parts found in ${PARTSDIR}" >&2
        exit 1
    fi

    echo "Final concatenation from ${#PARTS[@]} parts"

    bcftools concat \
     -a \
     -Oz \
     -o "${OUTPREFIX}_full_concatenated.vcf.gz" \
     "${PARTS[@]}"

    bcftools index -t "${OUTPREFIX}_full_concatenated.vcf.gz"

    # Only remove TMPDIR if final output exists
    if [[ -f "${OUTPREFIX}_full_concatenated.vcf.gz" ]]; then
    rm -rf "${TMPDIR}"
    fi
fi


# MODE: maf  (MAF filters only)
if [[ "${MODE}" == "maf" ]]; then

    bcftools view \
        --min-af 0.01:minor \
        -Oz \
        -o "${OUTPREFIX}_full_concatenated_maf01.vcf.gz" \
        "${OUTPREFIX}_full_concatenated.vcf.gz"

    bcftools index -t "${OUTPREFIX}_full_concatenated_maf01.vcf.gz"

    bcftools view \
        --min-af 0.05:minor \
        -Oz \
        -o "${OUTPREFIX}_full_concatenated_maf05.vcf.gz" \
        "${OUTPREFIX}_full_concatenated.vcf.gz"

    bcftools index -t "${OUTPREFIX}_full_concatenated_maf05.vcf.gz"
fi

echo "Done (${MODE})"