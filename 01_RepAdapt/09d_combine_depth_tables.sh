#!/usr/bin/env bash

# 09b_combine_depth_tables.sh DATASET INFO_DIR DEPTH_DIR FINAL_DIR 

set -euo pipefail

# variables
DATASET=$1
INFO_DIR=$2
DEPTH_DIR=$3
FINAL_DIR=$4

# make header row
echo -e "location\t$(tail -n +2 ${INFO_DIR}/datatable.txt | cut -f1 | sort -u | paste -s -d '\t')" \
    > ${DEPTH_DIR}/depthheader.txt

# get sample names
tail -n +2 ${INFO_DIR}/datatable.txt | cut -f1 | sort | uniq > ${INFO_DIR}/samples.txt

# WINDOWED DEPTH
while read samp; do
    cut -f2 ${DEPTH_DIR}/${samp}-windows.sorted.tsv \
        > ${DEPTH_DIR}/${samp}-windows.sorted.depthcol
done < ${INFO_DIR}/samples.txt

paste \
    $(sed "s|^|${DEPTH_DIR}/|" ${INFO_DIR}/samples.txt | sed 's/$/-windows.sorted.tsv/' | head -n 1) \
    $(sed "s|^|${DEPTH_DIR}/|" ${INFO_DIR}/samples.txt | sed 's/$/-windows.sorted.depthcol/' | tail -n +2) \
    > ${DEPTH_DIR}/combined-windows.temp

cat ${DEPTH_DIR}/depthheader.txt \
    ${DEPTH_DIR}/combined-windows.temp \
    > ${FINAL_DIR}/${DATASET}_combined_windows.tsv

# GENE DEPTH
while read samp; do
    cut -f2 ${DEPTH_DIR}/${samp}-genes.sorted.tsv \
        > ${DEPTH_DIR}/${samp}-genes.sorted.depthcol
done < ${INFO_DIR}/samples.txt

paste \
    $(sed "s|^|${DEPTH_DIR}/|" ${INFO_DIR}/samples.txt | sed 's/$/-genes.sorted.tsv/' | head -n 1) \
    $(sed "s|^|${DEPTH_DIR}/|" ${INFO_DIR}/samples.txt | sed 's/$/-genes.sorted.depthcol/' | tail -n +2) \
    > ${DEPTH_DIR}/combined-genes.temp

cat ${DEPTH_DIR}/depthheader.txt \
    ${DEPTH_DIR}/combined-genes.temp \
    > ${FINAL_DIR}/${DATASET}_combined_genes.tsv


# WHOLE-GENOME DEPTH
while read samp; do
    echo -e "${samp}\t$(cat ${DEPTH_DIR}/${samp}-wg.txt)"
done < ${INFO_DIR}/samples.txt \
    > ${FINAL_DIR}/${DATASET}_combined_wg.tsv

# Cleanup
rm -f ${DEPTH_DIR}/depthheader.txt
rm -f ${INFO_DIR}/samples.txt
rm -f ${DEPTH_DIR}/*-windows.sorted.depthcol
rm -f ${DEPTH_DIR}/*-genes.sorted.depthcol
rm -f ${DEPTH_DIR}/combined-windows.temp
rm -f ${DEPTH_DIR}/combined-genes.temp

echo "DONE! Combined depth tables created."
