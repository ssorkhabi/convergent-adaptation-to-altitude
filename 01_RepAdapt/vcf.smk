# RepAdapt – Variant Calling (Steps 7–9)

import os
import pandas as pd
from pathlib import Path
import sys

# CONFIG
species = config.get("species")
if species is None:
    raise ValueError("You must provide --config species=Hera (or other)")

print(f"Species set to: {species}", file=sys.stderr)

# PATHS
ROOT_DIR = Path("/rds/project/rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution")
PROJECT_DIR = ROOT_DIR / species
PIPE_DIR = PROJECT_DIR / "01_scripts"

INFO_DIR = PROJECT_DIR / "02_info_files"
GENOME_DIR = PROJECT_DIR / "03_genome"
BAM_DIR = PROJECT_DIR / "06_bam_files"

VCF_DIR = PROJECT_DIR / "07_raw_VCFs"
FILTERED_DIR = PROJECT_DIR / "08_filtered_VCFs"
FINAL_DIR = PROJECT_DIR / "09_final_vcf"
LOG_DIR = PROJECT_DIR / "98_log_files"
DEPTH_DIR = PROJECT_DIR / "97_Local_Depth"

GENOME = GENOME_DIR / f"{species}.fasta"
FAI = Path(str(GENOME) + ".fai")

BAMMAP = INFO_DIR / "bammap.txt"
PLOIDY = INFO_DIR / "ploidymap.txt"

# SAMPLE NAMES (from datatable.txt, column #SRA)
samples = (
    pd.read_csv(
        INFO_DIR / "datatable.txt",
        sep="\t",
        comment="#",   # skip header
        header=None,
        usecols=[0]
    )
    .iloc[:, 0]
    .dropna()        # DROP NaN FIRST
    .astype(str)
    .str.strip()
    .replace("", pd.NA)
    .dropna()
    .unique()
    .tolist()
)

print("SAMPLES (from datatable):", samples, file=sys.stderr)

# CHECK DIRS (and create if missing)
for d in [VCF_DIR, FILTERED_DIR, FINAL_DIR, LOG_DIR]:
    os.makedirs(d, exist_ok=True)

# SANITY CHECKS
for f in [GENOME, FAI, BAMMAP, PLOIDY]:
    if not f.exists():
        raise FileNotFoundError(f"Required file missing: {f}")

# SCAFFOLDS (from genome index)
scaffolds = pd.read_csv(
    FAI, sep="\t", header=None, usecols=[0]
)[0].tolist()

# SHELL
shell.executable("/bin/bash")

# FINAL TARGET
rule all:
    input:
        FINAL_DIR / f"{species}_full_concatenated.vcf.gz",
        FINAL_DIR / f"{species}_full_concatenated_maf01.vcf.gz",
        FINAL_DIR / f"{species}_full_concatenated_maf05.vcf.gz",
        FINAL_DIR / f"{species}_combined_windows.tsv",
        FINAL_DIR / f"{species}_combined_genes.tsv",
        FINAL_DIR / f"{species}_combined_wg.tsv"

############################################################
# RULE 7 — mpileup per scaffold
############################################################
rule mpileup:
    input:
        genome = GENOME,
        bammap = BAMMAP,
        ploidy = PLOIDY
    output:
        vcf = VCF_DIR / "{scaffold}.vcf.gz"
    threads: 4
    conda:
        PIPE_DIR / "vcf.yml"
    log:
        LOG_DIR / "{scaffold}_mpileup.log"
    shell:
        r"""
        bash {PIPE_DIR}/07_mpileup.sh \
            {wildcards.scaffold:q} \
            {threads} \
            {input.genome} \
            {input.bammap} \
            {input.ploidy} \
            {output.vcf:q} \
            &> {log:q}
        """


############################################################
# RULE 8 — Filter scaffold VCFs
############################################################
rule filter_vcfs:
    input:
        vcf = VCF_DIR / "{scaffold}.vcf.gz"
    output:
        vcf = FILTERED_DIR / "{scaffold}.filtered.vcf.gz"
    threads: 8
    conda:
        PIPE_DIR / "vcf.yml"
    log:
        LOG_DIR / "{scaffold}_filter.log"
    shell:
        r"""
        bash {PIPE_DIR}/08_scaffoldVCF_filtering.sh \
            {input.vcf:q} \
            {output.vcf:q} \
            &> {log:q}
        """


############################################################
# RULE 9a — Concatenate parts
############################################################
rule concat_parts:
    input:
        vcfs = expand(
            FILTERED_DIR / "{scaffold}.filtered.vcf.gz",
            scaffold=scaffolds
        )
    output:
        directory(FINAL_DIR / f"{species}_tmp/parts")
    threads: 1
    conda:
        PIPE_DIR / "vcf.yml"
    log:
        LOG_DIR / "concat_parts.log"
    shell:
        r"""
        bash {PIPE_DIR}/09_concat_VCFs.sh \
            {input.vcfs:q} \
            {FINAL_DIR}/{species} \
            parts \
            &> {log:q}
        """

############################################################
# RULE 9b — Build final concatenated VCF
############################################################
rule concat_final:
    input:
        parts = FINAL_DIR / f"{species}_tmp/parts"
    output:
        FINAL_DIR / f"{species}_full_concatenated.vcf.gz"
    threads: 1
    conda:
        PIPE_DIR / "vcf.yml"
    log:
        LOG_DIR / "concat_final.log"
    shell:
        r"""
        bash {PIPE_DIR}/09_concat_VCFs.sh \
            {FINAL_DIR}/{species} \
            final \
            &> {log:q}
        """


############################################################
# RULE 9c — MAF filtering
############################################################
rule maf_filters:
    input:
        FINAL_DIR / f"{species}_full_concatenated.vcf.gz"
    output:
        FINAL_DIR / f"{species}_full_concatenated_maf01.vcf.gz",
        FINAL_DIR / f"{species}_full_concatenated_maf05.vcf.gz"
    threads: 1
    conda:
        PIPE_DIR / "vcf.yml"
    log:
        LOG_DIR / "maf_filters.log"
    shell:
        r"""
        bash {PIPE_DIR}/09_concat_VCFs.sh \
            {input:q} \
            {FINAL_DIR}/{species} \
            maf \
            &> {log:q}
        """

############################################################
# RULE 9d — Combine windowed / gene / WG depth tables
############################################################
rule combine_depth_tables:
    input:
        windows = expand(
            DEPTH_DIR / "{sample}-windows.sorted.tsv",
            sample=samples
        ),
        genes = expand(
            DEPTH_DIR / "{sample}-genes.sorted.tsv",
            sample=samples
        ),
        wg = expand(
            DEPTH_DIR / "{sample}-wg.txt",
            sample=samples
        )
    output:
        windows = FINAL_DIR / f"{species}_combined_windows.tsv",
        genes   = FINAL_DIR / f"{species}_combined_genes.tsv",
        wg      = FINAL_DIR / f"{species}_combined_wg.tsv"
    threads: 1
    conda:
        PIPE_DIR / "vcf.yml"
    log:
        LOG_DIR / "combine_depth_tables.log"
    shell:
        r"""
        bash {PIPE_DIR}/09d_combine_depth_tables.sh \
            {species} \
            {INFO_DIR} \
            {DEPTH_DIR} \
            {FINAL_DIR} \
            &> {log:q}
        """
