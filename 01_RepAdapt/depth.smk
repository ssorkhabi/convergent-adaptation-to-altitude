# RepAdapt – final metrics + local depth ONLY

import os
import pandas as pd
from pathlib import Path
import sys

# CONFIG: species
species = config.get("species")
if species is None:
    raise ValueError("You must provide --config species=Hera (or other)")

print(f"Species set to: {species}", file=sys.stderr)

# PATHS AND VARIABLES
ROOT_DIR = Path("/rds/project/rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution")
PROJECT_DIR = ROOT_DIR / species
PIPE_DIR = PROJECT_DIR / "01_scripts"

GENOME_DIR = PROJECT_DIR / "03_genome"
INFO_DIR   = PROJECT_DIR / "02_info_files"
BAM_DIR    = PROJECT_DIR / "06_bam_files"
MERGED_METRICS_DIR = PROJECT_DIR / "99_metrics_merged"
SV_DIR     = PROJECT_DIR / "97_Local_Depth"
LOG_DIR    = PROJECT_DIR / "98_log_files"

GENOME = GENOME_DIR / f"{species}.fasta"
FAI = Path(str(GENOME) + ".fai")
DATATABLE = INFO_DIR / "datatable.txt"

GFF = GENOME_DIR / f"{species}.gff"
GFF3 = GENOME_DIR / f"{species}.gff3"

if GFF.exists():
    ANNOTATION = GFF
elif GFF3.exists():
    ANNOTATION = GFF3
else:
    raise FileNotFoundError(
        f"No annotation found: expected {species}.gff or {species}.gff3 in {GENOME_DIR}"
    )

# CREATE DIRS IF NEEDED
for d in [MERGED_METRICS_DIR, SV_DIR, LOG_DIR]:
    os.makedirs(d, exist_ok=True)

# SANITY CHECKS
if not GENOME.exists():
    raise FileNotFoundError(f"Genome not found: {GENOME}")

if not FAI.exists():
    raise FileNotFoundError(f"FASTA index not found: {FAI}")

if not DATATABLE.exists():
    raise FileNotFoundError(f"Datatable not found: {DATATABLE}")

# LOAD SAMPLES
df = pd.read_csv(DATATABLE, sep="\t")
df = df.rename(columns={"#SRA": "SRA"})
df = df.dropna(subset=["SRA"])
SAMPLES = sorted(df["SRA"].astype(str).unique())

# GLOBAL SHELL SETTINGS
shell.executable("/bin/bash")

# FINAL TARGETS
localrules: all

rule all:
    input:
        expand(MERGED_METRICS_DIR / "{sample}.done", sample=SAMPLES),
        expand(SV_DIR / "{sample}-windows.sorted.tsv", sample=SAMPLES),
        expand(SV_DIR / "{sample}-wg.txt", sample=SAMPLES)


# PREP — genome.bed + windows (once per species)
rule prepare_depth_files:
    input:
        fai = FAI,
        gff = ANNOTATION
    output:
        genome_bed   = INFO_DIR / "genome.bed",
        windows_bed  = INFO_DIR / "windows.bed",
        windows_list = INFO_DIR / "windows.list",
        genes_bed    = INFO_DIR / "genes.bed",
        genes_list   = INFO_DIR / "genes.list"
    conda:
        PIPE_DIR / "RepAdapt2.yml"
    log:
        LOG_DIR / "prepare_depth_files.log"
    shell:
        r"""
        # set -euo pipefail

        # genome.bed
        awk '{{print $1"\t"$2}}' {input.fai} > {output.genome_bed}

        # windows.bed (5 kb)
        awk -v w=5000 '{{chr=$1; len=$2;
            for(start=0; start<len; start+=w) {{
                end=((start+w)<len?start+w:len);
                print chr"\t"start"\t"end;
            }}
        }}' {input.fai} > {output.windows_bed}

        # windows.list
        awk -F"\t" '{{print $1":"$2"-"$3}}' {output.windows_bed} \
            | sort -k1,1 > {output.windows_list}

        # genes.bed
        awk '$3=="gene" {{print $1"\t"$4"\t"$5}}' {input.gff} | uniq > {output.genes_bed}

        # sort genes.bed by genome order
        cut -f1 {input.fai} | while read chr; do
            awk -v chr=$chr '$1==chr {{print}}' {output.genes_bed} | sort -k2,2n
        done > {output.genes_bed}.sorted

        mv {output.genes_bed}.sorted {output.genes_bed}

        # genes.list
        awk -F"\t" '{{print $1":"$2"-"$3}}' {output.genes_bed} \
            | sort -k1,1 > {output.genes_list}
        """


# RULE 6b — Final metrics (realigned BAMs)
rule collect_final_metrics:
    input:
        bam = BAM_DIR / "{sample}.realigned.bam",
        genome = GENOME
    output:
        alignment   = MERGED_METRICS_DIR / "{sample}_alignment_metrics.txt",
        insert_size = MERGED_METRICS_DIR / "{sample}_insert_size_metrics.txt",
        insert_pdf  = MERGED_METRICS_DIR / "{sample}_insert_size_histogram.pdf",
        wgs         = MERGED_METRICS_DIR / "{sample}_wgs_metrics.txt",
        wgs_pdf     = MERGED_METRICS_DIR / "{sample}_wgs_metrics.pdf",
        done        = MERGED_METRICS_DIR / "{sample}.done"
    threads: 4
    resources:
        total_cpus = 4
    conda:
        PIPE_DIR / "RepAdapt2.yml"
    log:
        LOG_DIR / "{sample}_final_metrics.log"
    shell:
        r"""
        bash {PIPE_DIR}/06b_collect_final_metrics.sh \
            {wildcards.sample} \
            {input.genome} \
            {input.bam} \
            {output.alignment} \
            {output.insert_size} \
            {output.insert_pdf} \
            {output.wgs} \
            {output.wgs_pdf} \
            &> {log}

        touch {output.done}
        """


# RULE 6c — Local depth
rule local_depth:
    input:
        bam          = BAM_DIR / "{sample}.realigned.bam",
        genome_bed   = INFO_DIR / "genome.bed",
        windows_bed  = INFO_DIR / "windows.bed",
        windows_list = INFO_DIR / "windows.list",
        genes_bed    = INFO_DIR / "genes.bed",
        genes_list   = INFO_DIR / "genes.list"
    output:
        genes_sorted   = SV_DIR / "{sample}-genes.sorted.tsv",
        windows_sorted = SV_DIR / "{sample}-windows.sorted.tsv",
        wg_depth       = SV_DIR / "{sample}-wg.txt"
    threads: 1
    resources:
        total_cpus = 1,
        mem_mb = 10000
    conda:
        PIPE_DIR / "RepAdapt2.yml"
    log:
        LOG_DIR / "{sample}_local_depth.log"
    shell:
        r"""
        bash {PIPE_DIR}/06c_local_depth.sh \
            {wildcards.sample} \
            {input.bam} \
            {input.genome_bed} \
            {input.windows_bed} \
            {input.windows_list} \
            {input.genes_bed} \
            {input.genes_list} \
            {output.genes_sorted} \
            {output.windows_sorted} \
            {output.wg_depth} \
            &> {log}
        """