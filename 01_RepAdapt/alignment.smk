# RepAdapt Pipeline (Steps 1–6c)

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
THISDIR = PROJECT_DIR / "01_scripts"

INFO_DIR = PROJECT_DIR / "02_info_files"
GENOME_DIR = PROJECT_DIR / "03_genome"
RAW_DIR = PROJECT_DIR / "04_raw_data"
TRIM_DIR = PROJECT_DIR / "05_trimmed_data"
BAM_DIR = PROJECT_DIR / "06_bam_files"
METRICS_DIR = PROJECT_DIR / "99_metrics"
MERGED_METRICS_DIR = PROJECT_DIR / "99_metrics_merged"
SV_DIR = PROJECT_DIR / "97_Local_Depth"
LOG_DIR = PROJECT_DIR / "98_log_files"

GENOME = GENOME_DIR / f"{species}.fasta"
FAI = Path(str(GENOME) + ".fai")           # e.g. Hera.fasta.fai
DATATABLE = INFO_DIR / "datatable.txt"
# ANNOTATION = GENOME_DIR / "Hera.gff3"

PIPE_DIR = THISDIR

# Checking key directories exist
for d in [TRIM_DIR, BAM_DIR, METRICS_DIR, MERGED_METRICS_DIR, SV_DIR, LOG_DIR, INFO_DIR]:
    os.makedirs(d, exist_ok=True)

# SANITY CHECKS
if not DATATABLE.exists():
    raise FileNotFoundError(f"Datatable not found at: {DATATABLE}")

if not GENOME.exists():
    raise FileNotFoundError(f"Genome not found at: {GENOME}")

# LOAD SAMPLES
datatable_df = pd.read_csv(DATATABLE, sep="\t")
datatable_df = datatable_df.rename(columns={"#SRA": "SRA"})
datatable_df = datatable_df.dropna(subset=["SRA"])
SAMPLES = sorted(datatable_df["SRA"].astype(str).unique().tolist())

# GLOBAL SHELL SETTINGS
shell.executable("/bin/bash")

# FINAL OUTPUTS
localrules: all

rule all:
    input:
        # Final analysis outputs
        expand(MERGED_METRICS_DIR / "{sample}.done", sample=SAMPLES),
        expand(SV_DIR / "{sample}-windows.sorted.tsv", sample=SAMPLES),
        expand(SV_DIR / "{sample}-wg.txt", sample=SAMPLES),

        # Per-sample metrics (force collection as otherwise it will be skipped)
        expand(METRICS_DIR / "{sample}_alignment_metrics.txt", sample=SAMPLES),
        expand(METRICS_DIR / "{sample}_insert_size_metrics.txt", sample=SAMPLES),
        expand(METRICS_DIR / "{sample}_wgs_metrics.txt", sample=SAMPLES)


############################
# RULE 1 — fastp trimming
############################
rule trim_reads:
    input:
        r1 = RAW_DIR / "{sample}_1.fastq.gz",
        r2 = RAW_DIR / "{sample}_2.fastq.gz",
    output:
        r1 = TRIM_DIR / "{sample}.1.trimmed.fastq.gz",
        r2 = TRIM_DIR / "{sample}.2.trimmed.fastq.gz",
    threads: 16
    resources:
        total_cpus = 16
    conda:
        THISDIR / "RepAdapt2.yml"
    log:
        LOG_DIR / "{sample}_fastp.log"
    shell:
        r"""
        bash {PIPE_DIR}/01_fastp.sh \
            {wildcards.sample} \
            {threads} \
            {input.r1} {input.r2} \
            {output.r1} {output.r2} \
            &> {log}
        """


############################
# RULE 2 — BWA alignment
############################
rule align_reads:
    input:
        r1 = TRIM_DIR / "{sample}.1.trimmed.fastq.gz",
        r2 = TRIM_DIR / "{sample}.2.trimmed.fastq.gz",
        genome = GENOME,
    output:
        bam = BAM_DIR / "{sample}.sorted.bam",
    threads: 32
    resources:
        total_cpus = 32
    conda:
        THISDIR / "RepAdapt2.yml"
    log:
        LOG_DIR / "{sample}_bwa.log"
    shell:
        r"""
        bash {PIPE_DIR}/02_bwa_alignments.sh \
            {wildcards.sample} \
            {threads} \
            {input.genome} \
            {input.r1} \
            {input.r2} \
            {output.bam} \
            &> {log}
        """

######################################
# RULE 3 — Collect initial metrics
######################################
rule collect_metrics:
    input:
        bam = BAM_DIR / "{sample}.sorted.bam",
        genome = GENOME,
    output:
        alignment   = METRICS_DIR / "{sample}_alignment_metrics.txt",
        insert_size = METRICS_DIR / "{sample}_insert_size_metrics.txt",
        insert_pdf  = METRICS_DIR / "{sample}_insert_size_histogram.pdf",
        wgs         = METRICS_DIR / "{sample}_wgs_metrics.txt",
        wgs_pdf     = METRICS_DIR / "{sample}_wgs_metrics.pdf",
    threads: 2
    resources:
        total_cpus = 2
    conda:
        THISDIR / "RepAdapt2.yml"
    log:
        LOG_DIR / "{sample}_collect_metrics.log"
    shell:
        r"""
        bash {PIPE_DIR}/03_collect_metrics.sh \
            {wildcards.sample} \
            {input.genome} \
            {input.bam} \
            {output.alignment} \
            {output.insert_size} \
            {output.insert_pdf} \
            {output.wgs} \
            {output.wgs_pdf} \
            &> {log}
        """

################################
# RULE 4 — Deduplicate BAMs
################################
rule remove_duplicates:
    input:
        bam = BAM_DIR / "{sample}.sorted.bam"
    output:
        bam = BAM_DIR / "{sample}.dedup.bam",
        dup_metrics = METRICS_DIR / "{sample}_DUP_metrics.txt"
    threads: 16
    resources:
        total_cpus = 16,
        mem_mb = 32000
    conda:
        THISDIR / "RepAdapt2.yml"
    log:
        LOG_DIR / "{sample}_dedup.log"
    shell:
        r"""
        bash {PIPE_DIR}/04_remove_duplicates.sh \
            {wildcards.sample} \
            {threads} \
            {resources.mem_mb} \
            {input.bam} \
            {output.bam} \
            {output.dup_metrics} \
            &> {log}
        """

################################
# RULE 5 — Add Read Groups
################################
rule change_read_group:
    input:
        bam = BAM_DIR / "{sample}.dedup.bam",
    output:
        bam = BAM_DIR / "{sample}_RG.bam",
    threads: 2
    resources:
        total_cpus = 2
    conda:
        THISDIR / "RepAdapt2.yml"
    log:
        LOG_DIR / "{sample}_RG.log"
    shell:
        r"""
        bash {PIPE_DIR}/05_change_RG.sh \
            {wildcards.sample} \
            {input.bam} \
            {output.bam} \
            &> {log}
        """

########################################################
# RULE 5b — Merge BAMs (effectively rename in our case)
########################################################
rule merge_bams:
    input:
        bam = BAM_DIR / "{sample}_RG.bam",
    output:
        bam = BAM_DIR / "{sample}.merged.bam",
    threads: 4
    resources:
        total_cpus = 4
    conda:
        THISDIR / "RepAdapt2.yml"
    log:
        LOG_DIR / "{sample}_merge.log"
    shell:
        r"""
        bash {PIPE_DIR}/05b_merge_bams.sh \
            {wildcards.sample} \
            {threads} \
            {input.bam} \
            {output.bam} \
            &> {log}
        """

##########################################
# RULE 6 — GATK Indel Realignment
##########################################
rule gatk_realignment:
    input:
        bam = BAM_DIR / "{sample}.merged.bam",
        genome = GENOME
    output:
        bam = BAM_DIR / "{sample}.realigned.bam"
    threads: 8
    resources:
        total_cpus = 8,
        mem_mb = 32000
    conda:
        THISDIR / "gatk3.8.yml"
    log:
        LOG_DIR / "{sample}_gatk.log"
    shell:
        r"""
        bash {PIPE_DIR}/06_gatk_realignments.sh \
            {wildcards.sample} \
            {threads} \
            {input.genome} \
            {input.bam} \
            {output.bam} \
            &> {log}
        """


########################################################
# RULE 6b — Final metrics on realigned BAMs (merged)
########################################################
rule collect_final_metrics:
    input:
        bam = BAM_DIR / "{sample}.realigned.bam",
        genome = GENOME,
    output:
        alignment   = MERGED_METRICS_DIR / "{sample}_alignment_metrics.txt",
        insert_size      = MERGED_METRICS_DIR / "{sample}_insert_size_metrics.txt",
        insert_pdf  = MERGED_METRICS_DIR / "{sample}_insert_size_histogram.pdf",
        wgs         = MERGED_METRICS_DIR / "{sample}_wgs_metrics.txt",
        wgs_pdf     = MERGED_METRICS_DIR / "{sample}_wgs_metrics.pdf",
        done        = MERGED_METRICS_DIR / "{sample}.done",
    threads: 4
    resources:
        total_cpus = 4
    conda:
        THISDIR / "RepAdapt2.yml"
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

################################################################################
# PREP FOR 06c — genome.bed, windows, genes, lists (once per species)
################################################################################
rule prepare_depth_files:
    input:
        genome = GENOME,
        fai = FAI
    output:
        genome_bed   = GENOME_DIR / "genome.bed",
        windows_bed  = GENOME_DIR / "windows.bed",
        windows_list = GENOME_DIR / "windows.list"
    conda:
        THISDIR / "RepAdapt2.yml"
    log:
        LOG_DIR / "prepare_depth_files.log"
    shell:
        r"""
        set -euo pipefail

        # genome.bed (chromosome + length)
        awk '{{print $1"\t"$2}}' {input.fai} > {output.genome_bed}

        # windows.bed (5 kb windows across genome)
        awk -v w=5000 '{{chr=$1; len=$2;
            for(start=0; start<len; start+=w) {{
                end=((start+w)<len?start+w:len);
                print chr"\t"start"\t"end;
            }}
        }}' {input.fai} > {output.windows_bed}

        # windows.list (chr:start-end format)
        awk -F"\t" '{{print $1":"$2"-"$3}}' {output.windows_bed} | sort -k1,1 > {output.windows_list}
        """

################################################################################
# RULE 6c — Local depth (windows, genes & genome depth outputs)
################################################################################
rule local_depth:
    input:
        bam         = BAM_DIR / "{sample}.realigned.bam",
        genome      = GENOME,
        fai         = FAI,
        genome_bed  = GENOME_DIR / "genome.bed",
        windows_bed = GENOME_DIR / "windows.bed",
        windows_list= GENOME_DIR / "windows.list"
    output:
        windows_sorted = SV_DIR / "{sample}-windows.sorted.tsv",
        wg_depth       = SV_DIR / "{sample}-wg.txt"
    threads: 8
    resources:
        total_cpus = 8
    conda:
        THISDIR / "RepAdapt2.yml"
    log:
        LOG_DIR / "{sample}_local_depth.log"
    shell:
        r"""
        bash {PIPE_DIR}/06c_local_depth.sh \
            {wildcards.sample} \
            {threads} \
            {input.bam} \
            {input.fai} \
            {input.genome} \
            {input.genome_bed} \
            {input.windows_bed} \
            {input.windows_list} \
            {output.windows_sorted} \
            {output.wg_depth} \
            &> {log}
        """
