import json
from os.path import join, basename, dirname
from os import getcwd
from subprocess import check_output
import csv

# Globals ---------------------------------------------------------------------
configfile: "config.yml"



logdir = os.path.join(os.getcwd(), "logs/slurm")
os.makedirs(logdir, exist_ok=True)

# Functions -------------------------------------------------------------------

def rstrip(text, suffix):
    # Remove a suffix from a string.
    if not text.endswith(suffix):
        return text
    return text[:len(text)-len(suffix)]

# Rules -----------------------------------------------------------------------

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)


units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels]) # enforce str in index



rule all:
    input:
        expand("fcounts/{sample}_counts.txt",sample=sample['sample'].tolist())
        #expand(["results/diffexp/{contrast}/diffexpr.tsv",
        #		"results/diffexp/{contrast}/ma-plot.svg"]
        #		contrast=config["diffexp"]["contrasts"]),
        #		"results/pca.svg"


def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

rule trim_galore:
    input:  get_fastq
    output:
        fastq1="trimmed/{sample}_1_trimmed.fq",
        fastq2="trimmed/{sample}_2_trimmed.fq"
        qc="trimmed/{sample}_1.fastq_trimming_report.txt","trimmed/{sample}_2.fastq_trimming_report.txt"
    log:    "logs/trim_galore/{sample}.fastqc"
    threads: 1
    params : jobname = "{sample}"
    message: "fastqc {input}: {threads}"
    shell:
        """
        trim_galore input[0] input[1] -o "trimmed/" 2> {log}
        """

rule align_star:
    input:
        "trimmed/{sample}_1_trimmed.fq",
        "trimmed/{sample}_2_trimmed.fq"
    output:
        "star/{sample}/Aligned.out.bam",
        "star/{sample}/ReadsPerGene.out.tab"
    log: "logs/star/{sample}.star"

    shell:
        """
        STAR --genomeDir config['ref']['index'] --sjdbGTFfile config['ref']['annotation'] --outSAMtype BAM Unsorted --runThreadN 20 --readFilesIn {input[0]} {input[1]} --outFileNamePrefix star/{wildcard.sample}/
        """

rule feature_counts:
    input:"star/{sample}/Aligned.out.bam"
    output:"fcounts/{sample}_counts.txt","fcounts/{sample}_results.bam"
    log: "logs/feature_counts/{sample}.fcounts"
    shell:
        ""
        module add apps/subRead
        featureCounts -T 5 -p -t exon -g gene_id -a config['ref]['annotation'] -o fcounts/{sample}_counts.txt fcounts/{sample}_results.bam
        """




