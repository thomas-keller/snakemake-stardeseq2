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
		expand(["results/diffexp/{contrast}/diffexpr.tsv",
				"results/diffexp/{contrast}/ma-plot.svg"]
				contrast=config["diffexp"]["contrasts"]),
				"results/pca.svg"


def get_fastq(wildcards):
	return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

rule trim_galore:
    input:  get_fastq
    output:
		fastq1="trimmed/{sample}-{unit}.1.fastq.gz",
		fastq2="trimmed/{sample}-{unit}.2.fastq.gz"
		qc="trimmed/{sample}-{unit}.qc.txt"
    log:    "logs/trim_galore/{sample}-{unit}.fastqc"
    threads: 1
    params : jobname = "{sample}-{unit}"
    message: "fastqc {input}: {threads}"
    shell:
        """
        trim_galore input[0] input[1] -o "trimmed/" 2> {log}
        """

rule align_star:
	input:
		"trimmed/{sample}-{unit}.1.fastq.gz",
		"trimmed/{sample}-{unit}.2.fastq.gz",

	log: "logs/star/{sample}-{unit}.star"

	shell:
		"""
		STAR --genomeDir config['ref']['index'] --sjdbGTFfile config['ref']['annotation'] --outSAMtype BAM Unsorted --runThreadN 20 --readFilesIn {input[0]} {input[1]} --outFileNamePrefix {wildcard.sample}-{wildcard-unit}-
		"""

rule feature_counts:
	input:"{sample}-{unit}-Aligned.out.bam




