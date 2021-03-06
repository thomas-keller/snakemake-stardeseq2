import json
import pandas as pd
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
        "fcounts/gene_counts.txt",
        "dex_deseq_star.csv",
        expand("star/{sample}/Aligned.out.bam",sample=samples['sample'].tolist()),
        ['meansd.png','sample_assoc.png','pcaplot.png','maplot.png','anno_heatmap.png']
        #expand(["results/diffexp/{contrast}/diffexpr.tsv",
        #		"results/diffexp/{contrast}/ma-plot.svg"]
        #		contrast=config["diffexp"]["contrasts"]),
        #		"results/pca.svg"


def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

rule trim_galore:
    input: 
        "{sample}_1.fastq","{sample}_2.fastq"
    output:
        fastq1="trimmed/{sample}_1_val_1.fq",
        qc1="trimmed/{sample}_1.fastq_trimming_report.txt",
        fastq2="trimmed/{sample}_2_val_2.fq",
        qc2="trimmed/{sample}_2.fastq_trimming_report.txt"
    log:    "logs/trim_galore/{sample}.fastqc"
    threads: 1
    params : jobname = "{sample}"
    message: "fastqc {input}: {threads}"
    shell:
        """
        trim_galore --paired -o trimmed/ {input[0]} {input[1]}  2> {log}
        """

rule align_star:
    input:
        "trimmed/{sample}_1_val_1.fq",
        "trimmed/{sample}_2_val_2.fq"
    output:
        "star/{sample}/Aligned.out.bam",
        "star/{sample}/ReadsPerGene.out.tab"
    log: "logs/star/{sample}.star"

    shell:
        """
        module purge
        module add apps/star
        STAR --genomeDir {config[ref][index]} --sjdbGTFfile {config[ref][annotation]} --outSAMtype BAM Unsorted --quantMode GeneCounts --runThreadN 12 --readFilesIn {input[0]} {input[1]} --outFileNamePrefix star/{wildcards.sample}/ 2> {log}
        """

rule feature_counts:
    input:
        expand("star/{sample}/Aligned.out.bam",sample=samples['sample'].tolist())
    output:"fcounts/gene_counts.txt"
    log: "logs/feature_counts/fcounts"
    shell:
        """
        module add apps/subRead
        featureCounts -T 5 -p -t exon -g gene_id -a {config[ref][annotation]} -o {output} {input} &> {log}
        """

rule multiqc:
    input:
        expand("star/{sample}/Log.final.out",sample=samples['sample'].tolist()),
    output:
        "multiqc_report.html"
    shell:
        """
        multiqc -f ./trimmed {input} fcounts/
        """
        
rule deseq:
    input:
        "fcounts/gene_counts.txt"
    output:
        "dex_deseq_sig.tsv","meansd.png","maplot.png","sample_assoc.png","pcaplot_sep.png","anno_heatmap.png"
        #report("sample_assoc.png,caption="report/heatmap.rst",category="PCA"),
        #report("pcaplot_sep.png,caption="report/pca.rst",category="PCA"),
        #report('maplot.png',caption="report/maplot.rst",category="maplot"),
        #report('anno_heatmap.png",caption="report/annheatmap.rst",category="ann heat" ),
        #report("dex_deseq_sig.tsv",caption="report/siggenes.rst",category="sig gene"),
    log: "logs/deseq2/deseq.txt"
    shell:
        """
        source deactivate
        source activate root
        module purge
        module add apps/R/3.5.0
        Rscript deseq_star.R &> {log}
        """
