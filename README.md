# snakemake-stardeseq2

This repo is somewhat ill-named, as it has grown to encompass both a star and kallisto alignment option. The before/after parts remain more or less the same (quality control) and using DESeq2.

The driver script is snakemake_sub.sh, you can submit this from a head node on your cluster.

There is some work to be done, but there should be a way to automate the differential expression analysis. That MAY be too much to hope for, we shall see.


