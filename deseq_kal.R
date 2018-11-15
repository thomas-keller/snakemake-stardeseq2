library(DESeq2)
library(tximport)
library(GenomicFeatures)
library(BiocParallel)
#can build this out to replace the stupid snakemake input thing
args <- commandArgs(trailingOnly=TRUE)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
register(MulticoreParam(8))

coldata<-read.table("samples.tsv",header=TRUE,row.names="sample")


files<-file.path("kallisto",rownames(coldata),"abundance.tsv")
#names(files)<-paste0("sample",1:8)
names(files)<-rownames(coldata)

txdb <- makeTxDbFromGFF('/work/t/tekeller/refhg38/gencode.v28.annotation.gtf',dataSource="genecode v28", organism='Homo sapiens')
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

txi <- tximport(files,type='kallisto',tx2gene=tx2gene,ignoreAfterBar=TRUE)

dds <- DESeqDataSetFromTximport(txi, coldata, ~cell+dex)

dds<-DESeq(dds)

res <- results(dds)
resOrdered <- res[order(res$pvalue),]


write.csv(as.data.frame(resOrdered), 
          file="dex_deseq_kal.csv")

