library(DESeq2)
library(tximport)

library(BiocParallel)
#can take this from snakemake@cluster or something
register(MulticoreParam(8))

coldata<-read.table("samples.tsv",header=TRUE,row.names="sample")


#deseq2 input and massaging modified slightly from this good stephen turner gist
#https://gist.github.com/stephenturner/f60c1934405c127f09a6
counts<-read.table("fcounts/gene_counts.txt",header=TRUE,skip=1)

#remove first 5 columns to get only count matrix
counts<-counts[,7:ncol(counts)]

colnames(counts)<-rownames(coldata)

library(magrittr)
library(dplyr)
coldata$dex %<>% relevel("untrt")

counts<-as.matrix(counts)



dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~cell+dex)

dds<-DESeq(dds)

print(resultsNames(dds))

print(nrow(dds))



dds <- dds[ rowSums(counts(dds)) > 1, ]

#print(nrow(dds))

res<-results(dds)
res<-lfcShrink(dds,coef="dex_trt_vs_untrt",type="apeglm")
resOrdered<-res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered),
          file="dex_deseq_star.csv")

resSig<-res %>% data.frame() %>% na.omit() %>% dplyr::filter(padj < 0.01)

write.table(resSig,
          file="dex_deseq_sig.tsv")


lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
library("vsn")
library(ggplot2)
p<-meanSdPlot(cts, ranks = FALSE)
ggsave("meansd.png",width=5,height=5)

vsd <- vst(dds, blind = FALSE)

rld <- rlog(dds, blind = FALSE)

library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
         mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
  
colnames(df)[1:2] <- c("x", "y")  

library(ggplot2)
p<-ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  

ggsave('sample_assoc.png',width=5,height=5)


library("pheatmap")
library("RColorBrewer")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png('heatmap.png')
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

dev.off()



p<-plotPCA(vsd,intgroup=c("dex","cell"))
ggsave("pcaplot.png",width=5,height=5)
pcaData<-plotPCA(vsd,intgroup=c("dex","cell"),returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

p<-ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()

ggsave("pcaplot_sep.png",width=5,height=5)


library(apeglm)

res<-lfcShrink(dds,coef="dex_trt_vs_untrt",type="apeglm")

png("maplot.png")
plotMA(res,ylim=c(-5,5))
dev.off()

library(genefilter)

topVarGenes <- (head(order(assay(vsd)),decreasing=TRUE,20))

mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("cell","dex")])

png('anno_heatmap.png')
pheatmap(mat, annotation_col = anno)
dev.off()

library(knitr)

knit2pdf("deseq.Rmd","deseq.tex")


