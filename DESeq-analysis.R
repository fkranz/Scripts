library(geneplotter)
library(ggplot2)
library(plyr)
library(LSD)
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(stringr)
library(topGO)
library(genefilter)
library(biomaRt)
library(dplyr)
library(EDASeq)
library(fdrtool)

data.dir <- file.path("/media/snfrbert/Ayana/0.Input/RNAs")
fastq <- list.files(data.dir, pattern = "*.fastq")
sampleNO <- str_sub(fastq, 1,19)
condition <- c(rep("dmog", 5), rep("norm", 4), rep(c(rep("dmog", 5), rep("norm", 5)), 2))
#condition <- c(rep("dmog", 5), rep("norm", 4))
#condition <- c(rep(c("dmog","norm"), 3))
libraryName = paste(condition, "-", sampleNO, sep= "")
output.dir = "/media/snfrbert/Ayana/7.RNADE/"


sampleTable <- data.frame(sampleName = libraryName, fileName = paste0(libraryName, "_s_no_DESeq.txt"), condition = condition, sampleNO = sampleNO, fastq = fastq)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory= output.dir, design = ~ condition)
#Filter counts with less then two counts in total
dds <- dds[rowSums(counts(dds)) > 1,]

#add conditions to dds object
dds$condition <- factor(dds$condition, levels=c("norm", "dmog"))

# Quality control
GeneCounts <- counts(dds)
idx.nz <- apply(GeneCounts, 1, function(x) {all(x > 0)})
sum(idx.nz) # should be several thousand

#check some random samples
nz.counts <- subset(GeneCounts, idx.nz)
sam <- sample(dim(nz.counts)[1], 5)
nz.counts[sam, ]


# Differential Expression Analysis
# perfoms: estimation of size factors, estimation of dispersion, model fitting
dds <- DESeq(dds)  
res <- results(dds)

#get some results
sum(res$padj < 0.1, na.rm = TRUE) #how many significantly differentially expressed genes
summary(res)

#create result files with gene names for mrna and ncrna
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'grch37.ensembl.org')
bm_mrna <- biomaRt::getBM(attributes = c("refseq_mrna", "external_gene_name"), filter = "refseq_mrna", values = rownames(res), mart = mart)
bm_mrna <- dplyr::rename(bm_mrna, refseq = refseq_mrna, gene_name = external_gene_name)
bm_mrna <- arrange(bm_mrna, refseq)

bm_ncrna <- biomaRt::getBM(attributes = c("refseq_ncrna", "external_gene_name"), filter = "refseq_ncrna", values = rownames(res), mart = mart)
bm_ncrna <- dplyr::rename(bm_ncrna, refseq = refseq_ncrna, gene_name = external_gene_name)
bm_ncrna <- arrange(bm_ncrna, refseq)

bm_rna <- rbind(bm_mrna, bm_ncrna)

rres <- data.frame(refseq=rownames(res), baseMean=res$baseMean, log2FoldChange=res$log2FoldChange, lfcSE=res$lfcSE, stat=res$stat, pvalue=res$pvalue, padj=res$padj)
rres <- merge(rres, bm_rna, by="refseq")
rres <- rres[order(-rres$padj),]
write.csv(rres, "deseq-lfc.csv", quote=FALSE, rownames=FALSE)

#generate some plots
pdf("multiecdf.pdf")
multiecdf(counts(dds, normalized = T)[idx.nz, ], xlab="mean counts", xlim=c(0,1000))
dev.off()

pdf("multidensity.pdf")
multidensity(counts(dds, normalized = T)[idx.nz, ], xlab="mean counts", xlim=c(0,1000))
dev.off()

pdf("MAplot.pdf")
plotMA(res, main="DESeq2-Result", ylim=c(-3,3))
dev.off()

pdf("pairwiseMAs.pdf")
#check wich combinations make sense!
maidx = t(combn(1:6, 2))	
for(i in 1:15) {
	MDPlot(counts(dds, normalized = T)[idx.nz ,], c(maidx[i,1], maidx[i,2]), main = paste(colnames(dds)[maidx[i,1]], " vs ", colnames(dds)[maidx[i,2]]), ylim = c(-3, 3) )
}
dev.off()

pdf("counts.pdf")
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
dev.off()

rld <- rlogTransformation(dds, blin=TRUE)
pdf("Heatmap.pdf")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13, 13))
dev.off()

pdf("PCAs.pdf")
DESeq2::plotPCA(rld, intgroup=c("condition"))
dev.off()