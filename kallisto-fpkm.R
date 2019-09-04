library(plyr)
library(rhdf5)


countToTpm <- function(counts, effLen)
{
    rate <- counts / effLen
    denom <- sum(rate)
    (rate / denom)* 1e6
}

countToFpkm <- function(counts, effLen)
{
    N <- sum(counts)
    exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

base_dir <- "/media/snfrbert/Ayana/5.RNA-Seq/all"
sample_id <- dir(file.path(base_dir, "lanes")) 
kallisto_files <- file.path(base_dir,"lanes",sample_id,"kallisto", "abundance.h5")
kallisto_dirs <- file.path(base_dir,"lanes",sample_id,"kallisto")

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'grch37.ensembl.org')
bm <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "external_gene_name"), mart = mart)
bm <- dplyr::rename(bm, target_id = ensembl_transcript_id, ext_gene = external_gene_name)
bm <- arrange(bm, target_id)

for(i in 1:length(kallisto_files)) {
	file <- H5Fopen(kallisto_files[[i]])
	fpkm <- countToFpkm(file$est_counts, file$"aux/eff_lengths")
	tpm <- countToTpm(file$est_counts, file$"aux/eff_lengths")
	df <- data.frame(target_id=file$"aux/ids", counts=file$est_counts, fpkm=fpkm, tpm=tpm)
	tfile <- file.path(kallisto_dirs[[i]], "transcript_abundance.csv")
	write.csv(df, file=tfile, row.names=FALSE, quote=FALSE)
	H5close()

	df <- arrange(df, target_id)
	df <- merge(bm, df, by="target_id")
	df <- arrange(df, ext_gene)
	ndf <- data.frame(gene=df$ext_gene, count_sum=df$counts, fpkm=df$fpkm, tpm=df$tpm, l2fpkm=l2fpkm)
	ndf <- aggregate(. ~ gene, ndf, sum)
	l2fpkm <- log2(ndf$fpkm)
	ndf$l2fpkm = l2fpkm

	gfile <- file.path(kallisto_dirs[[i]], "gene_abundance.csv")
	write.csv(ndf, file=gfile, row.names=FALSE, quote=FALSE)

}