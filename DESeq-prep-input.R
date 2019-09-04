library(stringr)

data.dir <- file.path("/media/snfrbert/Ayana/0.Input/RNAs/mix")
fastq <- list.files(data.dir, pattern = "*.fastq")
sampleNO <- str_sub(fastq, 1,19)
#condition <- c(rep("dmog", 5), rep("norm", 4), rep(c(rep("dmog", 5), rep("norm", 5)), 2))
#condition <- c(rep("dmog", 5), rep("norm", 4))
condition <- c(rep(c("dmog","norm"), 3))
libraryName = paste(condition, "-", sampleNO, sep= "")
metadata <- data.frame(sampleNO = sampleNO, condition = condition, fastq = fastq, libraryName = libraryName)
bowidx = "/home/snfrbert/Data/genomes/bowtie/hg19.ebwt/hg19"
gtf = "/home/snfrbert/Data/genomes/hg19_2.gtf"
output.dir = "/media/snfrbert/Ayana/7.RNADE/"

tophat.cmd = with(metadata, paste("tophat -G ", gtf, " -p 4 -o ", output.dir, libraryName, " ", bowidx, " ", data.dir, fastq, "\n\n", sep= "") )
sink(file = "tophat_commands.sh", type="output")
cat("#!/bin/bash \n\n")
cat(tophat.cmd)
sink()

sink(file = "sam-commands.sh")
cat("#!/bin/bash \n\n")
cat(paste("cd", output.dir, "\n\n"))
ob = file.path(output.dir, metadata$libraryName, "accepted_hits.bam")
for(i in seq_len(nrow(metadata))) {
	lib = metadata$libraryName[i]
	cat(paste0("samtools sort -o ", lib, "_s.bam ", ob[i]), "\n")
	cat(paste0("samtools index ", lib, "_s.bam"), "\n\n")
}
sink()

HTSeqCommands <- character(nrow(metadata)*2)
for(i in seq_len(nrow(metadata))) {
	lib = metadata$libraryName[i]
 	HTSeqCommands [2*(i-1) + 1] = paste0("samtools view ", lib, "_s.bam", " > ", lib, "_s.sam", "\n")
	HTSeqCommands [2*(i-1) + 2] = paste0("htseq-count ", " -s \'no\'", " ", lib, "_s.sam ", gtf, " > ", output.dir, lib, "_s_no_DESeq.txt", "\n")
}
sink(file = "countDESeq-no.sh", type="output")
cat('#!/bin/bash \n\n')
cat("cd ", output.dir, "\n\n")
cat(HTSeqCommands)
sink()
