##Analysis of WGS WES data
##Analysis of copy number data
library(copynumber)
library(tidyverse)
library(stringr)
library(plyr)
library(magrittr)
library(gplots)
seg.path <- "./seg/"
setwd(seg.path)
meta.rseq <- read.delim("./metadata_rseq_rms.txt", 
                        header = T, sep = "\t", row.names = 1)
meta <- read.delim("./rms_wgs.txt", header = T, sep = "\t") 
meta$subtype <- meta.rseq$subtype[match(meta$sample2, row.names(meta.rseq))]
in.files <- list.files(seg.path, pattern = "*called.seg")

#Read files and create a data frame with segments 
segments <- list()
for(source in unique(meta$source)){
  ids <- meta$id[meta$source == source]
  for(id in ids){
    all.file <- readLines(paste0(seg.path, in.files[grep(id, in.files)]))
    file.read <- read.delim(paste0(seg.path, in.files[grep(id, in.files)]), comment.char = "@", header = T)
    file.read$MEAN_LOG2_COPY_RATIO <- ifelse(file.read$CALL != 0,file.read$MEAN_LOG2_COPY_RATIO, 0)
    file.read$sampleID <- meta$sample[meta$id == id]
    file.read$arm <- "p"
    file.read$CONTIG <- gsub("chr", "", file.read$CONTIG)
    file.read <- file.read[,c(7, 1, 8, 2, 3, 4, 5)]
    colnames(file.read) <- c("sampleID", "chrom", "arm", "start.pos", "end.pos", "n.probes", "mean")
    file.read$chrom <- as.numeric(file.read$chrom)
    segments[[source]] <- rbind(segments[[source]], file.read)
  }
}
test <- ldply(segments)
test$type <- meta$subtype[match(test$sampleID, meta$sample)]

#Generate copy number frequency plots
copynumber::plotFreq(segments=test[test$.id %in% c("tumoroid", "tumoroidL", "tumoroid2") & 
                                     test$type %in% 
                         c("P3F aRMS", "P7F aRMS") , -c(1,9)], 
         thres.gain=0.1,thres.loss=-0.1, main = "Tumoroid Fusion Positive")
pdf("/./FreqPlot_0.5GL_Tumoroid_FN.pdf")
copynumber::plotFreq(segments=test[test$.id == "tumoroid" & 
                                     test$type %in% 
                                     c("eRMS", "FN aRMS") , -c(1,9)], 
                     thres.gain=0.32,thres.loss=-0.41, main = "Tumoroid Fusion Negative")
dev.off()
