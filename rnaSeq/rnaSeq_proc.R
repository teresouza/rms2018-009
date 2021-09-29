##Script to process RNA-seq RMS datasets
meta.rseq <- read.csv("./rnaseq_RMS.csv", header = T)
basedir <- "./rnaseq_data/"


##Function to transform read counts to transcripts per million 9TOM 
tpm <- function(counts, lengths) {
  rate <- counts /(lengths/1000)
  rate / sum(rate) * 1e6
}

in.files <- list.files(basedir, pattern = "*.RNA-Seq.gene_id.exon.counts.txt")
rseq <- list()
for(id in unique(meta.rseq$RNA)){
  all.file <- readLines(paste0(basedir, in.files[grep(id, in.files)])) 
  rseq[[row.names(meta.rseq)[meta.rseq$RNA == id]]] <- read.delim(paste0(basedir, in.files[grep(id, in.files)]), skip = grep("GeneID", all.file)-1, header = T) %>%
    dplyr::select("GeneID", "GeneName", "Counts", "CPM", "Length") %>% drop_na() %>% 
    add_column("log2CPM" = log((.$CPM)+1, base = 2)) %>%
    add_column("TPM" = tpm(.$Counts, .$Length)) %>% add_column("log2TPM" = log(.$TPM+1, base = 2)) %>%
    add_column("EnsemblGene" =gsub("\\..*","", .$GeneID)) %>% 
    add_column("sample" = row.names(meta.rseq)[meta.rseq$RNA == id]) %>% 
    add_column("PMCID" = meta.rseq$PMCID[meta.rseq$RNA == id])
}

#Flatten list to a data frame with all samples
rseq.df2 <- reshape2::dcast(rseq.df,  GeneName + EnsemblGene ~ sample, value.var = "log2TPM", fun.aggregate = sum) %>% 
  filter(rowMeans(.[,-c(1,2)]) >= 2) 
row.names(rseq.df2) <- make.names(rseq.df2$GeneName, unique = T)
rseq.df2 <- rseq.df2[,-c(1,2)]

save.image("./RMS_kidney_TPM.RData")
