###This script has the parameters that calls the wgs_proc_RMS file to process WGS GATK data

##Set paths and define metadata 
meta <- read.delim("./rms_wgs.txt", header = T, sep = "\t")
basedir <- "./RMS_2021/"
files.dir <- "./files_analysis/"
seg.path <- paste0(basedir, "./seg/")
vcf.path <- paste0(basedir, "./vcf_filtered/")
outdir <- paste0(basedir, "./results_genome/")


##Read annotation files
##Add bed annotation file for CNAs
refseq.sorted <- readRDS(paste0(files.dir, "refseq.sorted.RDS"))
clinvar <- readRDS(paste0(files.dir, "clinvar.RDS"))
cancer.genes <- readRDS(paste0(files.dir, "cancer.genes.RDS"))
sarcoma.genes <- readRDS(paste0(files.dir, "sarcoma.genes.RDS"))
cosmic <- readRDS(paste0(files.dir, "cmc_v92_V2.RDS"))
cosmic.genes <- readRDS(paste0(files.dir, "cosmic_genes.RDS"))
biomart <- readRDS(paste0(files.dir, "biomart.RDS"))
fusions.bl <- readRDS(paste0(files.dir, "fusions.bl.RDS"))
drivers <- readRDS(paste0(files.dir, "CG_drivers.RDS"))

##Source the file
source("./wgs_proc_RMS.R")

##Run for all samples
cnv <- list()
snv.proc <- list()
snv <- list()
fusion <- list()
for(sample in unique(meta$sample)){
  cnv[[sample]] <- cnv.analysis(sample, meta, outdir, T, c("tumor", "tumoroid"), F)
  snv.proc[[sample]] <- snv.processing(sample, meta, vcf.path)
  snv[[sample]] <- snv.analysis(snv.proc[[sample]])
  gc()
}

