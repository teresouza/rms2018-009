##Load libraries 
packages <- c("bedr", "stringr", "tidyverse", "gplots", "plyr", "magrittr","VariantAnnotation", 
              "UpSetR", "circlize", "ggsci", "reshape", "gridExtra", "vcfR")
invisible(lapply(packages, require, character.only = TRUE))
source(paste0(files.dir, "wgs_circos_plot.R"))
################ CNV analyses #######################
## Read called segment files, calculate jaccard similarities and generate circos plots from selected sample.names
cnv.analysis <- function(sample.name, meta, out.dir, circos, comparison, chrX){
  #circos: vector with the types of sample.names to create circos plots. Ex (c("tumor", "tumoroid"))
  #cut.ampl and cut.del specify the cutoff to select for gains or deletions
  ids <- meta$id[meta$sample == sample.name]
  in.files <- list.files(seg.path, pattern = "*called.seg")
  message(paste0("Reading segmentation files from sample ", sample.name))
  seg.list <- list()
  seg.sorted <- list()
  annot <- list()
  for(id in ids){
    source <- unique(meta$source[meta$id == id])
    all.file <- readLines(paste0(seg.path, in.files[grep(id, in.files)]))
    file.read <- read.delim(paste0(seg.path, in.files[grep(id, in.files)]), skip = grep("MEAN_LOG2", all.file)-1, header = T)
    file.read$sample.nameID <- id
    #file.read <- file.read[file.read$CALL != 0,]
    colnames(file.read) <- c("chr", "start", "end", "num.probes", "mean.log2", "call", "sample.nameID")
    ##Create segmentation files and annotate them
    seg.list[[source]] <- file.read[,c(7, 1:6)]
    seg.list[[source]]$sample.name <- sample.name
    seg.list[[source]]$source <- source
    seg.list[[source]]$region <- paste0(seg.list[[source]]$chr, ":", seg.list[[source]]$start, "-", seg.list[[source]]$end)
    annot[[source]] <- bedr.join.region(bedr.sort.region(seg.list[[source]][, c(2:4)], check.zero.based = T,check.chr = F,
                                                         method = "lexicographical", engine = "bedtools", verbose = F),
                                        refseq.sorted, report.n.overlap = TRUE, check.chr = F, verbose = F, check.zero.based = F, check.valid = F,
                                        check.sort = F)
    colnames(annot[[source]]) <- c('sample.CHROM', 'sample.POS', 'sample.END',
                                   'chr', 'start', 'end', 'Gene', 'Overlap')
    annot[[source]] <- annot[[source]] %>% add_column("region" = paste0(.$sample.CHROM, ":",.$sample.POS, "-", .$sample.END)) %>%
      add_column("mean.log2" = seg.list[[source]]$mean.log2[match(.$region,  seg.list[[source]]$region)],
                 "call" = seg.list[[source]]$call[match(.$region,  seg.list[[source]]$region)],
                 "sample2" = meta$sample2[meta$id == id]) %>%
      add_column(CN = ifelse(.$call == 0, 1,  round(2^.$mean.log2, digits = 1))) %>% dplyr::select(-c(sample.CHROM, sample.POS, sample.END))
    annot[[source]][,c("start", "end")] <- apply(annot[[source]][,c("start", "end")], 2, as.numeric)
  cnv.out <- ldply(annot, data.frame, .id = "source") #%>% dplyr::filter(ratio >= cut.ampl | ratio <= cut.del)
  #write.table(cnv.out, paste0(outdir, sample.name, "_CNAs_genes.txt"), quote = F, row.names =  F, sep =  "\t")
  if(isTRUE(circos)){
    circos.plot <- list()
    for(source in comparison){
      x <- annot[[source]]
      refseq.sorted2 <- refseq.sorted %>% dplyr::filter(gene %in% setdiff(.$gene, x$Gene)) %>%
        mutate("mean.log2" = 0)  %>%
        dplyr::select(-"gene")
      circos.plot[[source]] <- x[x$mean.log2 <= 5 &
                                   x$mean.log2 >= -2,
                                 c('chr', 'start', 'end', "mean.log2")] %>% #mutate("CN" = 2^.$mean.log2) %>%
        #dplyr::select(-"mean.log2") %>%
        rbind(refseq.sorted2, .) %>% arrange(chr, start)
    }
    pdf(paste0(out.dir, sample.name, "_circos_0.5GL_",  comparison[1], "_", comparison[2], ".pdf"))
    #Circos plots with points
    circos.plot.fun(circos.plot, chrX)
    text(0,0, sample.name)
    #text(0,0, paste0(sample.name, " - Late"))
    circos.clear()
    graphics.off()
  }
  }
}

############# SNPs and indels processing ############################
## Read filtered vcf files, filter and annotate indels and snps
snv.processing <- function(sample.name, meta, vcf.path){
  ids <- meta$id[meta$sample == sample.name]
  message(paste0("...Reading SNV file for ", sample.name))
  snvs.df <- list()
  for(id in ids){
    vcf.read <- read.vcfR(paste0(vcf.path, id, "_filter.vcf"), verbose = FALSE )
    vcf.df <- vcfR2tidy(vcf.read)
    vcf.info <- data.frame(vcf.df[["fix"]]) %>% 
      add_column(change = paste0(.$CHROM, ":", .$POS, "_", .$REF, "/", .$ALT), 
                 chromkey = paste0(.$ChromKey, ":", .$POS)) %>%
      separate(CSQ, into = paste0("col", c(1:20)), sep = "[|]", extra = "drop")
    assign("impact", vcf.info[,which(apply(vcf.info, 2, function(x) any(grepl("MODIFIER", x))))])
    assign("type", vcf.info[,which(apply(vcf.info, 2, function(x) any(grepl("*_variant", x))))])
    assign("SIFT", vcf.info[,which(apply(vcf.info, 2, function(x) any(grepl("tolerated", x))))])
    assign("PolyPhen2", vcf.info[,which(apply(vcf.info, 2, function(x) any(grepl("benign", x))))])
    assign("VarAnnot", vcf.info[,which(apply(vcf.info, 2, function(x) any(grepl("^COSM*",x))))])
    vcf.info <- vcf.info %>% dplyr::select("CHROM", "POS", "REF", "ALT", "col1","col4","col6", 
                    "DP", "change", "chromkey")
    vcf.info <- cbind(vcf.info, impact, type, SIFT, PolyPhen2, VarAnnot)
    colnames(vcf.info) <- c("chrom", "pos", "ref", "alt", "coord", "ensembl", "gene",
                            "DP","change", "chromkey", "impact","type", 
                            "SIFT",  "PolyPhen2", "VarAnnot")
    vcf.gt <- data.frame(vcf.df[["gt"]]) %>% add_column(chromkey = paste0(.$ChromKey, ":", .$POS)) %>% 
      dplyr::select(chromkey, Indiv, gt_AD, gt_AF, gt_GT) %>% 
      add_column(source = ifelse(.$Indiv == id, meta$source[meta$id ==id], "normal")) %>% 
      separate(gt_AD, into =c("AD_Ref", "AD_Alt"), sep = ",")
    vcf.all <- vcf.info %>% full_join(vcf.gt, by = "chromkey") %>%
      add_column(sample.name =sample.name, .before = "coord")
    snvs.df[[id]] <- vcf.all
    rm(list=c("vcf.info","vcf.all","vcf.gt"))
  }
  snvs.indels <- ldply(snvs.df, data.frame, .id = "sample")
}
  

##snvs.indels is a data frame processed in the previous step 
##type can be "intersection" or "difference"
##comparison is a vector that contains the type of comparisons to be performed (e.g., c("tumor", "tumoroid"))
snv.analysis <- function(snvs.indels){
  sample.name <- unique(snvs.indels$sample.name)
  non_syn <- snvs.indels[snvs.indels$type %in% c("missense_variant","stop_gained","stop_lost","start_lost","inframe_insertion", 
                                                 "inframe_deletion", "frameshift_variant"),]
  non_syn$SIFT <- str_replace(non_syn$SIFT,  " *\\(.*?\\) *", "")
  genes.AD <- non_syn$gene[non_syn$source == "normal" & non_syn$AD_Alt !=0]
  non_syn <- non_syn %>% filter(!gene %in% genes.AD)
  #non_syn$coord <- gsub("c(\"", "", fixed = T, non_syn$coord)
  non_syn$PolyPhen2 <- str_replace(non_syn$PolyPhen2,  " *\\(.*?\\) *", "")
  #non_syn$clinvar <- clinvar$INFO[match(non_syn$coord, clinvar$coord)]
  non_syn$sarcoma <- ifelse(non_syn$gene %in% cancer.genes$gene, "CG", NA)
  non_syn$cosmic <- ifelse(non_syn$gene %in% cosmic.genes$Gene.Symbol, "cosmic", NA)
  write.table(non_syn, file = paste0(outdir, sample.name, "_non_syn_mutect2_SNVs_filt.txt"), 
              quote = F, row.names = F, sep = "\t")
  return(non_syn)
}
