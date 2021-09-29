##Graphs for snvs
library(cowplot)
library(ggplot2)
library(plyr)

meta <- read.delim("./rms_wgs.txt", header = T, sep = "\t")
files.snvs <- list.files("./results_genome/", 
                         pattern = "*_non_syn_mutect2_SNVs_filt.txt")
#Create a list with all snvs identfied
snv <- list()
for(filenames in files.snvs){
  sample <- gsub("_non_syn_mutect2_SNVs_filt.txt", "", filenames)
  snv[[sample]] <- read.delim(paste0("./results_genome/",
                                     filenames), header = T, sep = "\t")
}

##Collate and clean all data
snvs.df <- ldply(snv, data.frame) %>% 
  add_column(cancer = ifelse(.$gene %in% drivers$SYMBOL, "CG", NA)) %>%
  add_column(id = meta$sample2[match(.$sample, meta$id)]) %>% 
  left_join(., cosmic, c("change" = "coord")) %>% filter(source %in% c("tumor", "tumoroid")) %>%
  #here filter for late and seccondary passages
         #filter(source != "tumoroidL", source != "tumoroid2")
snvs.df4 <- snvs.df %>% dplyr::filter(chrom != "chrX", chrom != "chrY", source != "normal",
                PolyPhen2 != "benign") %>%
  add_column(cancer = ifelse(.$gene.x %in% drivers$SYMBOL, "CG", NA)) %>%
  ## uset AF threshold below
  group_by(.id, gene.x, change) %>%
  dplyr::filter(any(gt_AF > 0.4)) %>% 
  ungroup(.id, gene.x) %>%
  dplyr::select(sample.name, source, gene.x, type, change, gt_AF, cancer) %>%
  add_column(subtype = meta.rseq$subtype[match(.$sample.name, meta.rseq$id2)])
snvs.df4$gt_AF <- as.numeric(snvs.df4$gt_AF)
snvs.df4$subtype <- factor(snvs.df4$subtype, levels =c("P7F aRMS", "P3F aRMS","FN aRMS","eRMS"))
colnames(snvs.df4) <- c("sample", "source", "Gene", "type", "change", "VAF","cancer", "subtype")
snvs.df4$Gene <- factor(snvs.df4$Gene)
bold.labels <- ifelse(levels(snvs.df4$Gene) %in% drivers$SYMBOL, yes = "bold", no = "plain")
data <- snvs.df4 %>%
  arrange(subtype, sample, source, Gene, type, change, VAF)
data$sample <- factor(data$sample, levels = unique(data$sample))
#start ggplot
#h1 has the matrix with VAF
h1 <- ggplot(data, aes(x = interaction(sample, source, lex.order = TRUE,sep = " "),
                      y = Gene, color = factor(type))) +
  #geom_point(aes(alpha = VAF), shape = 15, size =factor(VAF)) +
  geom_point(aes(size = VAF)) +
  #geom_vline(xintercept = 1:length(data$.id)*2 + 0.5, color = "grey") +
  geom_hline(yintercept = 1:length(data$Gene)+0.5, color = "grey") +
  scale_color_npg() + theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) + coord_flip() +
  scale_x_discrete(position = 'top') +
  scale_y_discrete(position = 'right') +
  scale_size_continuous(breaks=c(0.1, 0.25, 0.75, 0.9)) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, face=bold.labels, size = 11),
        axis.text.y = element_text(face="bold", hjust = 0.2, size = 13),
        legend.text = element_text(size = 16),
        legend.title=element_text(size=16)) + 
  labs(color = "SNV type")
data2 <- data %>% group_by(sample, source) %>% tally() %>% group_by(sample, source) %>% tally() %>%
  add_column(subtype = meta.rseq$subtype[match(.$sample, meta.rseq$id2)]) %>% 
  arrange(subtype)
data2$subtype <- factor(data2$subtype, levels =c("eRMS","FN aRMS", "P3F aRMS","P7F aRMS"))
colors <- c("#139174", "#f19f89", "#7ec7db", "#2f4175")
#h2 has color labels for the samples
h2 <- ggplot(data2)+
  geom_bar(mapping = aes(x = sample, y=n , fill = subtype), 
           stat = "identity", 
           width = 1)+ scale_fill_manual(values = colors) +
  theme_void()+theme(strip.background = element_blank(), strip.text = element_blank())+
  theme(panel.spacing.y = unit(0.02, "mm"), legend.text = element_text(size = 16), 
        legend.title=element_text(size=16))+
  facet_wrap(subtype~sample+source, nrow = nrow(data2),scale = "free_x") +
  labs(fill = "RMS subtype")
legend <- plot_grid(get_legend(h1), get_legend(h2), ncol = 1)
#extract the legends for both h1 and h2
h1 <- h1 + theme(legend.position = "none")
h2 <- h2 + theme(legend.position = "none")
#merge both in a griid side by side
plot <-plot_grid(h2, h1, align = "h", ncol = 2, axis = "tb", 
                 rel_widths = c(0.2, 10), rel_heights = 1)
#save pdf
pdf("./VAF04_all_stability_filt.pdf", width = 20, height =8)
plot_grid(plot, legend, nrow = 1, rel_widths = c(10,2))
dev.off()

