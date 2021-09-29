meta.rseq <- read.csv("./rnaseq_RMS.csv", header = T)
basedir <- "./seq_data/"
load("./RMS_kidney_TPM.RData")

##Heatmap
library("pheatmap")
#generate metadata
meta.hc <- data.frame(sample = colnames(rseq.df), Entity = factor(c(rep("RMS", times = 6), 
                                                                rep("Kidney", times = 5)), levels = c("RMS", "Kidney")),
                      Source = factor(c(rep(c("Tumor", "Tumoroid"), each = 3), rep("Tumor", times = 5)), 
                                      levels = c("Tumor", "Tumoroid")))
meta.hc <-  meta.rseq[meta.rseq$passage == "Early",]
meta.hc$type <- factor(meta.hc$type, levels=c("CCRCC","CMN","WT","eRMS", "FN aRMS", "P3F aRMS","P7F aRMS"))
meta.hc$source <- factor(meta.hc$source, levels = c("Tumor", "Tumoroid"))
meta.hc$origin <- factor(meta.hc$origin, levels = c("Kidney", "RMS"))
meta.hc <- meta.hc %>% select(sample, type, source,  origin) %>% arrange(origin, source, type)
row.names(meta.hc) <- meta.hc$sample
meta.hc <- meta.hc[,-1]
colnames(meta.hc) <- c("Type", "Source", "Entity")
meta.hc <- meta.hc[,c(1,2,3)]
#assign colors
ann_colors = list(
  Entity = c(RMS = "#f4a71d", Kidney= "#1760a9"),
  Source = c(Tumor="#f4a71d", Tumoroid ="#1760a9"), 
  Type = c(CCRCC = "#7e6148", CN = "#dc1f00", MN = "#3c5488", Nb = "#91d1c2", WT="#d24835",
           eRMS = "#139174", `FN aRMS` ="#f19f89", `P3F aRMS`="#7ec7db", `P7F aRMS` = "#2f4175"))
#calculate correlation
data.cor <- log2tpm[,colnames(log2tpm) %in% row.names(meta.hc)]
data.cor <- data.cor[,row.names(meta.hc)]
pearson.cor <- cor(data.cor)
#plot heatmap correlogram
pdf("dendrogram_kidney_RMS.pdf", height = 10, width = 10)
pheatmap(pearson.cor,annotation_colors =  ann_colors, 
         annotation_col = meta.hc, gaps_row = c(36, 52), gaps_col = c(36, 52),
         cluster_rows = F, cluster_cols =  F)
dev.off()


##Plot for models considering stability across passages (late and secondary). If using the original set, just skip the steps of data filtering
#select models and filter data
models <- c("RMS012", "RMS444", 'RMS335', "RMS006", "RMS007")
meta3 <- meta.rseq[meta.rseq$id2 %in% models,] %>% dplyr::filter(source != "Tumor")
meta.pca <- meta3
data <- log2tpm[,colnames(log2tpm) %in% meta3$sample]
#get principal component vectors of the data
pca.res <-  prcomp(t(data))
pca.resx <- pca.res$x
project.pca.proportionvariances <- ((pca.res$sdev^2) / (sum(pca.res$sdev^2)))*100
#Plot PCA
pca.resx <- merge(pca.resx, meta.pca, by.x = "row.names", by.y = "sample")
pdf("model_stability_kidney_RMS.pdf")
ggplot(pca.resx, aes(x= PC1, y = PC2)) + 
  geom_point(aes(shape = factor(passage), fill = factor(id2), size = 10), stroke =1) + 
  scale_shape_manual(values = c(21, 22, 23)) +
  xlab(paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%")) +
  ylab(paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%")) +
  guides(size = F,
         fill=guide_legend(override.aes=list(shape=21, size = 5)),
         shape = guide_legend(override.aes=list(size = 5))) +
  scale_color_manual(values = "black") +
  theme_few() + theme(legend.text = element_text(size = 14),
                      legend.title = element_text(size=14)) + 
  scale_fill_npg() + 
  labs(fill = "RMS model", shape = 'Passage')
dev.off()

