library(tximport)
library(DESeq2)
library(readr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

#imports
tx2gene <- read_delim("tx2gene_nover.tsv", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)

coldata <- read_delim("metadata.tsv", "\t", 
                      escape_double = FALSE, trim_ws = TRUE)

files <- file.path(coldata$Files)
names(files) <- coldata$Sample

coldata$Treatment <- as.factor(coldata$Treatment)

txi<- tximport(files, type="kallisto",txIn = TRUE, txOut = FALSE, countsFromAbundance="no",
               tx2gene = tx2gene,ignoreTxVersion = TRUE)


dds <- DESeqDataSetFromTximport(txi,
                                colData = coldata,
                                design = ~Treatment)

dds$Treatment <- relevel(dds$Treatment, ref = "Control")

##

dds <- DESeq(dds)
resultsNames(dds)

vst_dds <- vst(dds)

#PCA
plotPCA(vst_dds, intgroup="Treatment")
ggsave(filename = "bulk_iMGLs_PCA.png",device = "png")
pca_data <- plotPCA(vst_dds, intgroup="Treatment",returnData=TRUE)

#DEG spreadsheets
res_alpha_ctrl <- results(dds, contrast=c("Treatment","Alpha","Control"))
res_alpha_ctrl_ord <- subset(res_alpha_ctrl[order(res_alpha_ctrl$log2FoldChange,decreasing = TRUE),], is.na(log2FoldChange) == FALSE & padj < 0.01)
write.table(res_alpha_ctrl_ord,"alpha_vs_ctrl_DEGs.tsv",quote=FALSE,sep="\t")
write.table(res_alpha_ctrl,"alpha_vs_ctrl_all.tsv",quote=FALSE,sep="\t")

res_beta_ctrl <- results(dds, contrast=c("Treatment","Beta","Control"))
res_beta_ctrl_ord <- subset(res_beta_ctrl[order(res_beta_ctrl$log2FoldChange,decreasing = TRUE),], is.na(log2FoldChange) == FALSE & padj < 0.01)
write.table(res_beta_ctrl_ord,"beta_vs_ctrl_DEGs.tsv",quote=FALSE,sep="\t")
write.table(res_beta_ctrl,"beta_vs_ctrl_all.tsv",quote=FALSE,sep="\t")

res_gamma_ctrl <- results(dds, contrast=c("Treatment","Gamma","Control"))
res_gamma_ctrl_ord <- subset(res_gamma_ctrl[order(res_gamma_ctrl$log2FoldChange,decreasing = TRUE),], is.na(log2FoldChange) == FALSE & padj < 0.01)
write.table(res_gamma_ctrl_ord,"gamma_vs_ctrl_DEGs.tsv",quote=FALSE,sep="\t")
write.table(res_gamma_ctrl,"gamma_vs_ctrl_all.tsv",quote=FALSE,sep="\t")

#top up and down per comparison
alpha_up <- rownames(res_alpha_ctrl_ord[c(1:30),])
beta_up <- rownames(res_beta_ctrl_ord[c(1:30),])
gamma_up <- rownames(res_gamma_ctrl_ord[c(1:30),])

alpha_down <- rownames(tail(res_alpha_ctrl_ord,n=30))
beta_down <- rownames(tail(res_beta_ctrl_ord,n=30))
gamma_down <- rownames(tail(res_gamma_ctrl_ord,n=30))


all_up <- unique(c(alpha_up,beta_up,gamma_up))
all_down <- unique(c(alpha_down,beta_down,gamma_down))


sampleDists <- dist(t(assay(vst_dds.t0)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vst_dds.t0$Sample
colnames(sampleDistMatrix) <- vst_dds.t0$Sample
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,col=colors)


library(pheatmap)
ann_df <- as.data.frame(colData(dds)[,c("Treatment")])
colnames(ann_df) <- "Treatment"
rownames(ann_df) <- rownames(colData(dds))
top30_up <- subset(assay(vst_dds),rownames(assay(vst_dds)) %in% all_up)
top30_down <- subset(assay(vst_dds),rownames(assay(vst_dds)) %in% all_down)

#alpha/control colors flipped vs pca
library(scales)
ggplot_first4 <- hue_pal()(4) 

ann_colors <- list(Treatment=c(Control=ggplot_first4[1], Alpha=ggplot_first4[2], Beta=ggplot_first4[3], Gamma=ggplot_first4[4]))

pheatmap(top30_up,annotation_col=ann_df,annotation_colors = ann_colors)
#750x900
pheatmap(top30_down,annotation_col=ann_df,annotation_colors = ann_colors)


save(alpha_up,alpha_down,beta_up,beta_down,gamma_up,gamma_down,file="DEG_sets.Rdata")

#volcanos

library(tidyverse)
library(ggrepel)

alpha_volc <- data.frame(Gene = rownames(res_alpha_ctrl),logFC = res_alpha_ctrl$log2FoldChange, padj = res_alpha_ctrl$padj)

alpha_volc <- alpha_volc[is.finite(alpha_volc$padj),]

alpha_volc$signif <- ifelse(alpha_volc$padj < 0.01, 
                         ifelse(alpha_volc$logFC <= -2, 
                                "sig_low", 
                                ifelse(alpha_volc$logFC >= 2, 
                                       "sig_high",
                                       "sig_mid")), 
                         "NS")



alpha_volcplot <- ggplot(alpha_volc, aes(x = logFC, y = -log10(padj))) +
  geom_point(aes(color=signif), size=3) +
  scale_color_manual(values = c("sig_high" = "red", "sig_low" = "blue", "sig_mid" = "black", "NS" = "grey"),guide=FALSE) +
  theme_minimal()+ 
  #theme(legend.position=c(.9,.8),panel.grid = element_blank()) +
  theme(panel.grid = element_blank(),legend.key = element_blank()) +
  #geom_text_repel(
  #  data = subset(alpha_volc, signif %in% c("sig_low","sig_high")), aes(label = Gene), size = 3.5, segment.color = "gray60",segment.alpha = 0.3,
  #  box.padding = unit(0.35, "lines"), point.padding = unit(1, "lines"),min.segment.length = 0.2,max.overlaps = 100,bg.color = "white") +
  labs(x="log2(fold change)",y="-log10(p-value)",color="Significance",title="alpha vs. ctrl") 

ggsave(filename="alpha_volc.png",device="png")


beta_volc <- data.frame(Gene = rownames(res_beta_ctrl),logFC = res_beta_ctrl$log2FoldChange, padj = res_beta_ctrl$padj)

beta_volc <- beta_volc[is.finite(beta_volc$padj),]

beta_volc$signif <- ifelse(beta_volc$padj < 0.01, 
                            ifelse(beta_volc$logFC <= -2, 
                                   "sig_low", 
                                   ifelse(beta_volc$logFC >= 2, 
                                          "sig_high",
                                          "sig_mid")), 
                            "NS")



beta_volcplot <- ggplot(beta_volc, aes(x = logFC, y = -log10(padj))) +
  geom_point(aes(color=signif), size=3) +
  scale_color_manual(values = c("sig_high" = "red", "sig_low" = "blue", "sig_mid" = "black", "NS" = "grey"),guide=FALSE) +
  theme_minimal()+ 
  #theme(legend.position=c(.9,.8),panel.grid = element_blank()) +
  theme(panel.grid = element_blank(),legend.key = element_blank()) +
  #geom_text_repel(
  #  data = subset(beta_volc, signif %in% c("sig_low","sig_high")), aes(label = Gene), size = 3.5, segment.color = "gray60",segment.beta = 0.3,
  #  box.padding = unit(0.35, "lines"), point.padding = unit(1, "lines"),min.segment.length = 0.2,max.overlaps = 100,bg.color = "white") +
  labs(x="log2(fold change)",y="-log10(p-value)",color="Significance",title="beta vs. ctrl") 

ggsave(filename="beta_volc.png",device="png")

gamma_volc <- data.frame(Gene = rownames(res_gamma_ctrl),logFC = res_gamma_ctrl$log2FoldChange, padj = res_gamma_ctrl$padj)

gamma_volc <- gamma_volc[is.finite(gamma_volc$padj),]

gamma_volc$signif <- ifelse(gamma_volc$padj < 0.01, 
                            ifelse(gamma_volc$logFC <= -2, 
                                   "sig_low", 
                                   ifelse(gamma_volc$logFC >= 2, 
                                          "sig_high",
                                          "sig_mid")), 
                            "NS")



gamma_volcplot <- ggplot(gamma_volc, aes(x = logFC, y = -log10(padj))) +
  geom_point(aes(color=signif), size=3) +
  scale_color_manual(values = c("sig_high" = "red", "sig_low" = "blue", "sig_mid" = "black", "NS" = "grey"),guide=FALSE) +
  theme_minimal()+ 
  #theme(legend.position=c(.9,.8),panel.grid = element_blank()) +
  theme(panel.grid = element_blank(),legend.key = element_blank()) +
  #geom_text_repel(
  #  data = subset(gamma_volc, signif %in% c("sig_low","sig_high")), aes(label = Gene), size = 3.5, segment.color = "gray60",segment.gamma = 0.3,
  #  box.padding = unit(0.35, "lines"), point.padding = unit(1, "lines"),min.segment.length = 0.2,max.overlaps = 100,bg.color = "white") +
  labs(x="log2(fold change)",y="-log10(p-value)",color="Significance",title="gamma vs. ctrl") 

ggsave(filename="gamma_volc.png",device="png")
