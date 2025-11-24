library(Seurat)
library(ggplot2)
library(viridis)

#Dotplot of top 5 upregulated DEGs per cluster - Figure 3b

DefaultAssay(allmerge) <- "RNA"
Idents(allmerge) <- allmerge$Clusts
allmarkers_origclust <- FindAllMarkers(object = allmerge, min.pct = 0.1, logfc.threshold = 0.1, only.pos = F)

top5 <- data.frame()

#gene names and pct expressed in cluster

for (i in levels(allmarkers_origclust$cluster)) {
    tmp5 <- allmarkers_origclust[allmarkers_origclust[6] == i ,]
    top5 <- rbind(top5,tmp5[c(1:5) ,])
}

top5_genelist <- top5$gene

DotPlot(object = allmerge, features = top5_genelist, scale.by="size") + 
scale_color_viridis() +
theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))

ggsave(file="clust_top5_dotplot_origclust.pdf", device="pdf")

#Downsampled UMAPs split by genotype - Figure 3d

table(allmerge$Genotype)
#smallest is 5797

Idents(allmerge) <- allmerge$Genotype

allmerge_geno_downsample <- subset(allmerge,downsample=5797)
allmerge_geno_downsample$Genotype <- factor(allmerge_geno_downsample$Genotype,levels=levels(allmerge_geno_downsample$Genotype)[c(4,1,2,3)])
Idents(allmerge_geno_downsample) <- allmerge_geno_downsample$Genotype

npscols <- c("#6E8C6B","#AE062D","#6DA0E7","#4E3BA3")
DimPlot(allmerge_geno_downsample,split.by="Genotype",cols=npscols,pt.size=1.5)
ggsave(filename="genotype_umap_downsampled_split_v2.pdf",device="pdf", dpi=100, width=25, height=7)


#Top 20 heatmaps - Figure 3g

for (i in levels(allmarkers_orig$cluster)) {
    tmp20 <- allmarkers_orig[allmarkers_orig$cluster == i, ]
    top20 <- tmp20[1:min(20, nrow(tmp20)), ]
    top20_avexp <- AverageExpression(object = allmerge, assays = "RNA",features = unique(top20$Gene),slot = "data",verbose = TRUE)
    top20_avexp_scaled <- t(scale(t(as.matrix(top20_avexp$RNA)), center = TRUE))
    outname <- paste0("heatmap_top20_", i, "_v4.pdf")
    outnamepng <- paste0("heatmap_top20_", i, "_v4.png")
    ht <- Heatmap(top20_avexp_scaled,row_names_side="left",cluster_columns=FALSE,col=rev(brewer.pal(name="RdYlBu",n=7)),clustering_method_rows = "complete",clustering_distance_rows = "euclidean",
        column_title = i,name=" ",rect_gp = gpar(col = "gray60", lwd = 1))
    pdf(outname, height=10,width=8)
    draw(ht)
    dev.off()
}

