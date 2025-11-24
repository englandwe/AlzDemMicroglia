library(Vennerable)

markers_5x_vs_wt <- FindMarkers(allmerge,ident.1="5X",ident.2="WT",group.by="Genotype",logfc.threshold=0.1)
markers_5x_vs_wt_filtered <- subset(markers_5x_vs_wt,p_val_adj < 0.01)

markers_PS19_vs_wt <- FindMarkers(allmerge,ident.1="PS19",ident.2="WT",group.by="Genotype",logfc.threshold=0.1)
markers_PS19_vs_wt_filtered <- subset(markers_PS19_vs_wt,p_val_adj < 0.01)

markers_PS5x_vs_wt <- FindMarkers(allmerge,ident.1="PS5X",ident.2="WT",group.by="Genotype",logfc.threshold=0.1)
markers_PS5x_vs_wt_filtered <- subset(markers_PS5x_vs_wt,p_val_adj < 0.01)


markers_5X_filtered_up <- subset(markers_5x_vs_wt_filtered,X3 > 0)$X1
markers_5X_filtered_down <- subset(markers_5x_vs_wt_filtered,X3 < 0)$X1

markers_PS19_filtered_up <- subset(markers_PS19_vs_wt_filtered,X3 > 0)$X1
markers_PS19_filtered_down <- subset(markers_PS19_vs_wt_filtered,X3 < 0)$X1

markers_PS5X_filtered_up <- subset(markers_PS5x_vs_wt_filtered,X3 > 0)$X1
markers_PS5X_filtered_down <- subset(markers_PS5x_vs_wt_filtered,X3 < 0)$X1

geno_vennlist_up <- list(markers_5X_filtered_up,markers_PS19_filtered_up,markers_PS5X_filtered_up)
geno_vennlist_down <- list(markers_5X_filtered_down,markers_PS19_filtered_down,markers_PS5X_filtered_down)

names(geno_vennlist_up) <- c("5X","PS19","PS5X")
names(geno_vennlist_down) <- c("5X","PS19","PS5X")

geno_up_venn <- Venn(geno_vennlist_up)
plot(geno_up_venn,doWeights=FALSE)

dev.copy(pdf,'geno_venns/heatmap_cluster_markers_by_genotype_up.pdf',width=7,height=7)
dev.off()

geno_down_venn <- Venn(geno_vennlist_down)
plot(geno_down_venn,doWeights=FALSE)

dev.copy(pdf,'geno_venns/heatmap_cluster_markers_by_genotype_down.pdf',width=7,height=7)
dev.off()

for (i in seq(1,length(geno_down_venn@IntersectionSets))) {
  setnm <- names(geno_down_venn@IntersectionSets[i])
  for (j in geno_down_venn@IntersectionSets[i]) {
    outln <- paste(setnm,j,sep="\t" )
    write_lines(outln, append = TRUE,file = "geno_venns/down_genes_by_set.txt")
  }
}

for (i in seq(1,length(geno_up_venn@IntersectionSets))) {
  setnm <- names(geno_up_venn@IntersectionSets[i])
  for (j in geno_up_venn@IntersectionSets[i]) {
    outln <- paste(setnm,j,sep="\t" )
    write_lines(outln, append = TRUE,file = "geno_venns/up_genes_by_set.txt")
  }
}
