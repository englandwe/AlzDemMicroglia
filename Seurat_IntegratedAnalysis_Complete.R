# Detaching all packages
lapply(names(sessionInfo()$otherPkgs), function(pkgs)
  detach(paste0('package:', pkgs),character.only = T,unload = T,force = T))

suppressMessages({
  library(Seurat, lib.loc = "~/R/x86_64-pc-linux-gnu-library/3.6/")
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(sctransform)
  library(Rphenograph)
})

set.seed(030989)

# ================================================== Functions ========================================================

# Import preprocessed datasets for merging
readSample <- function(sampleID){
    
  # Read in the dataset
  obj <- readRDS(file = file.path(base.path, Exp, "Seurat", Analysis, sampleID, 
                                  paste0(Exp, "_", sampleID, "_Seurat_Object.rds")))
  
  # Add the original sample ID to the beginning of the cell barcodes
  obj <- RenameCells(obj, add.cell.id = gsub("_","", sampleID))
  Idents(obj) <- "orig.ident"
  
  return(obj)
  
}

# Import and merge between 2-4 datasets 
mergeDatasets <- function(sampleList){
  
  for (i in 1:length(sampleList)){
    
    print(paste0("Loading dataset ", sampleList[i], "..."))
    # Read in each dataset using the readSample function
    assign(paste0("scRNA_", as.character(sampleList[i])), readSample(sampleList[i]))
    
  }
  
  # Create a list of all the samples that were read in
  objList <- c(ls(pattern = "scRNA_*"))
  
  # Merging all of the loaded samples into a single sample
  for (i in seq(1, ((length(sampleList)*2)-2), by=2)){
    
    print(paste0("Merging ", objList[i], " and ", objList[i+1], "..."))
    
    assign(paste0("merge_", 
                  gsub("^.*_", "", objList[i]), 
                  gsub("^.*_", "", objList[i+1])), 
           merge(eval(parse(text = objList[[i]])),
                 eval(parse(text = objList[[i+1]]))))
    
    objList <- append(objList, paste0("merge_", 
                                      gsub("^.*_", "", objList[i]), 
                                      gsub("^.*_", "", objList[i+1])))
    
  }
  return(eval(parse(text=tail(objList, 1))))
}


# =========================== Import and merge up to four datasets from RDS objects ===================================
# Setting the object location variables
base.path <- "###"
Exp <- "###"
Addtl <- NULL
# Analysis where input samples are located
Analysis <- "Preliminary"
# Input samples IDs to be merged
samps <- c("###", "####")
ID <- "###"

# Import and merge the datasets listed above
scRNA <- mergeDatasets(samps)

# Determining the total cells in the dataset to confirm merging worked correctly
totCells <- ncol(scRNA)
print(paste0(totCells, " Total Cells"))

# Splittng the Seurat object for integrated analysis
scRNA.list <- SplitObject(scRNA, split.by = "orig.ident")


# ================================================================================================================
# =========================================== Begin Preprocessing ================================================
# ================================================================================================================

# Setting the output directory based on the input files
day <- format(Sys.time(), "%m%d%y")
outDir <- file.path(base.path, Exp, "Seurat/Integrated", paste0(day, "_Analysis"))

# Creating the output directories for the final analysis files
if(!(file.exists(file.path(base.path, Exp, "Seurat")))){
  dir.create(file.path(base.path, Exp, "Seurat"))
  dir.create(file.path(base.path, Exp, "Seurat/Integrated"))
  dir.create(outDir)
} else if(!(file.exists(file.path(base.path, Exp, "Seurat/Integrated")))){
  dir.create(file.path(base.path, Exp, "Seurat/Integrated"))
  dir.create(outDir)  
} else if(!(file.exists(outDir))){
  dir.create(outDir)  
}
if(!(file.exists(file.path(outDir, ID)))){
  dir.create(file.path(outDir, ID))
}


# ==========================================================================================================    
# ==================================== Seurat Integrated Analysis ==========================================
# ==========================================================================================================

# Normalize and Find Variable Genes for each dataset 
scRNA.list <- lapply(X = scRNA.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Removing mitochondiral and ribosomal genes from the variable gene list since these genes are more
# related to cellular damage than to interesting sources of variation
for (i in c(1:length(scRNA.list))){
  mito.genes <- grep(pattern = "^MT-", rownames(scRNA.list[[i]]), value = T)
  ribo.genes <- grep(pattern = "^RP[SL]", rownames(scRNA.list[[i]]), value = T)
  VariableFeatures(scRNA.list[[i]]) <- subset(VariableFeatures(scRNA.list[[i]]), !(VariableFeatures(scRNA.list[[i]]) %in% mito.genes))
  VariableFeatures(scRNA.list[[i]]) <- subset(VariableFeatures(scRNA.list[[i]]), !(VariableFeatures(scRNA.list[[i]]) %in% ribo.genes))
  scRNA.list[[i]][['RNA']]@meta.features$vst.variable <- rownames(HVFInfo(scRNA.list[[i]])) %in% VariableFeatures(scRNA.list[[i]])
  
  
}

# Determining genes to use as anchors across all datasets
# Using the first 25 PCs because this seems to work well and there aren't clear guidlines on
# validating this parameter
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNA.list, dims = 1:25)
scRNA.combined <- IntegrateData(anchorset = scRNA.anchors, dims = 1:25)

# Resetting the default assay to the integrated analysis
DefaultAssay(scRNA.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
# Scaling the integrated dataset and regressing out uninteresting sources of variation
regVars <- c("nCount_RNA", "percent.ribo", "percent.mito")
scRNA.combined <- ScaleData(scRNA.combined, vars.to.regress = regVars, verbose = FALSE)

# Performing PCA analysis
scRNA.combined <- RunPCA(scRNA.combined, npcs = 50, verbose = FALSE)

# Generating an Elbow plot to determine PC cutoffs for clustering
ElbowPlot(object = scRNA.combined, ndims = 50)

# Saving a copy of the Elbow plot
tiff(filename= file.path(outDir, ID, paste0("ElbowPlot_", Exp, "_", ID, ".tiff")), 
     width=15, height=10, units="in", res=300,
     compression="lzw")

ElbowPlot(object = scRNA.combined, ndims=50)

dev.off()

# =========================================================================================================
# =================================== Running Multicluster Analysis =======================================
# =========================================================================================================

# This portion of the pipeline is to determine the appropriate clustering parameters by clustering at
# multiple PC cutoffs and resolutions to identify the most stable clustering.

# Renaming the Seurat object because I'm too lazy to change the variable names in the code below
scRNA <- scRNA.combined

# Setting plot colors
Tcols <- c("#39b600", "#d89000", "#f8766d", "#00b0f6", "#a19f08", "orchid3", "#9590ff", "firebrick1", "gray40", 'purple')

# Creating the output directories for the final analysis files
if(!(file.exists(file.path(outDir, ID, "Iterations")))){
  dir.create(file.path(outDir, ID, "Iterations"))
  dir.create(file.path(outDir, ID, "Iterations/Genes"))
  dir.create(file.path(outDir, ID, "Iterations/Heatmaps"))
}


# Creating a list of the PCs to be used for Phenograph and MultiClustering
PCs <- list()
PCs[[1]] <- c(1:10) # Early cutoff from Elbow plot
PCs[[2]] <- c(1:18) # Late cutoff from Elbow plot
PCs[[3]] <- c(1:40) # Checking a distant cutoff to capture maximum complexity

names(PCs) <- c("Elbow1", "Elbow2", "Elbow3")

# =========================================================================================================
# ================================== Clustering at Multiple Resolutions ===================================
# =========================================================================================================

# Save a file containing the MultiCluster parameters
a <- list(c(PCs[[1]]), c(PCs[[2]]), c(PCs[[3]]))
dims <- data.frame("Cutoff" = c(names(PCs)), 
                   "PCs" = I(unlist(lapply(a,paste,collapse=","))),
                   "Date" = format(Sys.time(), "%m%d%y"),
                   "Time" = format(Sys.time(), "%H:%M:%S"))

write.table(dims, file=file.path(outDir, ID, "Iterations", paste0("IterationInfo_", ID, ".tsv")),
            sep = '\t', quote = F, row.names = F)


# Perform the MultiCluster analysis 
for (z in 1:length(PCs)){
  
  plot_list <- NULL
  remove(Res)
  
  # Find nearest neighbors
  scRNA <- FindNeighbors(scRNA, dims = PCs[[z]], force.recalc = F)
  
  for (Res in seq(0.1, 1.1, by=0.2)){
    # Find clusters and generate umap visualization
    scRNA <- FindClusters(object = scRNA, resolution = Res)
    scRNA <- RunUMAP(scRNA, dims = PCs[[z]])
    
    # Calculate and save the differentially expressed genes
    if (length(unique(Idents(scRNA))) > 1) {
      scRNA.markers <- FindAllMarkers(object = scRNA, min.pct = 0.1, logfc.threshold = 0.25, 
                                      only.pos = F)
      scRNA.markers <- subset(scRNA.markers, scRNA.markers$p_val_adj <= 0.01)
      scRNA.markers <- scRNA.markers[, c(7,1:6)]
      
      scRNA.markers <- scRNA.markers[order(scRNA.markers$avg_logFC, decreasing=T),]
      scRNA.markers <- scRNA.markers[order(scRNA.markers$cluster),]
      
      write.table(scRNA.markers, 
                  file= file.path(outDir, ID, "Iterations/Genes", paste0(names(PCs)[z], "_DiffGenes_PCs", max(PCs[[z]]), "_Res", Res , "_", ID, ".tsv")),
                  quote=F, sep='\t', row.names=F)
    } else {
      
      scRNA.markers <- data.frame(Genes = "NULL")
      
      write.table(scRNA.markers, 
                  file= file.path(outDir, ID, "Iterations/Genes", paste0(names(PCs)[z], "_NULLGenes_PCs", max(PCs[[z]]), "_Res", Res , "_", ID, ".tsv")),
                  quote=F, sep='\t', row.names=F)
      
    }
    
    
    # Creating heatmaps for each PC/Res pairing
    if (length(unique(Idents(scRNA))) > 1) {
      top10 <- scRNA.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
      
      heat <- DoHeatmap(object = scRNA, features = top10$gene) + NoLegend()
      
      ggsave(heat, file = file.path(outDir, ID, "Iterations/Heatmaps", paste0(names(PCs)[z], "_Heatmap_PCs", max(PCs[[z]]), "_Res", Res, "_", ID, "_Top10.tiff")),
             width = 15, height = 10, units = "in", dpi = 300)
      
    }
    
    plot_list[[ length(plot_list) + 1 ]] <- assign(paste0(names(PCs)[z], "_", Res), 
                                                   DimPlot(scRNA, reduction = "umap", pt.size = 0.8, label = T) + 
                                                     annotate('text', 
                                                              x = -Inf, 
                                                              y=Inf, 
                                                              hjust=0, 
                                                              vjust=1, 
                                                              label = paste0("Res=", Res, ", PCs=", max(PCs[[z]]))) + 
                                                     NoLegend())
    
  }
  
  a <- plot_grid(plotlist = plot_list, nrow = 2, ncol = 3)
  
  tiff(filename= file.path(outDir, ID, "Iterations", paste0(names(PCs)[z], "_UMAP_PCs", max(PCs[[z]]), "_Res0.1-1.1_", ID, "Combined.tiff")), 
       width=10, height=12, units="in", res=300, compression="lzw")
  
  print(a)
  
  dev.off()
  
}

# ===========================================================================================================    
# ======================================== Clustering the cells =============================================
# ===========================================================================================================

Clust_PCs <- PCs[[2]] # Cutoff that yielded the most appropriate clustering
Res <- 0.3 # Resolution the yielded the most appropriate clustering

scRNA <- FindNeighbors(scRNA, dims = Clust_PCs, force.recalc = T)
scRNA <- FindClusters(object = scRNA, resolution = Res)

# ===================================================================================================================   
# =========================================== Generating UMAP plots =================================================
# ===================================================================================================================

# Running UMAP with the previously calculated parameters
scRNA <- RunUMAP(scRNA, dims = Clust_PCs)
DimPlot(scRNA, reduction = "umap", label = T)

# Rename the clusters based on the defining markers for each cluster
# The heatmap from the multicluster analysis will provide the top markers
Idents(scRNA) <- scRNA$integrated_snn_res.0.3
new.cluster.ids <- c("Homeostatic", "MHCII", "IFN1", "DAM", "IFN2", "IL1B", "HumaninLike", "Secretory", "Macros")
names(new.cluster.ids) <- levels(scRNA)
scRNA$Clusts <- Idents(scRNA)
scRNA$Clusts <- new.cluster.ids[match(scRNA$Clusts, names(new.cluster.ids))]
scRNA$Clusts <- as.factor(scRNA$Clusts)
Idents(scRNA) <- scRNA$Clusts

DimPlot(scRNA, reduction = "umap", label = T)

# Save a heatmap showing top genes per cluster with correct cluster IDs
# DoHeatmap generates an expression heatmap for given cells and genes. In this case, we are plotting the 
# top 10 markers (or all markers if less than 10) for each cluster.
top10 <- scRNA.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
tiff(filename= file.path(outDir, ID, paste0("Heatmap_", Exp, "_", ID, ".tiff")), 
     width=15, height=10, units="in", res=300,
     compression="lzw")

DoHeatmap(object = scRNA, features = top10$gene, group.colors = Tcols) + NoLegend() 

dev.off()


# Save a UMAP plot with cluster ID labels
tiff(filename= file.path(outDir, ID, paste0("UMAP_", Exp, "_", ID, "_ClustersLabeled.tiff")), 
     width=5, height=5, units="in", res=300,
     compression="lzw")

DimPlot(object = scRNA, pt.size = 1, label = T, reduction = 'umap') + theme_classic() +
  scale_color_manual(values=Tcols) +
  theme(axis.line=element_line(size=1), 
        axis.ticks = element_line(color='black', size=3), 
        axis.ticks.length = unit(0, 'cm'), 
        axis.text = element_text(face='bold', color ='black', size=0),
        text = element_text(face='bold', color ='black', size=24),
        panel.border = element_blank(),
        plot.margin = unit(c(0,1,0,0), "cm"),
        legend.position='none') +
  xlab("UMAP 1") + ylab("UMAP 2")

dev.off()


# Save a UMAP the is split by the sample's orig.idents
if (length(unique(scRNA$orig.ident)) > 1){
  
  DimPlot(object = scRNA, split.by="orig.ident", pt.size = 1, label = T, reduction = 'umap') + 
    theme_classic() +
    scale_color_manual(values = Tcols) +
    theme(axis.line=element_line(size=1), 
          axis.ticks = element_line(color='black', size=3), 
          axis.ticks.length = unit(0, 'cm'), 
          axis.text = element_text(face='bold', color ='black', size=0),
          text = element_text(face='bold', color ='black', size=24),
          panel.border = element_blank(),
          plot.margin = unit(c(0,1,0,0), "cm"),
          legend.position='none') +
    xlab("UMAP 1") + ylab("UMAP 2")
  
  ggsave(filename = file.path(outDir, ID, paste0("UMAP_", Exp, "_", ID, "_SplitOrigIdent.tiff")),
         width=15, height=8, units="in", dpi=300, device = "tiff",
         compression = "lzw")
  
}


# Saving the Seurat object
saveRDS(scRNA, file = file.path(outDir, ID, paste0(Exp, "_", ID, "_Seurat_Object.rds")))




