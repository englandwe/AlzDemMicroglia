# Detaching all packages
  lapply(names(sessionInfo()$otherPkgs), function(pkgs)
  detach(paste0('package:', pkgs),character.only = T,unload = T,force = T))

# Loading new packages
suppressMessages({
  library(Seurat, lib.loc = "~/R/x86_64-pc-linux-gnu-library/3.6/")
  library(sctransform)
  library(future)
  library(dplyr)
  library(cowplot)
  library(ggplot2)
  library(Rphenograph)
  library(extrafont)
})

#plan("multisession", workers = )
set.seed(030989)
cols <- c("black", "red", "darkgreen", "blue", "purple", "orange", "yellow", "lightblue", "darkred", "magenta", "green", "yellow4", "wheat4", "gray48", "plum1", "black")

# ================================================== Functions ===============================================
# Loading CellRanger data and filtering out previously defined Mouse cells, doublets, and gene poor cells
LoadV2 <- function(base.path=NULL, ID=NULL, meta = NULL, metaname = NULL){
  
  # Load the CellRanger human read files
  if(dir.exists(file.path(base.path, "Cell_Ranger", ID, "####","###"))){
    CR.path <- file.path(base.path, "Cell_Ranger", ID)
    CR.data <- Read10X(paste0(CR.path, "#####"))
    CR.genes <- as.matrix(rownames(CR.data))
    rownames(CR.data) <- sapply(strsplit(CR.genes, split='_'), '[', 2)
  } else {
    CR.path <- file.path(base.path, "Cell_Ranger", ID)
    CR.data <- Read10X(paste0(CR.path, "###"))    
  }
  
  # Generate a list of all barcodes from the Seurat object
  Hu.bars <- as.matrix(colnames(CR.data))
  
  # Load the mouse cell barcodes
  MS.path <- file.path(base.path, "Contamination", "V3", ID, "Ms_Barcodes.tsv")
  
  # Check to see if a mouse barcode file exists in the contamination folder. If not, copy it from the 
  # CellRanger output for the sample
  if(!(file.exists(MS.path))){
    file.copy(from=file.path(base.path, "Cell_Ranger", ID, "outs", "filtered_gene_bc_matrices", "###", "barcodes.tsv"),
              to=file.path(base.path, "Contamination", "V3", ID))
    file.rename(from=file.path(base.path, "Contamination", "V3", ID, "barcodes.tsv"),
                to=MS.path)
  }
  
  MS.bars <- as.matrix(read.table(MS.path, header=T))
  MS.bars <- as.matrix(sapply(strsplit(MS.bars, split='-'), '[[', 1L))
  
  # Remove contaminating barcodes from the list of barcodes in the Seurat object
  Hu.bars <- subset(Hu.bars, !(Hu.bars %in% MS.bars))
  
  # Filter the Seurat object to only contain the Human barcodes
  CR.data <- CR.data[, which(colnames(CR.data) %in% Hu.bars)]
  
  # Create the Seurat object using only the human barcodes
  scRNA <- CreateSeuratObject(CR.data, 
                              project = ID, 
                              min.cells = 3, 
                              min.features = 0)
  
  if (!(is.null(meta)) && !(is.null(metaname))) {
    
    for (i in 1:length(meta)) {
      
      scRNA[[metaname[[i]]]] <- meta[[i]]
      
    }
  }
  
  
  return(scRNA)
}

LoadV3 <- function(base.path=NULL, ID=NULL, meta = NULL, metaname = NULL){
  
  # Load the CellRanger human read files
  CR.path <- file.path(base.path, "Cell_Ranger", ID)
  CR.data <- Read10X(paste0(CR.path, "/outs/filtered_feature_bc_matrix/"))
  
  Hu <- grep("GRCh", rownames(CR.data), value=T)
  if (length(Hu) > 0){
    CR.data <- CR.data[which(rownames(CR.data) %in% Hu), ]
    CR.genes <- as.matrix(rownames(CR.data))
    rownames(CR.data) <- sapply(strsplit(CR.genes, split='_'), '[', 2)
  }
  
  scRNA <- CreateSeuratObject(CR.data, 
                              min.cells = 3,
                              min.features = 0,
                              project = ID)
  
  if (!(is.null(meta)) && !(is.null(metaname))) {
    
    for (i in 1:length(meta)) {
      
      
      scRNA[[metaname[[i]]]] <- meta[[i]]
      
    }
  }
  
  return(scRNA)
  
}

# ===================================== Import Data from Cell Ranger =========================================

# Indentifying the location and sample ID for the file being processed
base.path <- "~###"
Exp <- "###"
ID <- "###"
meta <- "###"
metaname <- "###"
day <- format(Sys.time(), "%m%d%y")
emit <- format(Sys.time(), "%H:%M:%S")

# CellRanger file version (either "V2" or "V3")
Vers <- "V3"

# Setting the output file directories
con.path <- file.path(base.path, Exp, "Contamination", "V3")
pre.dir <- file.path(base.path, Exp, "Seurat", "Preliminary", "V3")

# Creating the output directories for the contaminating cell information
if(!(file.exists(con.path))){
  dir.create(file.path(base.path, Exp, "Contamination"))
  dir.create(file.path(con.path))
}
if(!(file.exists(file.path(con.path, ID)))){
  dir.create(file.path(con.path, ID))
  dir.create(file.path(con.path, ID, "Barcodes"))
  dir.create(file.path(con.path, ID, "Plots"))
}

# Creating the output directories for the preliminary analysis files
if(!(file.exists(file.path(base.path, Exp, "Seurat")))){
  dir.create(file.path(base.path, Exp, "Seurat"))
  dir.create(file.path(base.path, Exp, "Seurat", "Preliminary"))
  dir.create(pre.dir)
} else if(!(file.exists(file.path(pre.dir)))){
  dir.create(file.path(pre.dir))  
}
if(!(file.exists(file.path(pre.dir, ID)))){
  dir.create(file.path(pre.dir, ID))
  dir.create(file.path(pre.dir, ID, "Plots"))
}

# Loading the CellRanger dataset
if(Vers == "V3"){
  scRNA <- LoadV3(base.path=file.path(base.path, Exp), ID=ID, meta = meta, metaname = metaname)
  
  totCells <- ncol(scRNA)
  print(paste0(totCells, " Total Cells"))
  
  avgUMI <- mean(scRNA$nCount_RNA)
  medUMI <- median(scRNA$nCount_RNA)
  avgGene <- mean(scRNA$nFeature_RNA)
  medGene <- median(scRNA$nFeature_RNA)
  print(paste0(avgUMI, " Average UMI per cell"))
  print(paste0(medUMI, " Median UMI per cell"))
  print(paste0(avgGene, " Average Gene per cell"))
  print(paste0(medGene, " Median Gene per cell"))
} else if(Vers == "V2"){
  scRNA <- LoadV2(base.path=file.path(base.path, Exp), ID=ID, meta = meta, metaname = metaname)
  
  totCells <- ncol(scRNA)
  print(paste0(totCells, " Total Cells"))
  
  avgUMI <- mean(scRNA$nCount_RNA)
  medUMI <- median(scRNA$nCount_RNA)
  avgGene <- mean(scRNA$nFeature_RNA)
  medGene <- median(scRNA$nFeature_RNA)
  print(paste0(avgUMI, " Average UMI per cell"))
  print(paste0(medUMI, " Median UMI per cell"))
  print(paste0(avgGene, " Average Gene per cell"))
  print(paste0(medGene, " Median Gene per cell"))
}


# ================================================================================================================
# =========================================== Begin QC Filtering =================================================
# ================================================================================================================

if(TRUE){
  # The number of genes and UMIs (nGene and nUMI) are automatically calculated
  # for every object by Seurat.  For non-UMI data, nUMI represents the sum of
  # the non-normalized values within a cell We calculate the percentage of
  # mitochondrial genes here and store it in percent.mito using AddMetaData.
  # We use object@raw.data since this represents non-transformed and
  # non-log-normalized counts The % of UMI mapping to MT-genes is a common
  # scRNA-seq QC metric.
  # AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
  scRNA[["percent.mito"]] <- PercentageFeatureSet(object = scRNA, pattern = "^MT-")
  scRNA[["percent.ribo"]] <- PercentageFeatureSet(object = scRNA, pattern = "^RP[SL]")


  #cowplot::plot_grid(VlnPlot(object = scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), ncol = 4),
    #                 VlnPlot(object = scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), ncol = 4, pt.size = 0),
     #                nrow=2, ncol=1)
  
  
  # GenePlot is typically used to visualize gene-gene relationships, but can
  # be used for anything calculated by the object, i.e. columns in
  # object@meta.data, PC scores etc.  Since there is a rare subset of cells
  # with an outlier level of high mitochondrial percentage and also low UMI
  # content, we filter these as well
  cowplot::plot_grid(FeatureScatter(object = scRNA, feature1 = "nCount_RNA", feature2 = "percent.mito"),
                     FeatureScatter(object = scRNA, feature1 = "nCount_RNA", feature2 = "percent.ribo"),
                     FeatureScatter(object = scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"), nrow = 1, ncol=3)
  
}


# ===============================================================================================================
# ======================================= Scaling and Normalization =============================================
# ===============================================================================================================

if(TRUE){
  # Normalizing the data using default settings
  scRNA <- NormalizeData(scRNA, normalization.method = 'LogNormalize', scale.factor = 10000)
  
  # Scaling and centering the data using the full transcriptome
  scRNA <- ScaleData(scRNA, features = rownames(scRNA), block.size = 10000)
  
  # Performing PCA analysis on the variable genes (default)
  # This will display the genes that are driving each of the PCs
  # With regressed UMI data, this can also be performed with the full transcriptome although it may be slow
  
  scRNA <- RunPCA(object = scRNA, features = rownames(scRNA), do.print = TRUE, ndims.print = 1:5, 
                  nfeatures.print = 10)
  
  # JackStraw method to quantitatively score each PC
  # This method is very time consuming so comment out unless needed
  #scRNA <- JackStraw(scRNA, num.replicate = 100)
  #scRNA <- ScoreJackStraw(scRNA, dims = 1:20)
  #JackStrawPlot(scRNA, dims=1:20)
  
  
  # A faster option that requires the user to identify where the "elbow" in the plot is
  ElbowPlot(object = scRNA, ndims = 40)
  
}


# ===============================================================================================================    
# ========================================== Clustering the cells ===============================================
# ===============================================================================================================

# The FindClusters function implements the procedure, and contains a resolution parameter that sets the 
# 'granularity' of the downstream clustering, with increased values leading to a greater number of clusters. 

PCs <- 11
Res <- 0.6

scRNA <- FindNeighbors(scRNA, dims = 1:PCs, force.recalc = T)
scRNA <- FindClusters(object = scRNA, resolution = Res)

# ===================================================================================================================   
# =========================================== Generating UMAP plots =================================================
# ===================================================================================================================

# Running UMAP with the previously calculated parameters
scRNA <- RunUMAP(scRNA, dims = 1:PCs)
DimPlot(scRNA, reduction = "umap", label = T, )

# Generating a dataframe with the analysis details
preCon <- data.frame(day, emit, "V3", Vers, "Pre_Contamination", ID, totCells, avgUMI, medUMI, 
                     avgGene, medGene, PCs, Res, "All_Genes")
colnames(preCon) <- c("Analysis_Date", "Time", "Seurat_Version", "CellRanger_Version",
                      "Stage", "Sample", "Total_Cells", "Average_UMI", "Median_UMI", 
                      "Average_Gene", "Median_Gene", "PCs", "Resolution", "PCA_Genes")


# ============================================================================================================
# ======================================= Saving QC Anlysis Data =============================================
# ============================================================================================================

# =========================================== Save UMAP plot ================================================= 

  tiff(filename= file.path(con.path, ID, "Plots", paste0("UMAP_", Exp, "_", ID, "_Clusters.tiff")), 
       width=5, height=5, units="in", res=300,
       compression="lzw")
  
  DimPlot(object = scRNA, pt.size = 1, label = T, reduction = 'umap') + theme_classic() +
    scale_color_manual(values = cols) +
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
  
# ===================================== Save TSNE Signatures plot ============================================   

  # Creating color palette for heatmap
  colfunc <- colorRampPalette(c("midnightblue", "lightskyblue1", "orange", "red"))
  
  # Create module scores for DAM genes
  if(TRUE){
    if (!(exists("PVM"))){
      Homeo <- data.frame(c("LYVE1", "ABCG2", "VSIG4", "ITM2C"))
      MHCI <- read.table("###", header=T, sep='\t')    
      MHCII <- read.table("###", header=T, sep='\t')    
      IFN <- data.frame(Gene = c("IFI6", "IRF7", "STAT1", "ISG15", "IFIT3"))
      DAM <- data.frame(Gene = c("CD9", "TREM2", "LPL", "ITGAX"))
      Degran <- data.frame(Gene = c("AGR2", "ANXA3", "PLA2G7", "PLAC8"))
      Apop <- data.frame(Gene = c("MDM2", "BAX", "ZMAT3", "GADD45A"))
      Canon <- data.frame(Gene = c("P2RY12", "P2RY13", "CX3CR1", "TMEM119"))    
      Infl_Chemo <- read.table("###", header = T)
      Hom_Chemo <- read.table("###", header = T)
      Mono_Chemo <- read.table("###", header = T)
      Tcell_Chemo <- read.table("###", header = T)
      Cycle <- read.table("###", header = F)
      PVM <- read.table("###", header = T)
      Foam <- read.table("###", header = T)
      Neur <- data.frame(c("###", "###", "###", "###"))
    }
    if(TRUE){
      scRNA <- AddModuleScore(scRNA, features = Homeo, name = "Homeo", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = MHCI, name = "MHCI", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = MHCII, name = "MHCII", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = IFN, name = "IFN", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = DAM, name = "DAM", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = Degran, name = "Degran", seed = 543)    
      scRNA <- AddModuleScore(scRNA, features = Apop, name = "Apoptotic", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = Canon, name = "Canonical", seed = 543)    
      scRNA <- AddModuleScore(scRNA, features = Infl_Chemo, name = "Infl_Chemo", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = Hom_Chemo, name = "Hom_Chemo", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = Mono_Chemo, name = "Mono_Chemo", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = Tcell_Chemo, name = "Tcell_Chemo", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = Cycle, name = "Cycle", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = PVM, name = "PVM", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = Foam, name = "Foam", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = Neur, name = "Neuron", seed = 543)
    }
  }
  
  tiff(filename= file.path(con.path, ID, "Plots", paste0("UMAP_", Exp, "_", ID, "_Signatures.tiff")), 
       width=20, height=20, units="in", res=100,
       compression="lzw")
  
  FeaturePlot(object = scRNA, pt.size = 0.5,
              features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mito",
                                "Homeo1", "MHCI1", "MHCII1", "IFN1", "DAM1",
                                "Degran1", "Apoptotic1", "Canonical1", "Hom_Chemo1",
                                "Infl_Chemo1", "Mono_Chemo1", "Tcell_Chemo1", 
                                "Cycle1", "PVM1", "Foam1", "Neuron1"), 
              cols = colfunc(20), 
              reduction = "umap",
              order = T)
  
  dev.off()  


# ================================= Save the barcodes for each cluster =======================================

  for (i in unique(eval(parse(text=paste0("scRNA$RNA_snn_res.", Res))))){
    
    assign(paste0("Cells", i), data.frame(Cells = WhichCells(scRNA, idents=i)))
    
    write.table(eval(parse(text=paste0("Cells", i))), 
                file.path(con.path, ID, "Barcodes", paste0(Exp, "_", ID, "_Clust", i, "_Barcodes.tsv")), 
                quote=F, sep='\t', row.names = F)
  }


# ======================================== Save QC Scatter plots =============================================
  
  tiff(filename = file.path(con.path, ID, "Plots", paste0("Scatter_", Exp, "_", ID, ".tiff")), 
       width=9, height=4.5, units="in", res=300,
       compression="lzw")
  
  cowplot::plot_grid(FeatureScatter(object = scRNA, 
                                    feature1 = "nCount_RNA", 
                                    feature2 = "percent.mito", 
                                    cols = cols) + NoLegend(),
                     FeatureScatter(object = scRNA, 
                                    feature1 = "nCount_RNA", 
                                    feature2 = "percent.ribo",
                                    cols = cols) + NoLegend(),
                     FeatureScatter(object = scRNA, 
                                    feature1 = "nCount_RNA", 
                                    feature2 = "nFeature_RNA",
                                    cols = cols) + NoLegend(), nrow = 1, ncol=3)
  
  dev.off()

# ============================ Find Significant Genes for all clusters ========================================
# find markers for every cluster compared to all remaining cells, report
# only the positive ones
scRNA.markers <- FindAllMarkers(object = scRNA, min.pct = 0.1, logfc.threshold = 0.25, 
                                only.pos = F)
scRNA.markers <- subset(scRNA.markers, scRNA.markers$p_val_adj <= 0.05)
scRNA.markers <- scRNA.markers[, c(7,1:6)]

scRNA.markers <- scRNA.markers[order(scRNA.markers$avg_logFC, decreasing=T),]
scRNA.markers <- scRNA.markers[order(scRNA.markers$cluster),]

write.table(scRNA.markers, 
            file.path(con.path, ID, paste0("SigGenes_", Exp, "_", ID, "_PreFilter_AllClusters.tsv")), 
            quote=F, sep='\t', row.names=F)
  
# ============================================================================================================
# ==================================== Removing Contaminating Cells ==========================================
# ================================ Identify and save contaminating barcodes ==================================
# ============================================================================================================

# Dividing cells cluster numbers
div <- c(6, 10)

# Write dividing cluster barcodes to a file
if (is.null(div)){
  dividing <- data.frame(Cells="NULL")
  write.table(dividing, file.path(con.path, ID, "Dividing_Barcodes.tsv"), sep='\t', quote=F, row.names = F)
}else{
  dividing <- data.frame(Cells=WhichCells(scRNA, idents = div))
  write.table(dividing, file.path(con.path, ID, "Dividing_Barcodes.tsv"), sep='\t', quote=F, row.names = F)
}

# Damaged cells cluster numbers
dam <- c(4, 5, 7, 8, 9)
damaged <- data.frame(Cells=WhichCells(scRNA, idents = dam))
write.table(damaged, file.path(con.path, ID, "GenePoor_Barcodes.tsv"), sep='\t', quote=F, row.names = F)

# Doublet cells cluster numbers
doub <- NULL

# Write doublet cluster barcodes to a file
if (is.null(doub)){
  doublet <- data.frame(Cells="NULL")
  write.table(doublet, file.path(con.path, ID, "Doublet_Barcodes.tsv"), sep='\t', quote=F, row.names = F)
}else{
  doublet <- data.frame(Cells=WhichCells(scRNA, idents = doub))
  write.table(doublet, file.path(con.path, ID, "Doublet_Barcodes.tsv"), sep='\t', quote=F, row.names = F)
}

# Mouse cell cluster numbers
#ms <- c(NULL)

# Write mouse cluster barcodes to a file
#if (is.null(ms)){
#  mouse <- data.frame(Cells="NULL")
#  write.table(mouse, file.path(con.path, ID, "Mouse_Barcodes.tsv"), sep='\t', quote=F, row.names = F)
#}else{
#  mouse <- data.frame(Cells=WhichCells(scRNA, idents = ms))
#  write.table(mouse, file.path(con.path, ID, "Mouse_Barcodes.tsv"), sep='\t', quote=F, row.names = F)
#}


# ========================== Remove Contaminating barcodes and recalc cell totals =============================    

scRNA <- subset(scRNA, cells = WhichCells(scRNA, idents = c(div, dam, doub), invert = T))

# Printing details about the Seurat object
totCells <- ncol(scRNA)
print(paste0(totCells, " Total Cells"))
avgUMI <- mean(scRNA$nCount_RNA)
medUMI <- median(scRNA$nCount_RNA)
avgGene <- mean(scRNA$nFeature_RNA)
medGene <- median(scRNA$nFeature_RNA)
print(paste0(avgUMI, " Average UMI per cell"))
print(paste0(medUMI, " Median UMI per cell"))
print(paste0(avgGene, " Average Gene per cell"))
print(paste0(medGene, " Median Gene per cell"))

# Writing the post contamination details to the analysis dataframe
postCon <- data.frame(day, emit, "V3", Vers, "Post_Contamination", ID, totCells, avgUMI, medUMI, avgGene, medGene)
colnames(postCon) <- c("Analysis_Date", "Time", "Seurat_Version", "CellRanger_Version", "Stage", "Sample", "Total_Cells", 
                       "Average_UMI", "Median_UMI", "Average_Gene", "Median_Gene")
# Merging the postCon datafrome with the preCon details
SampleInfo <- merge(preCon, postCon, all = T)


# ============================ Additional filtering of remaining cells =======================================  

cowplot::plot_grid(FeatureScatter(object = scRNA, feature1 = "nCount_RNA", feature2 = "percent.mito"),
                   FeatureScatter(object = scRNA, feature1 = "nCount_RNA", feature2 = "percent.ribo"),
                   FeatureScatter(object = scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"), nrow = 1, ncol=3)


# We filter out cells that have unique gene counts over 2,500 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate'.  -Inf and Inf should be used if you don't want a lower or upper
# threshold.

geneLow <- 500
geneHigh <- (2*medGene)

UMILow <- 500
UMIHigh <- (2*medUMI)

MitoLow <- -Inf
MitoHigh <- 15

RiboLow <- -Inf
RiboHigh <- 20


#test <- subset(scRNA, subset = nFeature_RNA > geneLow)
#test <- subset(test, subset = nCount_RNA > UMILow)
#test <- subset(test, subset = percent.mito > MitoLow & percent.mito < MitoHigh)
#test <- subset(test, subset = percent.ribo > RiboLow & percent.ribo < RiboHigh)


test <- subset(scRNA, subset = nFeature_RNA > geneLow & nFeature_RNA < geneHigh)
test <- subset(test, subset = nCount_RNA > UMILow & nCount_RNA < UMIHigh)
test <- subset(test, subset = percent.mito > MitoLow & percent.mito < MitoHigh)
test <- subset(test, subset = percent.ribo > RiboLow & percent.ribo < RiboHigh)


cowplot::plot_grid(FeatureScatter(object = test, feature1 = "nCount_RNA", feature2 = "percent.mito"),
                   FeatureScatter(object = test, feature1 = "nCount_RNA", feature2 = "percent.ribo"),
                   FeatureScatter(object = test, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"), nrow = 1, ncol=3)

# ========================= Check to see if the cutoffs look right before continuing ===========================  

scRNA <- test

totCells <- ncol(scRNA)
print(paste0(totCells, " Total Cells"))

avgUMI <- mean(scRNA$nCount_RNA)
medUMI <- median(scRNA$nCount_RNA)
avgGene <- mean(scRNA$nFeature_RNA)
medGene <- median(scRNA$nFeature_RNA)
print(paste0(avgUMI, " Average UMI per cell"))
print(paste0(medUMI, " Median UMI per cell"))
print(paste0(avgGene, " Average Gene per cell"))
print(paste0(medGene, " Median Gene per cell"))

# ===============================================================================================================
# ======================================= Scaling and Normalization =============================================
# ===============================================================================================================

if(TRUE){
  # Normalizing the data using default settings
  scRNA <- NormalizeData(scRNA, normalization.method = 'LogNormalize', scale.factor = 10000)
  
  # Scaling and centering the data using the full transcriptome
  scRNA <- ScaleData(scRNA, features = rownames(scRNA), model.use = 'negbinom', block.size = 10000)
  
  # Performing PCA analysis on the variable genes (default)
  # This will display the genes that are driving each of the PCs
  # With regressed UMI data, this can also be performed with the full transcriptome although it may be slow
  
  scRNA <- RunPCA(object = scRNA, features = rownames(scRNA), do.print = TRUE, ndims.print = 1:5, 
                  nfeatures.print = 10)
  
  # JackStraw method to quantitatively score each PC
  # This method is very time consuming so comment out unless needed
  #scRNA <- JackStraw(scRNA, num.replicate = 100)
  #scRNA <- ScoreJackStraw(scRNA, dims = 1:20)
  #JackStrawPlot(scRNA, dims=1:20)
  
  
  # A faster option that requires the user to identify where the "elbow" in the plot is
  ElbowPlot(object = scRNA, ndims = 40)
  
}

# ===============================================================================================================    
# ========================================== Clustering the cells ===============================================
# ===============================================================================================================

# The FindClusters function implements the procedure, and contains a resolution parameter that sets the 
# 'granularity' of the downstream clustering, with increased values leading to a greater number of clusters. 

PCs <- 11
Res <- 0.35

scRNA <- FindNeighbors(scRNA, dims = 1:PCs, force.recalc = T)
scRNA <- FindClusters(object = scRNA, resolution = Res)

# ===================================================================================================================   
# =========================================== Generating UMAP plots =================================================
# ===================================================================================================================

# Running UMAP with the previously calculated parameters
scRNA <- RunUMAP(scRNA, dims = 1:PCs)
DimPlot(scRNA, reduction = "umap", label = T)


if (TRUE) {
postFilt <- data.frame(day, emit, "V3", Vers, "Post_Filtering", ID, totCells, avgUMI, medUMI, avgGene, 
                       medGene, PCs, Res, "All_Genes", geneLow, geneHigh, UMILow, 
                       UMIHigh, MitoLow, MitoHigh, RiboLow, RiboHigh)
colnames(postFilt) <- c("Analysis_Date", "Time", "Seurat_Version", "CellRanger_Version",
                        "Stage", "Sample", "Total_Cells", "Average_UMI", "Median_UMI", 
                        "Average_Gene", "Median_Gene", "PCs", "Resolution", 
                        "PCA_Genes", "GeneLow", "GeneHigh", "UIMLow", "UMIHigh",
                        "MitoLow", "MitoHigh", "RiboLow", "RiboHigh")

SampleInfo <- merge(SampleInfo, postFilt, all = T)

if(file.exists(file.path(base.path, Exp, paste0(Exp, "_Prelim_DatasetInfo.tsv")))){
  
  write.table(SampleInfo, file=file.path(base.path, Exp, paste0(Exp, "_Prelim_DatasetInfo.tsv")),
              sep='\t', quote = F, row.names=F, append = T)
  
} else{
  
  write.table(SampleInfo, file=file.path(base.path, Exp, paste0(Exp, "_Prelim_DatasetInfo.tsv")),
              sep='\t', quote = F, row.names=F)
  
}  

# ============================ Find Significant Genes for all clusters ========================================
# find markers for every cluster compared to all remaining cells, report
# only the positive ones
scRNA.markers <- FindAllMarkers(object = scRNA, min.pct = 0.1, logfc.threshold = 0.25, 
                                only.pos = F)
scRNA.markers <- subset(scRNA.markers, scRNA.markers$p_val_adj <= 0.01)

write.table(scRNA.markers, 
            file.path(pre.dir, ID, paste0("SigGenes_", Exp, "_", ID, "_AllGroups.tsv")), 
            quote=F, sep='\t', row.names=F)

# Heatmap showing top genes per cluster
# DoHeatmap generates an expression heatmap for given cells and genes. In this case, we are plotting the 
# top 10 markers (or all markers if less than 10) for each cluster.
top10 <- scRNA.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
DoHeatmap(object = scRNA, features = top10$gene) + NoLegend()

ggsave(filename= file.path(pre.dir, ID, "Plots", paste0("Heatmap_", Exp, "_", ID, "_Top10.tiff")), 
       width=15, height=10, units="in", dpi=300, compression="lzw")

# ====================================== Save the Seurat Object ==============================================  

saveRDS(scRNA, file = file.path(pre.dir, ID, paste0(Exp, "_", ID, "_Seurat_Object.rds")))

# ============================================================================================================
# ==================================== Saving PostFilt Anlysis Data ==========================================
# ============================================================================================================

# =========================================== Save TSNE plot ================================================= 

DimPlot(object = scRNA, pt.size = 1, label = T, reduction = 'umap') + theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.line=element_line(size=1), 
        axis.ticks = element_line(color='black', size=3), 
        axis.ticks.length = unit(0, 'cm'), 
        axis.text = element_text(face='bold', color ='black', size=0),
        text = element_text(face='bold', color ='black', size=24),
        panel.border = element_blank(),
        plot.margin = unit(c(0,1,0,0), "cm"),
        legend.position='none') +
  xlab("UMAP 1") + ylab("UMAP 2")

ggsave(filename= file.path(pre.dir, ID, "Plots", paste0("UMAP_", Exp, "_", ID, "_Clusters.tiff")), 
       width=5, height=5, units="in", dpi=300, compression="lzw")

# ===================================== Save TSNE Signatures plot ============================================   

if(TRUE){  
  # Creating color palette for heatmap
  colfunc <- colorRampPalette(c("midnightblue", "lightskyblue1", "orange", "red"))
  
  # Create module scores for DAM genes
  if(TRUE){
    if (!(exists("PVM"))){
      Homeo <- data.frame(c("LYVE1", "ABCG2", "VSIG4", "ITM2C"))
      MHCI <- read.table("~/Sequencing/References/SingleCell/MHCI_GeneList.tsv", header=T, sep='\t')    
      MHCII <- read.table("~/Sequencing/References/SingleCell/MHCII_GeneList.tsv", header=T, sep='\t')    
      IFN <- data.frame(Gene = c("IFI6", "IRF7", "STAT1", "ISG15", "IFIT3"))
      DAM <- data.frame(Gene = c("CD9", "TREM2", "LPL", "ITGAX"))
      Degran <- data.frame(Gene = c("AGR2", "ANXA3", "PLA2G7", "PLAC8"))
      Apop <- data.frame(Gene = c("MDM2", "BAX", "ZMAT3", "GADD45A"))
      Canon <- data.frame(Gene = c("P2RY12", "P2RY13", "CX3CR1", "TMEM119"))    
      Infl_Chemo <- read.table("###", header = T)
      Hom_Chemo <- read.table("###", header = T)
      Mono_Chemo <- read.table("###", header = T)
      Tcell_Chemo <- read.table("###", header = T)
      Cycle <- read.table("###", header = F)
      PVM <- read.table("###", header = T)
      Foam <- read.table("###", header = T)
    }
    if(TRUE){
      scRNA <- AddModuleScore(scRNA, features = Homeo, name = "Homeo", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = MHCI, name = "MHCI", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = MHCII, name = "MHCII", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = IFN, name = "IFN", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = DAM, name = "DAM", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = Degran, name = "Degran", seed = 543)    
      scRNA <- AddModuleScore(scRNA, features = Apop, name = "Apoptotic", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = Canon, name = "Canonical", seed = 543)    
      scRNA <- AddModuleScore(scRNA, features = Infl_Chemo, name = "Infl_Chemo", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = Hom_Chemo, name = "Hom_Chemo", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = Mono_Chemo, name = "Mono_Chemo", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = Tcell_Chemo, name = "Tcell_Chemo", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = Cycle, name = "Cycle", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = PVM, name = "PVM", seed = 543)
      scRNA <- AddModuleScore(scRNA, features = Foam, name = "Foam", seed = 543)
    }
  }
}
  
  FeaturePlot(object = scRNA, pt.size = 0.5,
              features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mito",
                           "Homeo1", "MHCI1", "MHCII1", "IFN1", "DAM1",
                           "Degran1", "Apoptotic1", "Canonical1", "Hom_Chemo1",
                           "Infl_Chemo1", "Mono_Chemo1", "Tcell_Chemo1", 
                           "Cycle1", "PVM1", "Foam1", "Neuron1"), 
              cols = colfunc(20), 
              reduction = "umap")

  ggsave(filename= file.path(pre.dir, ID, "Plots", paste0("UMAP_", Exp, "_", ID, "_Signatures.tiff")), 
         width=20, height=20, units="in", dpi=100, compression="lzw")

# ====================================== Save PostFilt Scatter plots ===========================================

  cowplot::plot_grid(FeatureScatter(object = scRNA, 
                                    feature1 = "nCount_RNA", 
                                    feature2 = "percent.mito", 
                                    cols = cols) + NoLegend(),
                     FeatureScatter(object = scRNA, 
                                    feature1 = "nCount_RNA", 
                                    feature2 = "percent.ribo",
                                    cols = cols) + NoLegend(),
                     FeatureScatter(object = scRNA, 
                                    feature1 = "nCount_RNA", 
                                    feature2 = "nFeature_RNA",
                                    cols = cols) + NoLegend(), nrow = 1, ncol=3)
  
  ggsave(filename= file.path(pre.dir, ID, "Plots", paste0("Scatter_", Exp, "_", ID, ".tiff")), 
         width=9, height=4.5, units="in", dpi=300, compression="lzw")

}


# ============================================================================================================
# ==================================== Save the barcodes for CITE Seq ========================================
# ============================================================================================================
if (TRUE){

  # Creating the output directories for the contaminating cell information
  if(!(file.exists(file.path(base.path, Exp, "CITE_Output")))){
    dir.create(file.path(base.path, Exp, "CITE_Output"))
    dir.create(file.path(base.path, Exp, "CITE_Output", "HTO"))
    dir.create(file.path(base.path, Exp, "CITE_Output", "ADT"))
    dir.create(file.path(base.path, Exp, "CITE_Output", "Whitelist_Barcodes"))
    dir.create(file.path(base.path, Exp, "CITE_Output", "Whitelist_Barcodes", ID))
  } else if(!(file.exists(file.path(base.path, Exp, "CITE_Output", "Whitelist_Barcodes", ID)))){
    dir.create(file.path(base.path, Exp, "CITE_Output", "Whitelist_Barcodes", ID)) 
  }
  if(!(file.exists(file.path(base.path, Exp, "CITE_Output", "HTO", ID)))){
    dir.create(file.path(base.path, Exp, "CITE_Output", "HTO", ID))
  }
  if(!(file.exists(file.path(base.path, Exp, "CITE_Output", "ADT", ID)))){
    dir.create(file.path(base.path, Exp, "CITE_Output", "ADT", ID))
  }
  
  codes <- colnames(scRNA)
  write.table(codes,
              file.path(base.path, Exp, "CITE_Output", "Whitelist_Barcodes", ID, paste0(ID, "_HighQualityBarcodes_", Exp, ".tsv")), 
              quote=F, sep='\t', row.names = F, col.names = F)
}




