#seurat object creation script to recreate figures from camp et al. 2015
# 26/04/2026 AKG last update

#used to start fresh
#rm(list = ls())
#graphics.off()

# Load libraries
#library(Seurat)
#library(ggplot2)
#library(dplyr)
#library(org.Hs.eg.db)
#library(AnnotationDbi)

# setting working directory(different for each user)
#setwd("~/steered_research_project")

# Load counts matrix (tabular in this case)
counts <- read.table(
  "counts_matrix_post2016.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)
# QC filtering 
counts <- as.matrix(counts)
mode(counts) <- "numeric"
counts <- round(counts)
#remove duplicates, and store matrix in counts
counts <- counts[!duplicated(rownames(counts)), ]
counts <- counts[rowSums(counts) > 0, ]

# Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "Post2016_UMAP",
  min.cells = 3,
  min.features = 200 #removes low-quality cells
)
#Quality control plots to show spread of data
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# quality control filtering
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)

# running PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
#elbow plot of PCA for standard devation
ElbowPlot(seurat_obj) #to determine optimal number of clusters (8 in this case)
DimPlot(seurat_obj, reduction = "pca") #calculate and display Principal components



# create clustering

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15) #find neighbour
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

#UMAP figure creation
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)

plot_umap <- DimPlot(
  seurat_obj,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.6
) + ggtitle("UMAP clustering")

plot_umap

# workflow of gene id to gene naming using biomrt package

markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) #mark top 15 markers in each cluster

print(top_markers)

#conversion gene id to symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(seurat_obj),
  column = "SYMBOL",
  keytype = "ENSEMBL", #ensembl database
  multiVals = "first"
)
#keep genes that successfully map from Ensembl,IDs to gene symbols
valid_genes <- !is.na(gene_symbols)
#create new value to use 
seurat_obj <- subset(seurat_obj, features = names(gene_symbols[valid_genes]))
gene_symbols <- gene_symbols[valid_genes]
gene_symbols <- make.unique(gene_symbols)

rownames(seurat_obj) <- gene_symbols

#update markers with new gene symbols
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) 

print(top_markers)

plotted_umap_symbols <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "seurat_clusters", #metadata clusters from seurat used
  label = TRUE,
  repel = TRUE,
  pt.size = 0.6
) + ggtitle("final UMAP")
plotted_umap_symbols

#naming of clusters 

seurat_obj <- RenameIdents(
  seurat_obj,
  "0" = "Early NPCs and Vental forebrain NPCs",
  "1" = "Cerebral cortex MN ",
  "2" = "Dorsal forebrain MCN",
  "3" = "Dorsal forebrain DN",
  "4" = "Dorsal Forerbain NPCs " ,
  "5" = "Mesenchymal cells Non-Cycling " ,
  "6" = "Mesenchymal cells Cycling" , 
  "7" = "Non-neural" 
)
#group new cluster names in cell type
seurat_obj$celltype <- Idents((seurat_obj))

#plot UMAP using new cluster names
plotted_umap_symbols <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "celltype", # Use the name of the column you just created
  label = TRUE,
  repel = TRUE,
  pt.size = 0.6
) + ggtitle("Post 2016 Figure 3D Recreation")

plotted_umap_symbols
#figure 3e creation based off top 15 significant features
FeaturePlot(
  seurat_obj,
  features = c("FOXG1", "ASPM", "LIN28A", "NEUROD6", "OTX2", "RSPO2", "MYT1L" , "DCN"),
  reduction = "umap",
  ncol = 4,
  order = TRUE
)
#violin plot displaying marker genes for clusters
#Figure 3F generation
VlnPlot(seurat_obj, features = c("FOXG1", "NEUROD6", "OTX2"), group.by = "celltype")
