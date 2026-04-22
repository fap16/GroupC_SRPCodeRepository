#seurat object creation script to recreate figures from camp et al. 2015
# 21/04/2026 AKG last update

#used to start fresh
#rm(list = ls())
#graphics.off()

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

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
  min.features = 0
)
#Quality control plots to show spread of data
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# quality control filterinng
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 500)
seurat_obj <- ScaleData(seurat_obj)

# running PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
#elbow plot of pca for standard devation
ElbowPlot(seurat_obj)
DimPlot(seurat_obj, reduction = "pca")

# create clustering

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

#UMAP figure creation
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

plot_umap <- DimPlot(
  seurat_obj,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.6
) + ggtitle("UMAP clustering")

plot_umap

# workflow ofr gene id to gene naming using bioconductor

markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

print(top_markers)

#conversion gene id to symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(seurat_obj),
  column = "SYMBOL",
  keytype = "ENSEMBL", #ensembl database
  multiVals = "first"
)
#input new gene symbols in valid genes with top 10 markers
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
  slice_max(avg_log2FC, n = 10)

print(top_markers)

plotted_umap_symbols <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "seurat_clusters", #metadata clasters from seurat used
  label = TRUE,
  repel = TRUE,
  pt.size = 0.6
) + ggtitle("final UMAP")
plotted_umap_symbols

#naming of clusters 

seurat_obj <- RenameIdents(
  seurat_obj,
  "0" = "Interneurons",
  "1" = "Cortical neurons",
  "2" = "Radial Glia cells",
  "3" = "Cycling apical Progenitors",
  "4" = "Newborn neurons " ,
  "5" = "Intermediate neurons " ,
  "6" = "Ventral forebrain Neurons" , 
  "7" = "Mesenchymal cells" 
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
) + ggtitle("Post 2016 figure 3d Recreation")

plotted_umap_symbols
#figure 3e creation based off significant features
FeaturePlot(
  seurat_obj,
  features = c("FOXG1", "NFIA", "NFIB", "NEUROD6", "OTX2", "RSPO2", "WNT2B"),
  reduction = "umap",
  ncol = 3,
  order = TRUE
)
#violin plot displaying marker genes for clusters
VlnPlot(seurat_obj, features = c("FOXG1", "NEUROD6", "OTX2"), group.by = "celltype")
