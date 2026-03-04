library(Seurat)

#featureCounts output has metadata row at the top so skip
raw_counts <- read.table("/home/jvk3/SRE/camp_gene_counts.txt", header=TRUE, row.names=1, skip=1)

# Keep only 734 sample columns (removing Chr, Start, End, etc.)
counts <- raw_counts[, 6:ncol(raw_counts)]

# New names: "SRRXXXXXX_sorted.bam" -> "SRRXXXXXX"
colnames(counts) <- gsub(".*SRR", "SRR", colnames(counts))
colnames(counts) <- gsub("_.*", "", colnames(counts))

# Create Seurat Object
camp_seurat <- CreateSeuratObject(counts = counts, project = "Camp_Replicate", min.cells = 3, min.features = 200)

# Calculate Mitochondrial percentage
camp_seurat[["percent.mt"]] <- PercentageFeatureSet(camp_seurat, pattern = "^MT-")

# Normalize and find the 2000 most variable genes
camp_seurat <- NormalizeData(camp_seurat)
camp_seurat <- FindVariableFeatures(camp_seurat, selection.method = "vst", nfeatures = 2000)

#Scale the data and run PCA
camp_seurat <- ScaleData(camp_seurat)
camp_seurat <- RunPCA(camp_seurat, features = VariableFeatures(object = camp_seurat))

#Cluster the cells
camp_seurat <- FindNeighbors(camp_seurat, dims = 1:10)
camp_seurat <- FindClusters(camp_seurat, resolution = 0.5)

# 4. Run t-SNE (to match Figure 3D)
camp_seurat <- RunTSNE(camp_seurat, dims = 1:10)

# 5. Plot result
DimPlot(camp_seurat, reduction = "tsne", label = TRUE)

# Find markers for every cluster compared to all remaining cells
all_markers <- FindAllMarkers(camp_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Look at the top 2 genes for each cluster to match with the paper's cell types
library(dplyr)
all_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


# Define the new names based on our marker mapping
new_labels <- c("Dorsal NPCs", "Dorsal Progenitors", "CR Neurons", "Early NPCs", 
                "Ventral NPCs", "IFIT1+ Cells", "Mesenchymal", "XIST+ Cells")

# Apply labels to the Seurat object
names(new_labels) <- levels(camp_seurat)
camp_seurat <- RenameIdents(camp_seurat, new_labels)

# Generate final "Figure 3D" style plot
library(ggplot2)
DimPlot(camp_seurat, reduction = "tsne", label = TRUE, label.size = 5) + 
  ggtitle("Replication of Camp et al. Figure 3D") +
  NoLegend()