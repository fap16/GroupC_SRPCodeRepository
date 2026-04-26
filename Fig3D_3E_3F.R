# Load required libraries
library(Seurat)      # single-cell analysis framework
library(FactoMineR)  # advanced PCA for feature extraction
library(tidyverse)   # data manipulation and visualisation

#1. SETUP & METADATA
raw_counts <- read.csv("/home/jvk3/SRE/camp_organoid_matrix_V2.csv", row.names = 1)
meta_data <- read.csv("/home/jvk3/SRE/SraRunTable.csv", row.names = 1)

# initialise Seurat object
# min.cells = 3: Filters out low-abundance genes (technical noise)
# min.features = 200: Removes low-quality cells (viability filter)
organoid <- CreateSeuratObject(counts = raw_counts, project = "Camp_Replication", 
                               min.cells = 3, min.features = 200)

head(rownames(organoid))

# add sample-specific metadata (SRR IDs, stage, etc.)
organoid <- AddMetaData(object = organoid, metadata = meta_data)

head(organoid@meta.data)

# 2. STANDARD SEURAT PROCESSING
# adjusts for sequencing depth; makes gene expression comparable across all cells
organoid <- NormalizeData(organoid)

# identifies the top 2000 genes with high variation to focus on biological signal
organoid <- FindVariableFeatures(organoid, selection.method = "vst", nfeatures = 2000)

# centers and scales data (Z-score); ensures all genes have equal weight for PCA
organoid <- ScaleData(organoid)

# 3. FACTOMINER PCA 
# authors used FactoMineR to identify specific genes significantly 
# correlated with Principal Components rather than standard Seurat PCA.
scaled_matrix <- t(as.matrix(GetAssayData(organoid, layer = "scale.data")[VariableFeatures(organoid), ]))
res.pca <- PCA(scaled_matrix, ncp = 15, graph = FALSE)

# 4. DIMDESC FEATURE EXTRACTION
# extract genes significantly correlated (p < 0.001) with the first 15 PCs.
# matches the "top 50 genes per PC" constraint mentioned.
pca_stats <- dimdesc(res.pca, proba = 0.001) 
paper_features <- c()
for (i in 1:15) {
  pc_name <- paste0("Dim.", i)
  if (!is.null(pca_stats[[pc_name]]$quanti)) {
    genes_pc <- rownames(pca_stats[[pc_name]]$quanti)[1:min(50, nrow(pca_stats[[pc_name]]$quanti))]
    paper_features <- c(paper_features, genes_pc)
  }
}
paper_features <- unique(paper_features)

# 5. INTEGRATE PCA & RUN t-SNE
# perplexity = 5: matches the authors choice for low-density mapping.
# embeddings are manually injected from FactoMineR into the Seurat object.
pca_embeddings <- res.pca$ind$coord
colnames(pca_embeddings) <- paste0("PC_", 1:15)
organoid[["pca"]] <- CreateDimReducObject(embeddings = pca_embeddings, key = "PC_", assay = "RNA")

organoid <- RunTSNE(organoid, features = paper_features, dims = 1:15, 
                    perplexity = 5, check_duplicates = FALSE)

# 6. CLUSTERING
# builds Shared Nearest Neighbor (SNN) graph using Euclidean distance in PCA space
# defines the neighborhood of each cell based on its 15 most significant PCs
organoid <- FindNeighbors(organoid, dims = 1:15)

# applies the Louvain algorithm to partition the graph into clusters
# resolution 1.0 is a standard setting that captures the diverse cell types in the organoid
organoid <- FindClusters(organoid, resolution = 1.0)

# organises clusters into a hierarchical tree based on average gene expression 
# reveals developmental relationships between different cell groups
organoid <- BuildClusterTree(organoid, dims = 1:15, reorder = TRUE)

# 7. MARKER ANALYSIS, ROC TEST
# test.use = "roc" calculates the Area Under the Curve (AUC) for each gene,
# identifying markers based on predictive power rather than just fold-change.
Idents(organoid) <- organoid$seurat_clusters # ensure we are testing clusters, not samples
markers_roc <- FindAllMarkers(organoid, test.use = "roc", only.pos = TRUE, min.pct = 0.25)

# 8. BIOLOGICAL MAPPING
# assigns biological identities to clusters based on differential expression (ROC test) 
# and marker gene signatures (e.g. FOXG1, OTX2, GAD1) described in paper.
organoid <- RenameIdents(organoid, "0" = "Midbrain Progenitors","1" = "Midbrain-Hindbrain", 
"2" = "Forebrain NPCs", "3" = "Dorsal Forebrain Neurons", "4" = "Midbrain Neurons",          
"5" = "Ventral Forebrain NPCs", "6" = "Early Forebrain Progs", "7" = "Mesenchymal Cells",         
"8" = "Forebrain-Specific NPCs", "9" = "Late Midbrain Progs")

#9. FIGURE GENERATION 

# GENE MAP 
gene_map <- c("ENSG00000176165" = "FOXG1", "ENSG00000164600" = "NEUROD6",
"ENSG00000165588" = "OTX2", "ENSG00000147655" = "RSPO2", "ENSG00000011465" = "DCN",
"ENSG00000066279" = "ASPM", "ENSG00000131914" = "LIN28A", "ENSG00000186487" = "MYT1L")

# FIGURE 3D (t-SNE)
p3d <- DimPlot(organoid, reduction = "tsne", label = TRUE, repel = TRUE, pt.size = 1.2) +
  ggtitle("Figure 3D Replication") + theme_minimal()

print(p3d)

# FIGURE 3E (FeaturePlots Grid)
p3e <- FeaturePlot(organoid, features = names(gene_map), ncol = 4, 
                   cols = c("lightgrey", "firebrick"), pt.size = 0.5)
for(i in 1:length(gene_map)) { p3e[[i]] <- p3e[[i]] + ggtitle(gene_map[i]) }

print(p3e)

# FIGURE 3F (Violin plots)
# map specific stages to the 'r1-r4' and 'fetal' groups used in Figure 3F
organoid$Paper_Identity <- NA 

# split Day 53 cells into two groups to recreate 'r1' and 'r2' replicates
d53 <- which(organoid$stage == "53 days"); m53 <- floor(length(d53)/2)
organoid$Paper_Identity[d53[1:m53]] <- "r1"
organoid$Paper_Identity[d53[(m53+1):length(d53)]] <- "r2"

# split Day 58 cells ... the 'r3' and 'r4' replicates
d58 <- which(organoid$stage == "58 days"); m58 <- floor(length(d58)/2)
organoid$Paper_Identity[d58[1:m58]] <- "r3"
organoid$Paper_Identity[d58[(m58+1):length(d58)]] <- "r4"

# identify fetal samples by searching for "weeks" in the stage metadata and labeling as 'fetal'
organoid$Paper_Identity[grep("weeks", organoid$stage)] <- "fetal"

# subset and order data for the timeline (r1 -> r4 -> fetal)
organoid_subset <- subset(organoid, subset = Paper_Identity %in% c("r1", "r2", "r3", "r4", "fetal"))
organoid_subset$Paper_Identity <- factor(organoid_subset$Paper_Identity, 
                                         levels = c("r1", "r2", "r3", "r4", "fetal"))

genes_3f <- c("ENSG00000176165", "ENSG00000164600", "ENSG00000165588")
gene_titles <- c("FOXG1", "NEUROD6", "OTX2")

p3f <- VlnPlot(organoid_subset, features = genes_3f, group.by = "Paper_Identity", 
               ncol = 3, pt.size = 0.1, cols = rep("lightgrey", 5)) 

for(i in 1:3) {
  p3f[[i]] <- p3f[[i]] + ggtitle(gene_titles[i]) + theme_bw() + 
    theme(legend.position = "none", axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 10))
}

print(p3f)


