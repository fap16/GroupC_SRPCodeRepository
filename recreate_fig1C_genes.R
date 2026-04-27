### Recreation of Fig 1C-style clustering using Seurat, with marker-based cluster labelling and PC loading gene extraction for monocle ordering
### This script is for recreation of figure 2A from the paper using the original dataset


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

infile <- "GSE75140_hOrg.fetal.master.data.frame.txt.gz"

if (!file.exists(infile)) {
  stop("Missing file: ", normalizePath(infile, mustWork = FALSE))
}

if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("outputs")) dir.create("outputs")

message("Reading GEO matrix...")
raw_df <- read.table(
  gzfile(infile),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

colnames(raw_df)[1] <- "cell_id"

gene_cols <- colnames(raw_df)[-1]
gene_cols <- sub("^X\\.", "", gene_cols)
gene_cols <- sub("\\.+$", "", gene_cols)
gene_cols <- gsub('^"|"$', "", gene_cols)

cell_ids <- raw_df$cell_id

expr_matrix <- data.matrix(raw_df[, -1])
rownames(expr_matrix) <- cell_ids
colnames(expr_matrix) <- gene_cols

if (anyNA(expr_matrix)) {
  warning("NAs detected; replacing with 0.")
  expr_matrix[is.na(expr_matrix)] <- 0
}

# Seurat expects genes x cells, therefore ensuring rows and columns are correctly assigned
expr_matrix <- t(expr_matrix)

message("Creating Seurat object...")
seu <- CreateSeuratObject(
  counts = expr_matrix,
  project = "GSE75140_fetal",
  min.cells = 3,
  min.features = 200
)

rm(expr_matrix, raw_df)

DefaultAssay(seu) <- "RNA"

message("Preprocessing...")
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2500)
seu <- ScaleData(seu, features = VariableFeatures(seu))
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 20, verbose = FALSE)

message("Clustering...")
seu <- FindNeighbors(seu, dims = 1:10)

seu <- FindClusters(seu, resolution = 0.50)

png("figures/Fig1C_style_PCA_clusters.png", width = 7, height = 5, units = "in", res = 300)
print(
  DimPlot(seu, reduction = "pca", group.by = "seurat_clusters") +
    ggtitle("Fig 1C-style PCA clustering")
)
dev.off()

message("Cluster sizes:")
print(table(Idents(seu)))

message("Finding markers...")
markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.20,
  logfc.threshold = 0.25
)

write.csv(markers, "outputs/fig1C_cluster_markers_all.csv", row.names = FALSE)

# Using top 40 markers per cluster to get a stronger classification gene set
top_markers <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 40) %>%
  ungroup()

write.csv(top_markers, "outputs/fig1C_cluster_markers_top40.csv", row.names = FALSE)

classification_genes <- sort(unique(top_markers$gene))
writeLines(classification_genes, "outputs/fig1C_classification_genes.txt")
saveRDS(classification_genes, "outputs/fig1C_classification_genes.rds")

canonical_markers <- c(
  "PAX6", "GLI3", "SOX2", "HES1", "VIM", "PROM1",
  "EOMES", "INSM1", "HES6", "ASPM", "NEUROD4",
  "MYT1L", "TBR1", "BCL11B", "NEUROD6"
)

present_markers <- canonical_markers[canonical_markers %in% rownames(seu)]

# Scaling again including these markers so heatmap does not omit them
seu <- ScaleData(seu, features = unique(c(VariableFeatures(seu), present_markers)))

avg_expr <- AverageExpression(seu, features = present_markers, assays = "RNA", layer = "data")$RNA
write.csv(avg_expr, "outputs/fig1C_cluster_average_expression.csv")

saveRDS(seu, "outputs/fig1C_seurat_object.rds")

message(" - outputs/fig1C_cluster_markers_all.csv")
message(" - outputs/fig1C_cluster_markers_top40.csv")
message(" - outputs/fig1C_classification_genes.rds")
message(" - outputs/fig1C_cluster_average_expression.csv")
message(" - outputs/fig1C_seurat_object.rds")
cluster_map <- data.frame(
  cell_id = colnames(seu),
  fig1c_cluster = as.character(Idents(seu)),
  stringsAsFactors = FALSE
)

write.csv(cluster_map, "outputs/fig1C_cell_clusters.csv", row.names = FALSE)

