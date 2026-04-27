### Recreation of Fig 1C-style clustering using Seurat, with marker-based cluster labelling and PC loading gene extraction for monocle ordering
### This script is for recreation of figure 2A from the paper using the STAR alternate pipeline (post-2016) analysis data
# INPUTS:  deseq2_counts_matrix.tsv
# OUTPUTS: outputs/STAR_fig1C_cell_clusters.csv
#          outputs/STAR_fig1C_cluster_markers_top40.csv
#          outputs/STAR_fig1C_seurat_object.rds
#          figures/STAR_Fig1C_style_PCA_paper_classes.png

## Loading required packages

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(org.Hs.eg.db)
})

if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("outputs")) dir.create("outputs")

## Loading the STAR dataset

message("Loading STAR alignment data...")
full_data <- read.table("deseq2_counts_matrix.tsv",
                        header           = TRUE,
                        sep              = "\t",
                        check.names      = FALSE,
                        row.names        = 1)

if ("gene_id" %in% colnames(full_data)) {
  full_data <- full_data[, colnames(full_data) != "gene_id"]
}
message("Raw: ", nrow(full_data), " genes x ", ncol(full_data), " cells")

## Converting the Ensembl IDs to gene symbols

message("Converting Ensembl IDs to gene symbols...")
gene_map <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys    = rownames(full_data),
                                  columns = c("ENSEMBL", "SYMBOL"),
                                  keytype = "ENSEMBL")
gene_map <- gene_map[!is.na(gene_map$SYMBOL) & gene_map$SYMBOL != "", ]
gene_map <- gene_map[!duplicated(gene_map$ENSEMBL), ]

dup_symbols     <- unique(gene_map$SYMBOL[duplicated(gene_map$SYMBOL)])
non_dup_ensembl <- gene_map$ENSEMBL[!gene_map$SYMBOL %in% dup_symbols]
ensembl_to_keep <- non_dup_ensembl

for (sym in dup_symbols) {
  ens_ids    <- gene_map$ENSEMBL[gene_map$SYMBOL == sym]
  ens_ids    <- ens_ids[ens_ids %in% rownames(full_data)]
  if (length(ens_ids) == 0) next
  total_expr <- rowSums(full_data[ens_ids, , drop = FALSE])
  ensembl_to_keep <- c(ensembl_to_keep, names(which.max(total_expr)))
}

full_data_mapped <- full_data[rownames(full_data) %in% ensembl_to_keep, ]
rownames(full_data_mapped) <- gene_map$SYMBOL[match(rownames(full_data_mapped),
                                                    gene_map$ENSEMBL)]
message("After symbol conversion: ",
        nrow(full_data_mapped), " genes x ", ncol(full_data_mapped), " cells")

## Filtering to identify fetal cells
-
message("Filtering to fetal cells...")
srr_nums  <- as.integer(sub("SRR", "", colnames(full_data_mapped)))
fetal_srr <- colnames(full_data_mapped)[srr_nums >= 2967101 &
                                          srr_nums <= 2967326]
message("Fetal cells: ", length(fetal_srr))
fetal_data <- full_data_mapped[, fetal_srr]

## Creating the Seurat object

message("Creating Seurat object...")
seu <- CreateSeuratObject(counts       = fetal_data,
                          project      = "STAR_fetal",
                          min.cells    = 3,
                          min.features = 200)
message("Cells after QC: ", ncol(seu))

## Preprocessing the data through normalisation, variable feature identification and data scaling

message("Preprocessing...")
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2500)
seu <- ScaleData(seu, features = VariableFeatures(seu))
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 20, verbose = FALSE)

## Clustering at resolution 2.0

message("Clustering...")
seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, resolution = 2.0)
message("Clusters: ", length(unique(Idents(seu))))
print(table(Idents(seu)))

## Attempting biologically driven cluster labelling

message("Assigning cluster labels...")

ap_markers  <- c("PAX6", "GLI3", "SOX2", "HES1", "VIM", "PROM1")
ap1_markers <- c("MKI67", "KIF11", "TOP2A", "PCNA")
bp_markers  <- c("EOMES", "INSM1", "HES6", "ASPM", "NEUROD4")
bp2_markers <- c("NEUROD6")
n_markers   <- c("MYT1L", "TBR1", "BCL11B", "NEUROD6")
n1_markers  <- c("NRP1", "ROBO2", "SEMA5A")
n3_markers  <- c("NRCAM", "FAT3", "CDH6")

all_score_genes     <- unique(c(ap_markers, ap1_markers, bp_markers,
                                bp2_markers, n_markers, n1_markers, n3_markers))
present_score_genes <- intersect(all_score_genes, rownames(seu))

avg_expr <- AverageExpression(seu, features = present_score_genes,
                              assays = "RNA", layer = "data")$RNA
colnames(avg_expr) <- sub("^g", "", colnames(avg_expr))

score_group <- function(genes, avg_mat) {
  g <- intersect(genes, rownames(avg_mat))
  if (length(g) == 0) return(rep(0, ncol(avg_mat)))
  colMeans(avg_mat[g, , drop = FALSE])
}

ap_score  <- score_group(ap_markers,  avg_expr)
ap1_score <- score_group(ap1_markers, avg_expr)
bp_score  <- score_group(bp_markers,  avg_expr)
bp2_score <- score_group(bp2_markers, avg_expr)
n_score   <- score_group(n_markers,   avg_expr)
n1_score  <- score_group(n1_markers,  avg_expr)
n3_score  <- score_group(n3_markers,  avg_expr)

score_df <- data.frame(
  cluster  = colnames(avg_expr),
  ap_score = ap_score, ap1_cc   = ap1_score,
  bp_score = bp_score, bp2_comm = bp2_score,
  n_score  = n_score,  n1_mig   = n1_score,  n3_mat = n3_score,
  stringsAsFactors = FALSE
)

score_df$dominant <- apply(score_df[, c("ap_score","bp_score","n_score")], 1,
                           function(x) c("AP","BP","N")[which.max(x)])

score_df$bp_ap_ratio <- score_df$bp_score / (score_df$ap_score + 0.01)
score_df$bp_n_ratio  <- score_df$bp_score / (score_df$n_score  + 0.01)
score_df$bp_index    <- score_df$bp_ap_ratio * score_df$bp_n_ratio

bp_threshold <- quantile(score_df$bp_index, 0.35)

score_df$dominant2 <- ifelse(
  score_df$n_score > score_df$ap_score * 1.2, "N",
  ifelse(score_df$bp_index >= bp_threshold &
           score_df$bp_score > score_df$ap_score * 0.7,
         "BP", score_df$dominant)
)

message("Dominant lineage per cluster:")
print(score_df[, c("cluster","dominant2","ap_score","bp_score","n_score","bp_index")])

assign_labels <- function(score_df) {
  labels      <- character(nrow(score_df))
  ap_clusters <- score_df[score_df$dominant2 == "AP", ]
  bp_clusters <- score_df[score_df$dominant2 == "BP", ]
  n_clusters  <- score_df[score_df$dominant2 == "N",  ]
  
  if (nrow(ap_clusters) >= 2) {
    ap_labs <- rep("AP2", nrow(ap_clusters))
    ap_labs[which.max(ap_clusters$ap1_cc)] <- "AP1"
    labels[match(ap_clusters$cluster, score_df$cluster)] <- ap_labs
  } else if (nrow(ap_clusters) == 1) {
    labels[match(ap_clusters$cluster, score_df$cluster)] <- "AP1"
  }
  
  if (nrow(bp_clusters) >= 2) {
    bp_labs <- rep("BP1", nrow(bp_clusters))
    bp_labs[which.max(bp_clusters$bp2_comm)] <- "BP2"
    labels[match(bp_clusters$cluster, score_df$cluster)] <- bp_labs
  } else if (nrow(bp_clusters) == 1) {
    labels[match(bp_clusters$cluster, score_df$cluster)] <- "BP1"
  }
  
  if (nrow(n_clusters) >= 3) {
    n1_row <- which.max(n_clusters$n1_mig)
    n3_row <- which.max(n_clusters$n3_mat)
    if (n1_row == n3_row) n3_row <- order(-n_clusters$n3_mat)[2]
    n_labs <- rep("N2", nrow(n_clusters))
    n_labs[n1_row] <- "N1"
    n_labs[n3_row] <- "N3"
    labels[match(n_clusters$cluster, score_df$cluster)] <- n_labs
  } else if (nrow(n_clusters) == 2) {
    n_labs <- rep("N2", 2)
    n_labs[which.max(n_clusters$n1_mig)] <- "N1"
    labels[match(n_clusters$cluster, score_df$cluster)] <- n_labs
  } else if (nrow(n_clusters) == 1) {
    labels[match(n_clusters$cluster, score_df$cluster)] <- "N1"
  }
  
  labels[labels == ""] <- "Other"
  labels
}

score_df$paper_label <- assign_labels(score_df)
message("Cluster to paper label:")
print(score_df[, c("cluster","dominant2","paper_label")])

score_df$cluster <- sub("^g", "", score_df$cluster)
cluster_to_label <- setNames(score_df$paper_label, score_df$cluster)
seu$paper_class  <- unname(cluster_to_label[as.character(trimws(Idents(seu)))])
seu$paper_class  <- factor(seu$paper_class,
                           levels = c("AP1","AP2","BP1","BP2","N1","N2","N3","Other"))

message("NA check: ", sum(is.na(seu$paper_class)))
message("Cells per paper class:")
print(table(seu$paper_class, useNA = "ifany"))

## Saving the required outputs 

cluster_map_STAR <- data.frame(
  cell_id       = colnames(seu),
  fig1c_cluster = as.character(Idents(seu)),
  paper_class   = as.character(seu$paper_class),
  stringsAsFactors = FALSE
)
write.csv(cluster_map_STAR, "outputs/STAR_fig1C_cell_clusters.csv",
          row.names = FALSE)
message("Saved outputs/STAR_fig1C_cell_clusters.csv")

message("Finding cluster markers...")
markers <- FindAllMarkers(seu, only.pos = TRUE,
                          min.pct = 0.20, logfc.threshold = 0.25)
write.csv(markers, "outputs/STAR_fig1C_cluster_markers_all.csv",
          row.names = FALSE)

top_markers <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 40) %>%
  ungroup()
write.csv(top_markers, "outputs/STAR_fig1C_cluster_markers_top40.csv",
          row.names = FALSE)
message("Saved outputs/STAR_fig1C_cluster_markers_top40.csv")

saveRDS(seu, "outputs/STAR_fig1C_seurat_object.rds")
message("Outputs: outputs/STAR_fig1C_cell_clusters.csv")
message("         outputs/STAR_fig1C_cluster_markers_top40.csv")

