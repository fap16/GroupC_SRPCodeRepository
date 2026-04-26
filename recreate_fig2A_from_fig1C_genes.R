### Recreate Fig 2A using recreated Fig 1C classification genes AND 7 Fig 1C clusters

## Importing required packages
suppressPackageStartupMessages({
  library(monocle)
  library(Matrix)
  library(igraph)
  library(ggplot2)
})

## Assigning required input and output files
infile <- "GSE75140_hOrg.fetal.master.data.frame.txt.gz"
gene_file <- "outputs/fig1C_classification_genes.rds"
cluster_file <- "outputs/fig1C_cell_clusters.csv"

if (!file.exists(infile)) stop("Missing file: ", normalizePath(infile, mustWork = FALSE))
if (!file.exists(gene_file)) stop("Missing file: ", normalizePath(gene_file, mustWork = FALSE))
if (!file.exists(cluster_file)) stop("Missing file: ", normalizePath(cluster_file, mustWork = FALSE))

if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("outputs")) dir.create("outputs")

# Loading in the original data GEO file matrix
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

expr_matrix <- t(expr_matrix)  # genes x cells - assigning data to rows and columns

if (anyNA(expr_matrix)) {
  warning("NAs detected in matrix; replacing with 0.")
  expr_matrix[is.na(expr_matrix)] <- 0
}

expr_matrix <- as(expr_matrix, "dgCMatrix")
rm(raw_df)

## Building the Monocle dataset for trajectory analysis and figure recreation of the original Figure 2A. 
message("Creating Monocle CellDataSet...")
pd_df <- data.frame(
  cell_id = colnames(expr_matrix),
  row.names = colnames(expr_matrix),
  stringsAsFactors = FALSE
)

fd_df <- data.frame(
  gene_short_name = rownames(expr_matrix),
  row.names = rownames(expr_matrix),
  stringsAsFactors = FALSE
)

pd <- new("AnnotatedDataFrame", data = pd_df)
fd <- new("AnnotatedDataFrame", data = fd_df)

cds <- newCellDataSet(
  expr_matrix,
  phenoData = pd,
  featureData = fd,
  expressionFamily = negbinomial.size()
)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

available_genes <- rownames(exprs(cds))

## Generating a robust gene resolver
resolve_one <- function(g, available) {
  g_clean <- toupper(trimws(gsub('^"|"$', "", gsub("\\.+$", "", sub("^X\\.", "", g)))))
  avail_clean <- toupper(trimws(gsub('^"|"$', "", gsub("\\.+$", "", sub("^X\\.", "", available)))))
  idx <- match(g_clean, avail_clean)
  if (!is.na(idx)) return(available[idx])
  NA_character_
}

resolve_genes <- function(gene_vec, available) {
  out <- character(0)
  for (g in gene_vec) {
    hit <- resolve_one(g, available)
    if (!is.na(hit)) out <- c(out, hit)
  }
  unique(out)
}

## loading recreated Fig 1C classification genes 
classification_genes_raw <- readRDS(gene_file)
classification_genes <- resolve_genes(classification_genes_raw, available_genes)
classification_genes <- classification_genes[!is.na(classification_genes)]

message("Classification genes loaded: ", length(classification_genes))

if (length(classification_genes) < 20) {
  stop("Too few classification genes survived matching.")
}

## Loading Fig 1C cluster assignments based on classification genes

cluster_map <- read.csv(cluster_file, stringsAsFactors = FALSE)
rownames(cluster_map) <- cluster_map$cell_id

pData(cds)$fig1c_cluster <- cluster_map[colnames(exprs(cds)), "fig1c_cluster"]

cluster_levels <- sort(unique(na.omit(pData(cds)$fig1c_cluster)))
pData(cds)$fig1c_cluster <- factor(pData(cds)$fig1c_cluster, levels = cluster_levels)

message("Fig 1C clusters loaded:")
print(table(pData(cds)$fig1c_cluster, useNA = "ifany"))

## Filtering to identify fetal cells and attaching labels

# Striping quotes from cluster_map
cluster_map$cell_id <- gsub('^"|"$', "", cluster_map$cell_id)

# Striping quotes from cds cell IDs via the assay slot
clean_ids <- gsub('^"|"$', "", colnames(exprs(cds)))
cds@assayData$exprs@Dimnames[[2]] <- clean_ids
rownames(pData(cds)) <- clean_ids
pData(cds)$cell_id  <- clean_ids

# Syncing protocolData with cleaned cell IDs
pd_clean <- pData(cds)
proto    <- new("AnnotatedDataFrame",
                data = data.frame(row.names = rownames(pd_clean),
                                  labelDescription = rep(NA, nrow(pd_clean))))
protocolData(cds) <- proto

# Filtering to fetal cells using cluster_map cell IDs directly
keep_in_matrix <- intersect(cluster_map$cell_id, colnames(exprs(cds)))
message("Fetal cells matched: ", length(keep_in_matrix))
cds2 <- cds[, keep_in_matrix]

# Attaching paper_class labels
rownames(cluster_map) <- cluster_map$cell_id
pData(cds2)$paper_class <- cluster_map[colnames(exprs(cds2)), "paper_class"]

# Adding broad class
broad_map <- c(AP1="AP", AP2="AP", BP1="BP", BP2="BP", N1="N", N2="N", N3="N")
pData(cds2)$cell_class <- broad_map[pData(cds2)$paper_class]

# Setting factor levels
pData(cds2)$paper_class <- factor(pData(cds2)$paper_class,
                                  levels = c("AP1","AP2","BP1","BP2","N1","N2","N3"))
pData(cds2)$cell_class <- factor(pData(cds2)$cell_class,
                                 levels = c("AP","BP","N"))

message("Cells per paper class:")
print(table(pData(cds2)$paper_class, useNA = "ifany"))
message("Cells per broad class:")
print(table(pData(cds2)$cell_class, useNA = "ifany"))

## Removing unlabelled cells

labelled_cells <- colnames(exprs(cds2))[!is.na(pData(cds2)$paper_class)]
cds2 <- cds2[, labelled_cells]
message("Cells after removing unlabelled: ", ncol(exprs(cds2)))

## Setting ordering genes by using marker-based genes for ordering to provide an ideal trajectory shape
# PC loading genes are saved as a fallback but identified marker genes produce an ICA space that is closer to original paper

# Loading marker genes from Fig1C script output
marker_genes_df <- read.csv("outputs/fig1C_cluster_markers_top40.csv",
                            stringsAsFactors = FALSE)
ordering_genes_markers <- unique(marker_genes_df$gene)
ordering_genes <- ordering_genes_markers[ordering_genes_markers %in% 
                                           rownames(exprs(cds2))]

message("Marker-based ordering genes matched: ", length(ordering_genes))

# Fallback to PC genes if markers insufficient
if (length(ordering_genes) < 20) {
  message("Falling back to PC loading genes...")
  pc_genes_raw   <- readRDS(gene_file)
  ordering_genes <- pc_genes_raw[pc_genes_raw %in% rownames(exprs(cds2))]
  message("PC ordering genes matched: ", length(ordering_genes))
}

if (length(ordering_genes) < 20) {
  stop("Too few ordering genes remain after subsetting.")
}

cds2 <- setOrderingFilter(cds2, ordering_genes)

## Generating igraph compatibility patch

orig_dfs <- get("dfs", envir = asNamespace("igraph"))
patched_dfs <- function(graph, root,
                        mode = c("out","in","all","total"),
                        neimode = mode, unreachable = TRUE,
                        order = TRUE, order.out = FALSE,
                        father = FALSE, parent = father,
                        dist = FALSE, in.callback = NULL,
                        out.callback = NULL, extra = NULL,
                        rho = parent.frame(), ...) {
  mode   <- if (!missing(neimode)) neimode else mode
  parent <- if (!missing(father))  father  else parent
  orig_dfs(graph=graph, root=root, mode=mode, unreachable=unreachable,
           order=order, order.out=order.out, parent=parent, dist=dist,
           in.callback=in.callback, out.callback=out.callback,
           extra=extra, rho=rho, ...)
}
unlockBinding("dfs", asNamespace("igraph"))
assign("dfs", patched_dfs, envir = asNamespace("igraph"))
lockBinding("dfs", asNamespace("igraph"))

## Performing ICA and orderCells
# orderCells assigns pseudotime values to each cell by traversing
# the minimum spanning tree from the root state outward.
# tryCatch prevents the script from stopping if a stack overflow occurs.

set.seed(42)
cds2 <- reduceDimension(cds2, method = "ICA", max_components = 2)
message("reduceDimension complete")

cds2 <- tryCatch(
  orderCells(cds2),
  error = function(e) {
    message("orderCells() failed: ", conditionMessage(e))
    cds2
  }
)

# Re-rooting at the state most enriched for AP cells
if ("State" %in% colnames(pData(cds2))) {
  state_tab  <- table(pData(cds2)$State, pData(cds2)$cell_class)
  
  if ("AP" %in% colnames(state_tab)) {
    ap_prop    <- state_tab[, "AP"] / rowSums(state_tab)
    root_state <- names(which.max(ap_prop))
    message("Re-rooting at state ", root_state)
    cds2 <- orderCells(cds2, root_state = root_state)
  }
}

message("Pseudotime assigned: ", "Pseudotime" %in% colnames(pData(cds2)))
print(table(pData(cds2)$paper_class, useNA = "ifany"))

## Saving the figures

# Broad AP/BP/N plot (closest to paper Fig 2A)
p <- plot_cell_trajectory(cds2, color_by = "cell_class") +
  scale_color_manual(
    values = c("AP" = "#5F8F4A", "BP" = "#7DC9C9", "N" = "#2E4E9B"),
    na.value = "grey70",
    name = "Cell class"
  ) +
  ggtitle("Fig 2A — AP/BP/Neuron lineage (recreated)")

png("figures/Fig2A_broad_classes.png", width = 7, height = 5,
    units = "in", res = 300)
print(p)
dev.off()
message("Saved figures/Fig2A_broad_classes.png")

# 7-subtype version
p7 <- plot_cell_trajectory(cds2, color_by = "paper_class") +
  scale_color_manual(
    values = c(
      "AP1" = "#5F8F4A", "AP2" = "#9BC67E",
      "BP1" = "#7DC9C9", "BP2" = "#BFE6E6",
      "N1"  = "#9DB8E6", "N2"  = "#5F86C9", "N3" = "#2E4E9B"
    ),
    na.value = "grey70",
    name = "Cell type"
  ) +
  ggtitle("Monocle trajectory — fetal neocortex (Camp et al. data)") +
  xlab("Component 1") + ylab("Component 2") +
  theme_classic(base_size = 12)

png("figures/Fig2A_final_original_data.png", width = 7, height = 5,
    units = "in", res = 300)
print(p7)
dev.off()
message("Saved figures/Fig2A_final_original_data.png")
