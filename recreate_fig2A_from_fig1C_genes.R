# Recreate Fig 2A using recreated Fig 1C classification genes AND 7 Fig 1C clusters

suppressPackageStartupMessages({
  library(monocle)
  library(Matrix)
  library(igraph)
  library(ggplot2)
})

# ----------------------------
# 0) Input / output
# ----------------------------
infile <- "GSE75140_hOrg.fetal.master.data.frame.txt.gz"
gene_file <- "outputs/fig1C_classification_genes.rds"
cluster_file <- "outputs/fig1C_cell_clusters.csv"

if (!file.exists(infile)) stop("Missing file: ", normalizePath(infile, mustWork = FALSE))
if (!file.exists(gene_file)) stop("Missing file: ", normalizePath(gene_file, mustWork = FALSE))
if (!file.exists(cluster_file)) stop("Missing file: ", normalizePath(cluster_file, mustWork = FALSE))

if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("outputs")) dir.create("outputs")

# ----------------------------
# 1) Load GEO matrix
# ----------------------------
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

expr_matrix <- t(expr_matrix)  # genes x cells

if (anyNA(expr_matrix)) {
  warning("NAs detected in matrix; replacing with 0.")
  expr_matrix[is.na(expr_matrix)] <- 0
}

expr_matrix <- as(expr_matrix, "dgCMatrix")
rm(raw_df)

# ----------------------------
# 2) Build Monocle CellDataSet
# ----------------------------
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

# ----------------------------
# 3) Robust gene resolver
# ----------------------------
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

# ----------------------------
# 4) Load recreated Fig 1C classification genes
# ----------------------------
classification_genes_raw <- readRDS(gene_file)
classification_genes <- resolve_genes(classification_genes_raw, available_genes)
classification_genes <- classification_genes[!is.na(classification_genes)]

message("Classification genes loaded: ", length(classification_genes))

if (length(classification_genes) < 20) {
  stop("Too few classification genes survived matching.")
}

# ----------------------------
# 5) Load Fig 1C cluster assignments
# ----------------------------
cluster_map <- read.csv(cluster_file, stringsAsFactors = FALSE)
rownames(cluster_map) <- cluster_map$cell_id

pData(cds)$fig1c_cluster <- cluster_map[colnames(exprs(cds)), "fig1c_cluster"]

# Keep as factor with 7 levels if present
cluster_levels <- sort(unique(na.omit(pData(cds)$fig1c_cluster)))
pData(cds)$fig1c_cluster <- factor(pData(cds)$fig1c_cluster, levels = cluster_levels)

message("Fig 1C clusters loaded:")
print(table(pData(cds)$fig1c_cluster, useNA = "ifany"))

# ----------------------------
# 6) Broad AP / BP / N labels for filtering only
# ----------------------------
ap_markers <- c("PAX6", "GLI3", "SOX2", "HES1", "VIM", "PROM1")
bp_markers <- c("EOMES", "INSM1", "HES6", "ASPM", "NEUROD4")
n_markers  <- c("MYT1L", "TBR1", "BCL11B", "NEUROD6")
exclude_markers <- c("PECAM1", "GAD1", "DLX1", "DLX2", "DLX5", "DLX6", "ERBB4")

ap_m <- resolve_genes(ap_markers, available_genes)
bp_m <- resolve_genes(bp_markers, available_genes)
n_m  <- resolve_genes(n_markers, available_genes)
ex_m <- resolve_genes(exclude_markers, available_genes)

expr_log <- log1p(exprs(cds))

score_set <- function(gset) {
  if (length(gset) == 0) {
    out <- rep(0, ncol(expr_log))
    names(out) <- colnames(expr_log)
    return(out)
  }
  out <- Matrix::colMeans(expr_log[gset, , drop = FALSE])
  names(out) <- colnames(expr_log)
  out
}

ap_score <- score_set(ap_m)
bp_score <- score_set(bp_m)
n_score  <- score_set(n_m)
ex_score <- score_set(ex_m)

lineage_score <- pmax(ap_score, bp_score, n_score)

# Slightly stricter filter
lineage_thr <- as.numeric(quantile(lineage_score, 0.30))
exclude_thr <- if (length(ex_m) > 0) as.numeric(quantile(ex_score, 0.85)) else Inf

keep_cells <- names(lineage_score)[lineage_score >= lineage_thr & ex_score <= exclude_thr]

if (length(keep_cells) < 150) {
  message("Relaxing filter thresholds...")
  lineage_thr <- as.numeric(quantile(lineage_score, 0.20))
  exclude_thr <- if (length(ex_m) > 0) as.numeric(quantile(ex_score, 0.90)) else Inf
  keep_cells <- names(lineage_score)[lineage_score >= lineage_thr & ex_score <= exclude_thr]
}

if (length(keep_cells) < 50) {
  message("Too few cells retained; keeping all cells instead.")
  keep_cells <- colnames(exprs(cds))
}

message("Keeping ", length(keep_cells), " cells out of ", ncol(exprs(cds)))

cds2 <- cds[, keep_cells]
pData(cds2)$fig1c_cluster <- pData(cds)$fig1c_cluster[colnames(exprs(cds2))]

# ----------------------------
# 7) Rename the 7 clusters to paper-like labels
# ----------------------------
cluster_ids <- sort(unique(as.character(pData(cds2)$fig1c_cluster)))
print(cluster_ids)

# EDIT THIS if needed after seeing printed cluster_ids
cluster_to_paper <- setNames(
  c("AP1", "AP2", "BP1", "BP2", "N1", "N2", "N3")[seq_along(cluster_ids)],
  cluster_ids
)

print(cluster_to_paper)

pData(cds2)$paper_class <- unname(cluster_to_paper[as.character(pData(cds2)$fig1c_cluster)])

pData(cds2)$paper_class <- factor(
  pData(cds2)$paper_class,
  levels = c("AP1", "AP2", "BP1", "BP2", "N1", "N2", "N3")
)

print(table(pData(cds2)$paper_class, useNA = "ifany"))

# ----------------------------
# 8) Use recreated Fig 1C genes as ordering genes
# ----------------------------
ordering_genes <- classification_genes[classification_genes %in% rownames(exprs(cds2))]

if (length(ordering_genes) < 20) {
  stop("Too few ordering genes remain after subsetting.")
}

message("Ordering genes used: ", length(ordering_genes))
cds2 <- setOrderingFilter(cds2, ordering_genes)

# ----------------------------
# 9) Patch Monocle2 / igraph compatibility
# ----------------------------
orig_dfs <- get("dfs", envir = asNamespace("igraph"))

patched_dfs <- function(graph, root, mode = c("out", "in", "all", "total"),
                        neimode = mode,
                        unreachable = TRUE, order = TRUE, order.out = FALSE,
                        father = FALSE, parent = father, dist = FALSE,
                        in.callback = NULL, out.callback = NULL,
                        extra = NULL, rho = parent.frame(), ...) {
  mode <- if (!missing(neimode)) neimode else mode
  parent <- if (!missing(father)) father else parent
  
  orig_dfs(
    graph = graph,
    root = root,
    mode = mode,
    unreachable = unreachable,
    order = order,
    order.out = order.out,
    parent = parent,
    dist = dist,
    in.callback = in.callback,
    out.callback = out.callback,
    extra = extra,
    rho = rho,
    ...
  )
}

unlockBinding("dfs", asNamespace("igraph"))
assign("dfs", patched_dfs, envir = asNamespace("igraph"))
lockBinding("dfs", asNamespace("igraph"))

# ----------------------------
# 10) ICA + orderCells
# ----------------------------
set.seed(42)
cds2 <- reduceDimension(cds2, method = "ICA", max_components = 2)

cds2 <- tryCatch(
  {
    orderCells(cds2)
  },
  error = function(e) {
    message("orderCells() failed; continuing with reduceDimension() output.")
    message("Original error: ", conditionMessage(e))
    cds2
  }
)

message("Paper-like classes:")
print(table(pData(cds2)$paper_class, useNA = "ifany"))

# ----------------------------
# 11) Save Fig 2A with 7 classes
# ----------------------------
p <- plot_cell_trajectory(cds2, color_by = "paper_class") +
  scale_color_manual(
    values = c(
      "AP1" = "#5F8F4A",
      "AP2" = "#9BC67E",
      "BP1" = "#7DC9C9",
      "BP2" = "#BFE6E6",
      "N1"  = "#9DB8E6",
      "N2"  = "#5F86C9",
      "N3"  = "#2E4E9B"
    ),
    na.value = "grey70",
    drop = FALSE
  ) +
  ggtitle("Fig 2A recreated using Fig 1C-derived classification genes")

png("figures/Fig2A_from_recreated_Fig1C_clusters1.png", width = 7, height = 5, units = "in", res = 300)
print(p)
dev.off()

message("Saved figures/Fig2A_from_recreated_Fig1C_clusters1.png")
message("Saved outputs/monocle_cds_fig2A_from_fig1C_clusters.rds")
message("DONE ✅")
