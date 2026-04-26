### Recreating Fig 2A using the preprocessed Fig 1C classification genes and clusters from star_fig1c_genes.R
# Monocle trajectory — STAR alternate (post-2016) analysis data
# Based on recreate_fig2A_from_fig1C_genes.R (working script)
# INPUTS:
#   deseq2_counts_matrix.tsv
#   outputs/STAR_fig1C_cell_clusters.csv
#   outputs/STAR_fig1C_cluster_markers_top40.csv
#
# OUTPUTS:
#   figures/Fig3_STAR_broad_classes.png
#   figures/Fig3_STAR_seven_classes.png

## Loading required packages

suppressPackageStartupMessages({
  library(monocle)
  library(Matrix)
  library(igraph)
  library(ggplot2)
  library(org.Hs.eg.db)
})

## Assigning input and output files

infile       <- "deseq2_counts_matrix.tsv"
cluster_file <- "outputs/STAR_fig1C_cell_clusters.csv"
markers_file <- "outputs/STAR_fig1C_cluster_markers_top40.csv"

if (!file.exists(infile))       stop("Missing: ", infile)
if (!file.exists(cluster_file)) stop("Missing: ", cluster_file,
                                     "\nRun star_fig1c_genes_FIXED.R first.")
if (!file.exists(markers_file)) stop("Missing: ", markers_file,
                                     "\nRun star_fig1c_genes_FIXED.R first.")

if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("outputs")) dir.create("outputs")

## Loading in the STAR data matrix

message("Reading STAR matrix...")
full_data <- read.table(infile,
                        header           = TRUE,
                        sep              = "\t",
                        check.names      = FALSE,
                        row.names        = 1)

if ("gene_id" %in% colnames(full_data)) {
  full_data <- full_data[, colnames(full_data) != "gene_id"]
}
message("Raw: ", nrow(full_data), " genes x ", ncol(full_data), " cells")

## Converting Ensembl IDs to gene symbols

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
message("After conversion: ",
        nrow(full_data_mapped), " genes x ", ncol(full_data_mapped), " cells")

## Filtering to identify fetal cells

message("Filtering to fetal cells...")
srr_nums  <- as.integer(sub("SRR", "", colnames(full_data_mapped)))
fetal_srr <- colnames(full_data_mapped)[srr_nums >= 2967101 &
                                          srr_nums <= 2967326]
message("Fetal cells: ", length(fetal_srr))
fetal_data <- full_data_mapped[, fetal_srr]

expr_matrix <- as(as.matrix(fetal_data), "dgCMatrix")
rm(full_data, full_data_mapped, fetal_data)

## Building the Monocle CellDataSet

message("Creating Monocle CellDataSet...")
pd_df <- data.frame(cell_id   = colnames(expr_matrix),
                    row.names = colnames(expr_matrix),
                    stringsAsFactors = FALSE)
fd_df <- data.frame(gene_short_name = rownames(expr_matrix),
                    row.names       = rownames(expr_matrix),
                    stringsAsFactors = FALSE)

pd <- new("AnnotatedDataFrame", data = pd_df)
fd <- new("AnnotatedDataFrame", data = fd_df)

cds <- newCellDataSet(expr_matrix,
                      phenoData        = pd,
                      featureData      = fd,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
message("CellDataSet ready")

## Loading in STAR cluster assignments

cluster_map <- read.csv(cluster_file, stringsAsFactors = FALSE)
rownames(cluster_map) <- cluster_map$cell_id

message("STAR cluster map loaded: ", nrow(cluster_map), " cells")
message("Cells matched in matrix: ",
        length(intersect(colnames(exprs(cds)), cluster_map$cell_id)))

## Filtering again to fetal cells and attaching labels

keep_in_matrix <- intersect(cluster_map$cell_id, colnames(exprs(cds)))
message("Fetal cells matched: ", length(keep_in_matrix))
cds2 <- cds[, keep_in_matrix]

pData(cds2)$paper_class <- cluster_map[colnames(exprs(cds2)), "paper_class"]

broad_map <- c(AP1="AP", AP2="AP", BP1="BP", BP2="BP",
               N1="N",  N2="N",  N3="N")
pData(cds2)$cell_class <- broad_map[pData(cds2)$paper_class]

pData(cds2)$paper_class <- factor(pData(cds2)$paper_class,
                                  levels = c("AP1","AP2","BP1","BP2","N1","N2","N3"))
pData(cds2)$cell_class  <- factor(pData(cds2)$cell_class,
                                  levels = c("AP","BP","N"))

message("Cells per paper class:")
print(table(pData(cds2)$paper_class, useNA = "ifany"))
message("Cells per broad class:")
print(table(pData(cds2)$cell_class, useNA = "ifany"))

## Removing unlabelled cells

labelled_cells <- colnames(exprs(cds2))[!is.na(pData(cds2)$paper_class)]
cds2 <- cds2[, labelled_cells]
message("Cells after removing unlabelled: ", ncol(exprs(cds2)))

## Setting the ordering genes

marker_genes_df <- read.csv(markers_file, stringsAsFactors = FALSE)
ordering_genes  <- unique(marker_genes_df$gene)
ordering_genes  <- ordering_genes[ordering_genes %in% rownames(exprs(cds2))]

message("Ordering genes matched: ", length(ordering_genes))
if (length(ordering_genes) < 20) stop("Too few ordering genes.")

cds2 <- setOrderingFilter(cds2, ordering_genes)

## Applying an igraph compatibility patch

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

## Performing ICA dimension reduction
## Applying orderCells and re-rooting at AP1
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

## Re-rooting specifically at AP1 state

if ("State" %in% colnames(pData(cds2))) {
  
  # rooting at state with most AP1 cells
  ap1_cells <- colnames(exprs(cds2))[
    !is.na(pData(cds2)$paper_class) &
      as.character(pData(cds2)$paper_class) == "AP1"]
  
  if (length(ap1_cells) > 0) {
    state_ap1  <- pData(cds2)[ap1_cells, "State"]
    root_state <- names(which.max(table(state_ap1)))
    message("Re-rooting at state ", root_state, " (most AP1 cells)")
    cds2 <- tryCatch(
      orderCells(cds2, root_state = root_state),
      error = function(e) {
        message("AP1 re-rooting failed: ", conditionMessage(e)); cds2
      }
    )
  }
  
  # if AP1 rooting does not work, use most AP-enriched state
  if (!"Pseudotime" %in% colnames(pData(cds2)) ||
      all(pData(cds2)$Pseudotime == 0)) {
    state_tab  <- table(pData(cds2)$State, pData(cds2)$cell_class)
    if ("AP" %in% colnames(state_tab)) {
      ap_prop    <- state_tab[, "AP"] / rowSums(state_tab)
      root_state <- names(which.max(ap_prop))
      message("Fallback: re-rooting at state ", root_state,
              " (most AP cells broadly)")
      cds2 <- tryCatch(
        orderCells(cds2, root_state = root_state),
        error = function(e) {
          message("Fallback re-rooting failed: ", conditionMessage(e)); cds2
        }
      )
    }
  }
}

message("Pseudotime assigned: ", "Pseudotime" %in% colnames(pData(cds2)))
print(table(pData(cds2)$paper_class, useNA = "ifany"))

## Saving required figure outputs

# Broad AP/BP/N
p_broad <- plot_cell_trajectory(cds2, color_by = "cell_class") +
  scale_color_manual(
    values = c("AP" = "#5F8F4A", "BP" = "#7DC9C9", "N" = "#2E4E9B"),
    na.value = "grey70",
    name = "Cell class"
  ) +
  ggtitle("Monocle trajectory — fetal neocortex (STAR alternative pipeline)")

png("figures/Fig3_STAR_broad_classes.png", width = 7, height = 5,
    units = "in", res = 300)
print(p_broad)
dev.off()
message("Saved figures/Fig3_STAR_broad_classes.png")

# 7-subtype version
p_seven <- plot_cell_trajectory(cds2, color_by = "paper_class") +
  scale_color_manual(
    values = c("AP1"="#5F8F4A", "AP2"="#9BC67E",
               "BP1"="#7DC9C9", "BP2"="#BFE6E6",
               "N1"="#9DB8E6",  "N2"="#5F86C9", "N3"="#2E4E9B"),
    na.value = "grey70",
    name = "Cell type"
  ) +
  ggtitle("Monocle trajectory — fetal neocortex (STAR alternative pipeline)") +
  xlab("Component 1") + ylab("Component 2") +
  theme_classic(base_size = 12)

png("figures/Fig3_STAR_seven_classes.png", width = 7, height = 5,
    units = "in", res = 300)
print(p_seven)
dev.off()
message("Saved figures/Fig3_STAR_seven_classes.png")

