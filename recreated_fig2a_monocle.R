### Recreating Fig 2A using the preprocessed Fig 1C classification genes and clusters from recreated_fig1c_genes.R
# Monocle trajectory — recreated analysis data
# Based on recreate_fig2A_from_fig1C_genes.R (working script)

## Applying dplur compatibility patch before loading required packages

select_ <- function(.data, ...) {
  dots <- unlist(list(...))
  dplyr::select(.data, dplyr::all_of(dots))
}
group_by_ <- function(.data, ..., add = FALSE) {
  dots <- unlist(list(...))
  dplyr::group_by(.data, dplyr::across(dplyr::all_of(dots)), .add = add)
}
filter_ <- function(.data, ...) {
  dots <- list(...)
  dplyr::filter(.data, !!!lapply(dots, rlang::parse_expr))
}
mutate_ <- function(.data, ...) {
  dots <- list(...)
  dplyr::mutate(.data, !!!lapply(dots, rlang::parse_expr))
}
arrange_ <- function(.data, ...) {
  dots <- list(...)
  dplyr::arrange(.data, !!!lapply(dots, rlang::parse_expr))
}
summarise_ <- function(.data, ...) {
  dots <- list(...)
  dplyr::summarise(.data, !!!lapply(dots, rlang::parse_expr))
}
tryCatch({
  assignInNamespace("select_",    select_,    ns = "dplyr")
  assignInNamespace("group_by_",  group_by_,  ns = "dplyr")
  assignInNamespace("filter_",    filter_,    ns = "dplyr")
  assignInNamespace("mutate_",    mutate_,    ns = "dplyr")
  assignInNamespace("arrange_",   arrange_,   ns = "dplyr")
  assignInNamespace("summarise_", summarise_, ns = "dplyr")
  message("dplyr patch applied")
}, error = function(e) {
  message("dplyr patch blocked - global definitions in place")
})

## Loading required packages
suppressPackageStartupMessages({
  library(monocle)
  library(Matrix)
  library(igraph)
  library(ggplot2)
  library(org.Hs.eg.db)
})

## Applying igraph dfs() patch

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

## Patching monocle nei() calls

monocle_fns_to_patch <- c("assign_cell_lineage", "buildBranchCellDataSet",
                          "count_leaf_descendents", "cth_classifier_cds",
                          "cth_classifier_cell", "extract_good_branched_ordering",
                          "extract_good_ordering", "make_canonical",
                          "measure_diameter_path", "project2MST")

for (fn_name in monocle_fns_to_patch) {
  tryCatch({
    fn      <- get(fn_name, envir = asNamespace("monocle"))
    fn_text <- deparse(body(fn))
    if (!any(grepl("\\bnei\\(", fn_text))) next
    fn_text_fixed <- gsub("\\bnei\\(", "igraph::neighbors(", fn_text)
    if (fn_name == "project2MST") {
      fn_text_fixed <- gsub(
        "igraph::neighbors\\(closest_vertex_names\\[i\\],",
        "igraph::neighbors(dp_mst, closest_vertex_names[i],",
        fn_text_fixed)
    }
    body(fn) <- parse(text = paste(fn_text_fixed, collapse="\n"))[[1]]
    assignInNamespace(fn_name, fn, ns = "monocle")
    message("Patched: ", fn_name)
  }, error = function(e) message("Could not patch ", fn_name))
}

## Creating a custom plot function that uses reducedDimS, draws MST from principal graph

plot_cell_trajectory <- function(cds, color_by = "State",
                                 show_tree = TRUE, ...) {
  reduced <- t(reducedDimS(cds))
  if (ncol(reduced) < 2) stop("reducedDimS has fewer than 2 dimensions")
  
  plot_df <- as.data.frame(reduced[, 1:2, drop = FALSE])
  colnames(plot_df) <- c("Component 1", "Component 2")
  plot_df$cell_id <- rownames(plot_df)
  
  for (col in c("Component 1", "Component 2")) {
    m  <- mean(plot_df[[col]])
    s  <- sd(plot_df[[col]])
    plot_df <- plot_df[abs(plot_df[[col]] - m) <= 3 * s, ]
  }
  message("Cells after outlier removal: ", nrow(plot_df))
  
  pd <- pData(cds)
  if (color_by %in% colnames(pd)) {
    plot_df$color_var <- as.character(pd[rownames(plot_df), color_by])
  } else {
    plot_df$color_var <- "unknown"
  }
  
  p <- ggplot(plot_df,
              aes(x = .data[["Component 1"]],
                  y = .data[["Component 2"]],
                  colour = color_var)) +
    geom_point(size = 1.5) +
    labs(colour = color_by, x = "Component 1", y = "Component 2") +
    theme_classic(base_size = 12)
  
  if (show_tree && !is.null(cds@minSpanningTree)) {
    tryCatch({
      mst       <- cds@minSpanningTree
      edge_list <- igraph::as_edgelist(mst)
      
      if (length(cds@reducedDimK) > 0) {
        graph_coords <- t(cds@reducedDimK)
        if (nrow(graph_coords) >= 2) {
          colnames(graph_coords) <- c("Component 1", "Component 2")
          rownames(graph_coords) <- igraph::V(mst)$name
          edge_rows <- lapply(seq_len(nrow(edge_list)), function(i) {
            n1 <- edge_list[i, 1]; n2 <- edge_list[i, 2]
            if (n1 %in% rownames(graph_coords) && n2 %in% rownames(graph_coords)) {
              data.frame(x=graph_coords[n1,"Component 1"],
                         y=graph_coords[n1,"Component 2"],
                         xend=graph_coords[n2,"Component 1"],
                         yend=graph_coords[n2,"Component 2"],
                         stringsAsFactors=FALSE)
            } else NULL
          })
          edge_df <- do.call(rbind, Filter(Negate(is.null), edge_rows))
          if (!is.null(edge_df) && nrow(edge_df) > 0) {
            p <- p + ggplot2::geom_segment(
              data=edge_df, aes(x=x,y=y,xend=xend,yend=yend),
              colour="black", linewidth=0.75, inherit.aes=FALSE)
            message("Principal graph drawn: ", nrow(edge_df), " edges")
            return(p)
          }
        }
      }
      
      # Fall back to cell-based MST edges (ICA)
      edge_rows <- lapply(seq_len(nrow(edge_list)), function(i) {
        n1 <- edge_list[i, 1]; n2 <- edge_list[i, 2]
        if (n1 %in% rownames(plot_df) && n2 %in% rownames(plot_df)) {
          data.frame(x=plot_df[n1,"Component 1"], y=plot_df[n1,"Component 2"],
                     xend=plot_df[n2,"Component 1"], yend=plot_df[n2,"Component 2"],
                     stringsAsFactors=FALSE)
        } else NULL
      })
      edge_df <- do.call(rbind, Filter(Negate(is.null), edge_rows))
      if (!is.null(edge_df) && nrow(edge_df) > 0) {
        p <- p + ggplot2::geom_segment(
          data=edge_df, aes(x=x,y=y,xend=xend,yend=yend),
          colour="black", linewidth=0.75, inherit.aes=FALSE)
        message("MST drawn: ", nrow(edge_df), " edges")
      }
    }, error = function(e) message("Could not draw tree: ", conditionMessage(e)))
  }
  p
}

## Assigning input and output files

infile       <- "camp_organoid_matrix_V2.csv"
cluster_file <- "outputs/recreated_fig1C_cell_clusters.csv"
markers_file <- "outputs/recreated_fig1C_cluster_markers_top40.csv"

if (!file.exists(infile))       stop("Missing: ", infile)
if (!file.exists(cluster_file)) stop("Missing: ", cluster_file,
                                     "\nRun Fig2_recreated_fig1C.R first.")
if (!file.exists(markers_file)) stop("Missing: ", markers_file,
                                     "\nRun Fig2_recreated_fig1C.R first.")

if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("outputs")) dir.create("outputs")

## Loading the recreated analysis pipeline data matrix

message("Reading recreated matrix...")
full_data <- read.csv(infile, header=TRUE, row.names=1, check.names=FALSE)
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
rm(full_data)

## Filtering to identify fetal cells

message("Filtering to fetal cells...")
srr_nums  <- as.integer(sub("SRR", "", colnames(full_data_mapped)))
fetal_srr <- colnames(full_data_mapped)[srr_nums >= 2967101 &
                                          srr_nums <= 2967326]
message("Fetal cells: ", length(fetal_srr))
fetal_data  <- full_data_mapped[, fetal_srr]
expr_matrix <- as(as.matrix(fetal_data), "dgCMatrix")
rm(full_data_mapped, fetal_data)

## Building the Monocle CellDataSet

message("Creating Monocle CellDataSet...")
pd <- new("AnnotatedDataFrame",
          data = data.frame(cell_id   = colnames(expr_matrix),
                            row.names = colnames(expr_matrix),
                            stringsAsFactors = FALSE))
fd <- new("AnnotatedDataFrame",
          data = data.frame(gene_short_name = rownames(expr_matrix),
                            row.names       = rownames(expr_matrix),
                            stringsAsFactors = FALSE))

cds <- newCellDataSet(expr_matrix, phenoData=pd, featureData=fd,
                      expressionFamily=negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
message("CellDataSet ready")

## Loading in the assigned clusters

cluster_map <- read.csv(cluster_file, stringsAsFactors=FALSE)
rownames(cluster_map) <- cluster_map$cell_id
message("Cluster map: ", nrow(cluster_map), " cells")

## Filtering and attaching labels

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

message("Cells per class:")
print(table(pData(cds2)$paper_class, useNA="ifany"))

labelled_cells <- colnames(exprs(cds2))[!is.na(pData(cds2)$paper_class)]
cds2 <- cds2[, labelled_cells]
message("Cells after removing unlabelled: ", ncol(exprs(cds2)))

# Ordering genes (capped at 50)

marker_genes_df <- read.csv(markers_file, stringsAsFactors=FALSE)
ordering_genes  <- unique(marker_genes_df$gene)
ordering_genes  <- ordering_genes[ordering_genes %in% rownames(exprs(cds2))]
ordering_genes  <- ordering_genes[1:min(50, length(ordering_genes))]
message("Ordering genes: ", length(ordering_genes))
if (length(ordering_genes) < 20) stop("Too few ordering genes")
cds2 <- setOrderingFilter(cds2, ordering_genes)

## Adding Gaussian jitter to resolve near-identical cell profiles which cause ICA kmeans initialisation to fail with this dataset

message("Adding jitter to expression matrix...")
set.seed(42)
expr_dense    <- as.matrix(exprs(cds2))
expr_jittered <- expr_dense + matrix(
  rnorm(nrow(expr_dense) * ncol(expr_dense), 0, 0.1),
  nrow = nrow(expr_dense))
expr_jittered[expr_jittered < 0] <- 0

cds3 <- cds2
cds3@assayData$exprs <- as(expr_jittered, "dgCMatrix")
message("Jitter applied")

## Performing ICA dimension reduction

message("Running ICA...")
options(expressions = 50000)
set.seed(42)
cds3 <- reduceDimension(cds3, method = "ICA", max_components = 2)
message("ICA complete")
message("reducedDimS: ", paste(dim(reducedDimS(cds3)), collapse = " x "))

## Applying orderCells and re-rooting at AP1
# orderCells assigns pseudotime values to each cell by traversing
# the minimum spanning tree from the root state outward.
# tryCatch prevents the script from stopping if a stack overflow occurs.

message("Ordering cells...")
cds3 <- tryCatch(
  orderCells(cds3),
  error = function(e) {
    message("orderCells failed: ", conditionMessage(e)); cds3
  }
)

if ("State" %in% colnames(pData(cds3))) {
  ap1_cells <- colnames(exprs(cds3))[
    !is.na(pData(cds3)$paper_class) &
      as.character(pData(cds3)$paper_class) == "AP1"]
  if (length(ap1_cells) > 0) {
    root_state <- names(which.max(table(pData(cds3)[ap1_cells, "State"])))
    message("Re-rooting at state ", root_state)
    cds3 <- tryCatch(
      orderCells(cds3, root_state = root_state),
      error = function(e) { message("Re-root failed: ", conditionMessage(e)); cds3 }
    )
  }
}

message("Pseudotime assigned: ", "Pseudotime" %in% colnames(pData(cds3)))
print(table(pData(cds3)$paper_class, useNA = "ifany"))

## Saving figure outputs

seven_colours <- c("AP1"="#5F8F4A", "AP2"="#9BC67E",
                   "BP1"="#7DC9C9", "BP2"="#BFE6E6",
                   "N1"="#9DB8E6",  "N2"="#5F86C9", "N3"="#2E4E9B")

p7 <- plot_cell_trajectory(cds3, color_by = "paper_class") +
  scale_color_manual(values = seven_colours, na.value = "grey70",
                     name = "Cell type") +
  ggtitle("Monocle trajectory - fetal neocortex (recreated analysis)") +
  xlab("Component 1") + ylab("Component 2") +
  theme_classic(base_size = 12)

png("figures/Fig2_recreated_seven_classes.png", width = 7, height = 5,
    units = "in", res = 300)
print(p7)
dev.off()
message("Saved figures/Fig2_recreated_seven_classes.png")

saveRDS(cds3, "outputs/monocle_cds_fig2_recreated.rds")
message("Saved outputs/monocle_cds_fig2_recreated.rds")
