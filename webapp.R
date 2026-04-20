# ================================
# Organoid Explorer
# ================================
# Updated to include 5 new tabs:
#   Tab 1: Study Comparison (Fig 3D)
#   Tab 2: Recreated Gene Query
#   Tab 3: Alternate Pipeline Gene Query
#   Tab 4: Monocle Figure Comparison
#   Tab 5: Heatmaps / Alternate Figures
#
# NOTE:
# Search for "# TODO:" comments throughout this file.
# These mark where you need to plug in the data files and figures.
# ================================

setwd("/home/fap16/steered_research_project/r_shiny")
getwd()

library(shiny)
library(bslib)
library(ggplot2)
library(plotly)
library(dplyr)
library(Seurat)

# ================================
# LOAD DATA — ORIGINAL DATASET
# ================================
umap_df    <- readRDS("data/umap_coords.rds")
expr_mat   <- readRDS("data/expression_matrix.rds")
seurat_obj <- readRDS("data/seurat_processed.rds")

# Extract t-SNE
tsne_df <- as.data.frame(Embeddings(seurat_obj, "tsne"))
tsne_df$cell_id <- rownames(tsne_df)
tsne_df$cluster <- Idents(seurat_obj)

# Clean column names
colnames(umap_df) <- tolower(colnames(umap_df))
colnames(tsne_df) <- c("tsne_1", "tsne_2", "cell_id", "cluster")

umap_df$cluster <- as.factor(umap_df$cluster)
tsne_df$cluster <- as.factor(tsne_df$cluster)

genes <- rownames(expr_mat)

# ================================
# LOAD DATA — RECREATED DATASET
# TODO: Replace these with the correct paths to the recreated dataset RDS files
# The recreated dataset uses the same format as the original
# ================================
# recreated_umap_df  <- readRDS("data/recreated_umap_coords.rds")
# recreated_expr_mat <- readRDS("data/recreated_expression_matrix.rds")
# recreated_genes    <- rownames(recreated_expr_mat)

# ================================
# LOAD DATA — STAR ALTERNATIVE PIPELINE DATASET
# TODO: Replace these with the correct paths to the STAR dataset RDS files
# ================================
# star_umap_df  <- readRDS("data/star_umap_coords.rds")
# star_expr_mat <- readRDS("data/star_expression_matrix.rds")
# star_genes    <- rownames(star_expr_mat)

# ================================
# UI
# ================================
ui <- fluidPage(
  
  theme = bs_theme(version = 5, bootswatch = "darkly"),
  
  titlePanel("Cerebral Organoid scRNA-seq Explorer"),
  
  tabsetPanel(
    
    # ----------------------------
    # EXISTING TAB: Gene query (original dataset)
    # ----------------------------
    tabPanel(
      "Gene Query (Original)",
      br(),
      
      fluidRow(
        
        column(
          4,
          card(
            card_header("Search"),
            
            selectizeInput(
              "gene",
              "Gene symbol",
              choices  = genes,
              selected = genes[1]
            ),
            
            radioButtons(
              "plot_type",
              "Plot type",
              choices  = c("UMAP" = "umap", "Violin" = "violin"),
              selected = "umap"
            ),
            
            actionButton("go", "Search")
          )
        ),
        
        column(
          8,
          card(
            card_header("Results"),
            
            fluidRow(
              
              column(
                7,
                plotOutput("gene_plot", height = "420px")
              ),
              
              column(
                5,
                div(
                  style = "height:420px; overflow-y:auto;",
                  tableOutput("gene_table")
                )
              )
            )
          )
        )
      )
    ),
    
    # ----------------------------
    # EXISTING TAB: UMAP clusters
    # ----------------------------
    tabPanel(
      "UMAP Clusters",
      plotlyOutput("umap_plot", height = "600px")
    ),
    
    # ----------------------------
    # EXISTING TAB: t-SNE clusters
    # ----------------------------
    tabPanel(
      "t-SNE Clusters",
      plotlyOutput("tsne_plot", height = "600px")
    ),
    
    # ----------------------------
    # NEW TAB 1: Study Comparison — Fig 3D
    # Shows the two comparison figures from the original paper analysis
    # TODO: Replace the img src paths with the correct paths to your
    #       Fig 3D comparison figures once they are saved to the www/ folder.
    #       Shiny serves static files from a folder called www/ in the same
    #       directory as app.R. Copy your PNG files there.
    # ----------------------------
    tabPanel(
      "Study Comparison",
      br(),
      
      h4("Comparison of Cell Composition: Fetal Neocortex vs Cerebral Organoids"),
      p("The figures below show the comparison of single-cell transcriptome 
        data between human fetal neocortex samples and cerebral organoids, 
        replicating the analyses from Camp et al. (2015)."),
      br(),
      
      fluidRow(
        
        column(
          6,
          card(
            card_header("Fetal Neocortex — Fig 3D"),
            # TODO: Replace 'placeholder_fig3D_fetal.png' with your actual
            #       figure filename in the www/ folder
            tags$img(
              src   = "placeholder_fig3D_fetal.png",
              style = "width:100%; height:auto;"
            ),
            p("Figure 3D: Cell clustering from fetal neocortex samples.",
              style = "font-size:0.85em; color:#aaa; margin-top:8px;")
          )
        ),
        
        column(
          6,
          card(
            card_header("Cerebral Organoids — Fig 3D"),
            # TODO: Replace 'placeholder_fig3D_organoid.png' with your actual
            #       figure filename in the www/ folder
            tags$img(
              src   = "placeholder_fig3D_organoid.png",
              style = "width:100%; height:auto;"
            ),
            p("Figure 3D: Cell clustering from cerebral organoid samples.",
              style = "font-size:0.85em; color:#aaa; margin-top:8px;")
          )
        )
      )
    ),
    
    # ----------------------------
    # NEW TAB 2: Recreated Gene Query
    # Same layout as the original gene query but uses the recreated dataset
    # TODO: Uncomment the server section for this tab once the recreated
    #       RDS files are available (search for "recreated_gene_data" below)
    # ----------------------------
    tabPanel(
      "Gene Query (Recreated)",
      br(),
      
      # Placeholder notice — remove once data is loaded
      div(
        class = "alert alert-warning",
        style = "margin:20px;",
        strong("Note: "),
        "This tab requires the recreated dataset RDS files. ",
        "Once available, uncomment the data loading lines at the top of app.R ",
        "and the server section labelled 'recreated_gene_data'."
      ),
      
      fluidRow(
        
        column(
          4,
          card(
            card_header("Search"),
            
            # TODO: Change choices to recreated_genes once data is loaded
            selectizeInput(
              "recreated_gene",
              "Gene symbol",
              choices  = genes,   # TODO: replace with recreated_genes
              selected = genes[1] # TODO: replace with recreated_genes[1]
            ),
            
            radioButtons(
              "recreated_plot_type",
              "Plot type",
              choices  = c("UMAP" = "umap", "Violin" = "violin"),
              selected = "umap"
            ),
            
            actionButton("recreated_go", "Search")
          )
        ),
        
        column(
          8,
          card(
            card_header("Results"),
            
            fluidRow(
              
              column(
                7,
                plotOutput("recreated_gene_plot", height = "420px")
              ),
              
              column(
                5,
                div(
                  style = "height:420px; overflow-y:auto;",
                  tableOutput("recreated_gene_table")
                )
              )
            )
          )
        )
      )
    ),
    
    # ----------------------------
    # NEW TAB 3: STAR Alternative Pipeline Gene Query
    # Same layout as the original gene query but uses the STAR dataset
    # TODO: Uncomment the server section for this tab once the STAR
    #       RDS files are available (search for "star_gene_data" below)
    # ----------------------------
    tabPanel(
      "Gene Query (STAR Pipeline)",
      br(),
      
      # Placeholder notice — remove once data is loaded
      div(
        class = "alert alert-warning",
        style = "margin:20px;",
        strong("Note: "),
        "This tab requires the STAR alternative pipeline dataset RDS files. ",
        "Once available, uncomment the data loading lines at the top of app.R ",
        "and the server section labelled 'star_gene_data'."
      ),
      
      fluidRow(
        
        column(
          4,
          card(
            card_header("Search"),
            
            # TODO: Change choices to star_genes once data is loaded
            selectizeInput(
              "star_gene",
              "Gene symbol",
              choices  = genes,   # TODO: replace with star_genes
              selected = genes[1] # TODO: replace with star_genes[1]
            ),
            
            radioButtons(
              "star_plot_type",
              "Plot type",
              choices  = c("UMAP" = "umap", "Violin" = "violin"),
              selected = "umap"
            ),
            
            actionButton("star_go", "Search")
          )
        ),
        
        column(
          8,
          card(
            card_header("Results"),
            
            fluidRow(
              
              column(
                7,
                plotOutput("star_gene_plot", height = "420px")
              ),
              
              column(
                5,
                div(
                  style = "height:420px; overflow-y:auto;",
                  tableOutput("star_gene_table")
                )
              )
            )
          )
        )
      )
    ),
    
    # ----------------------------
    # NEW TAB 4: Monocle Figure Comparison
    # Shows all three Monocle trajectory figures as static images
    # Static images keep the app lightweight — no RDS files needed
    # TODO: Copy your three Monocle PNG figures to the www/ folder:
    #         www/Fig2A_final_original_data.png     (Figure 1 — original)
    #         www/Fig2_recreated_seven_classes.png  (Figure 2 — recreated)
    #         www/Fig3_STAR_seven_classes.png       (Figure 3 — STAR)
    # ----------------------------
    tabPanel(
      "Monocle Trajectories",
      br(),
      
      h4("Pseudotemporal Trajectory Analysis — Comparison Across Pipelines"),
      p("Monocle 2 was used to reconstruct the apical progenitor (AP) to 
        basal progenitor (BP) to neuron differentiation lineage across three 
        datasets. ICA dimensionality reduction was applied and cells were 
        ordered along a minimum spanning tree."),
      br(),
      
      fluidRow(
        
        column(
          4,
          card(
            card_header("Figure 1 — Original GEO Data (Camp et al.)"),
            tags$img(
              src   = "Fig2A_final_original_data.png",
              style = "width:100%; height:auto;"
            ),
            p("Monocle 2 trajectory reconstructed from the original GEO 
              supplementary count data (GSE75140). Seven cell populations 
              resolved: AP1, AP2, BP1, BP2, N1, N2, N3.",
              style = "font-size:0.85em; color:#aaa; margin-top:8px;")
          )
        ),
        
        column(
          4,
          card(
            card_header("Figure 2 — Recreated Analysis"),
            tags$img(
              src   = "Fig2_recreated_seven_classes.png",
              style = "width:100%; height:auto;"
            ),
            p("Monocle 2 trajectory reconstructed from the recreated analysis 
              dataset. Four cell populations resolved: AP1, AP2, BP1, N1. 
              Reduced subtype resolution reflects differences in upstream 
              data processing.",
              style = "font-size:0.85em; color:#aaa; margin-top:8px;")
          )
        ),
        
        column(
          4,
          card(
            card_header("Figure 3 — STAR Alternative Pipeline"),
            tags$img(
              src   = "Fig3_STAR_seven_classes.png",
              style = "width:100%; height:auto;"
            ),
            p("Monocle 2 trajectory reconstructed from STAR-aligned raw 
              counts (SRP066834). Five cell populations resolved: AP1, AP2, 
              BP1, BP2, N1. Raw count data produced lower subtype resolution 
              than pre-processed FPKM values.",
              style = "font-size:0.85em; color:#aaa; margin-top:8px;")
          )
        )
      )
    ),
    
    # ----------------------------
    # NEW TAB 5: Heatmaps / Alternate Figures
    # Placeholder tab — add heatmap PNG files to www/ and update img src
    # TODO: Decide which heatmaps/figures to show here and copy the
    #       PNG files to the www/ folder, then replace the placeholder
    #       img tags below with the correct filenames
    # ----------------------------
    tabPanel(
      "Heatmaps & Figures",
      br(),
      
      h4("Marker Gene Expression Heatmaps"),
      p("Heatmaps showing canonical marker gene expression across identified 
        cell clusters for each analysis pipeline."),
      br(),
      
      # Placeholder notice — remove once figures are ready
      div(
        class = "alert alert-info",
        style = "margin:20px;",
        strong("Note for group member: "),
        "Add heatmap PNG files to the www/ folder and replace the ",
        "placeholder sections below with the correct filenames. ",
        "Use the same tags$img() format as the Monocle tab above."
      ),
      
      fluidRow(
        
        column(
          6,
          card(
            card_header("Original Dataset — Marker Heatmap"),
            # TODO: Replace with correct filename
            # tags$img(src = "Fig1C_style_marker_heatmap.png",
            #          style = "width:100%; height:auto;")
            div(
              style = "height:300px; display:flex; align-items:center; 
                       justify-content:center; color:#aaa; border:1px dashed #555;",
              "Heatmap figure placeholder — add PNG to www/ folder"
            )
          )
        ),
        
        column(
          6,
          card(
            card_header("STAR Pipeline — Marker Heatmap"),
            # TODO: Replace with correct filename
            # tags$img(src = "STAR_Fig1C_style_marker_heatmap.png",
            #          style = "width:100%; height:auto;")
            div(
              style = "height:300px; display:flex; align-items:center; 
                       justify-content:center; color:#aaa; border:1px dashed #555;",
              "Heatmap figure placeholder — add PNG to www/ folder"
            )
          )
        )
      ),
      
      br(),
      
      fluidRow(
        
        column(
          12,
          card(
            card_header("Additional Figures"),
            # TODO: Add any additional figures here
            div(
              style = "height:200px; display:flex; align-items:center; 
                       justify-content:center; color:#aaa; border:1px dashed #555;",
              "Additional figures placeholder"
            )
          )
        )
      )
    )
  )
)

# ================================
# SERVER
# ================================
server <- function(input, output, session) {
  
  # ----------------------------
  # Original dataset — Gene data
  # ----------------------------
  gene_data <- eventReactive(input$go, {
    
    req(input$gene)
    
    g    <- input$gene
    expr <- expr_mat[g, , drop = TRUE]
    
    df <- data.frame(
      cell_id = colnames(expr_mat),
      expr    = as.numeric(expr)
    )
    
    idx <- match(df$cell_id, umap_df$cell_id)
    
    df$umap_1  <- umap_df$umap_1[idx]
    df$umap_2  <- umap_df$umap_2[idx]
    df$cluster <- umap_df$cluster[idx]
    
    df <- df[complete.cases(df), ]
    df
  })
  
  output$gene_plot <- renderPlot({
    
    df <- gene_data()
    req(df)
    
    g <- input$gene
    
    if (input$plot_type == "violin") {
      
      ggplot(df, aes(cluster, expr)) +
        geom_violin(fill = "steelblue", alpha = 0.7) +
        geom_boxplot(width = 0.15, outlier.size = 0.3) +
        theme_minimal() +
        labs(title = paste("Expression of", g),
             x = "Cluster", y = "Expression")
      
    } else {
      
      df <- df[order(df$expr), ]
      
      ggplot(df, aes(umap_1, umap_2)) +
        geom_point(aes(color = expr), size = 1.5) +
        scale_color_viridis_c(option = "magma") +
        theme_minimal() +
        labs(title = paste("UMAP:", g), color = "Expression")
    }
  })
  
  output$gene_table <- renderTable({
    
    df <- gene_data()
    req(df)
    
    df %>%
      group_by(cluster) %>%
      summarise(avg_expr = mean(expr),
                pct_expr = mean(expr > 0) * 100,
                n_cells  = n(),
                .groups  = "drop")
  })
  
  # ----------------------------
  # Original dataset — UMAP clusters
  # ----------------------------
  output$umap_plot <- renderPlotly({
    
    p <- ggplot(umap_df, aes(umap_1, umap_2, color = cluster)) +
      geom_point(size = 1.5) +
      theme_minimal()
    
    ggplotly(p)
  })
  
  # ----------------------------
  # Original dataset — t-SNE clusters
  # ----------------------------
  output$tsne_plot <- renderPlotly({
    
    p <- ggplot(tsne_df, aes(tsne_1, tsne_2, color = cluster)) +
      geom_point(size = 1.5) +
      theme_minimal()
    
    ggplotly(p)
  })
  
  # ----------------------------
  # Recreated dataset — Gene query
  # TODO: Uncomment this entire block once recreated_expr_mat and
  #       recreated_umap_df are loaded at the top of this script
  # ----------------------------
  # recreated_gene_data <- eventReactive(input$recreated_go, {
  #
  #   req(input$recreated_gene)
  #
  #   g    <- input$recreated_gene
  #   expr <- recreated_expr_mat[g, , drop = TRUE]
  #
  #   df <- data.frame(
  #     cell_id = colnames(recreated_expr_mat),
  #     expr    = as.numeric(expr)
  #   )
  #
  #   idx <- match(df$cell_id, recreated_umap_df$cell_id)
  #
  #   df$umap_1  <- recreated_umap_df$umap_1[idx]
  #   df$umap_2  <- recreated_umap_df$umap_2[idx]
  #   df$cluster <- recreated_umap_df$cluster[idx]
  #
  #   df <- df[complete.cases(df), ]
  #   df
  # })
  #
  # output$recreated_gene_plot <- renderPlot({
  #
  #   df <- recreated_gene_data()
  #   req(df)
  #
  #   g <- input$recreated_gene
  #
  #   if (input$recreated_plot_type == "violin") {
  #
  #     ggplot(df, aes(cluster, expr)) +
  #       geom_violin(fill = "steelblue", alpha = 0.7) +
  #       geom_boxplot(width = 0.15, outlier.size = 0.3) +
  #       theme_minimal() +
  #       labs(title = paste("Expression of", g),
  #            x = "Cluster", y = "Expression")
  #
  #   } else {
  #
  #     df <- df[order(df$expr), ]
  #
  #     ggplot(df, aes(umap_1, umap_2)) +
  #       geom_point(aes(color = expr), size = 1.5) +
  #       scale_color_viridis_c(option = "magma") +
  #       theme_minimal() +
  #       labs(title = paste("UMAP:", g), color = "Expression")
  #   }
  # })
  #
  # output$recreated_gene_table <- renderTable({
  #
  #   df <- recreated_gene_data()
  #   req(df)
  #
  #   df %>%
  #     group_by(cluster) %>%
  #     summarise(avg_expr = mean(expr),
  #               pct_expr = mean(expr > 0) * 100,
  #               n_cells  = n(),
  #               .groups  = "drop")
  # })
  
  # ----------------------------
  # STAR pipeline dataset — Gene query
  # TODO: Uncomment this entire block once star_expr_mat and
  #       star_umap_df are loaded at the top of this script
  # ----------------------------
  # star_gene_data <- eventReactive(input$star_go, {
  #
  #   req(input$star_gene)
  #
  #   g    <- input$star_gene
  #   expr <- star_expr_mat[g, , drop = TRUE]
  #
  #   df <- data.frame(
  #     cell_id = colnames(star_expr_mat),
  #     expr    = as.numeric(expr)
  #   )
  #
  #   idx <- match(df$cell_id, star_umap_df$cell_id)
  #
  #   df$umap_1  <- star_umap_df$umap_1[idx]
  #   df$umap_2  <- star_umap_df$umap_2[idx]
  #   df$cluster <- star_umap_df$cluster[idx]
  #
  #   df <- df[complete.cases(df), ]
  #   df
  # })
  #
  # output$star_gene_plot <- renderPlot({
  #
  #   df <- star_gene_data()
  #   req(df)
  #
  #   g <- input$star_gene
  #
  #   if (input$star_plot_type == "violin") {
  #
  #     ggplot(df, aes(cluster, expr)) +
  #       geom_violin(fill = "steelblue", alpha = 0.7) +
  #       geom_boxplot(width = 0.15, outlier.size = 0.3) +
  #       theme_minimal() +
  #       labs(title = paste("Expression of", g),
  #            x = "Cluster", y = "Expression")
  #
  #   } else {
  #
  #     df <- df[order(df$expr), ]
  #
  #     ggplot(df, aes(umap_1, umap_2)) +
  #       geom_point(aes(color = expr), size = 1.5) +
  #       scale_color_viridis_c(option = "magma") +
  #       theme_minimal() +
  #       labs(title = paste("UMAP:", g), color = "Expression")
  #   }
  # })
  #
  # output$star_gene_table <- renderTable({
  #
  #   df <- star_gene_data()
  #   req(df)
  #
  #   df %>%
  #     group_by(cluster) %>%
  #     summarise(avg_expr = mean(expr),
  #               pct_expr = mean(expr > 0) * 100,
  #               n_cells  = n(),
  #               .groups  = "drop")
  # })
}

# ================================
# RUN APP
# ================================
shinyApp(ui, server)