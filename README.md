# Monocle 2 Pseudotemporal Trajectory Analysis

This repository implements a pseudotemporal trajectory analysis pipeline based 
on Camp et al. (2015), aiming to replicate the AP-BP-neuron differentiation 
lineage from Fig 2A using three separate datasets. Where exact replication was 
not possible (e.g. software compatibility constraints, package deprecation), 
programmatic patches and reasonable approximations were implemented while 
maintaining consistency with the original study.

The workflow consists of two main stages applied to three datasets:

    1. Cell clustering and biological classification using Seurat v5
    2. Pseudotemporal trajectory inference using Monocle 2

The deployed web application presenting the results is accessible at:
https://www122.lamp.le.ac.uk/

---

## Fig1C-style Preprocessing — Original GEO Data  `recreate_fig1C_genes.R`

Performs cell clustering and biological classification of the original GEO 
supplementary count data (GSE75140) using Seurat v5, replicating the approach 
of Camp et al. Fig 1C. Raw count data is loaded directly with gene symbols as 
row identifiers. Fetal cells are identified by cell ID pattern matching for 
strings containing "wpc", "fetal" or "ctx". Data is log-normalised, the top 
2,500 variable features are identified by VST, and Louvain clustering is 
performed at resolution 2.0 using the first 10 principal components. Clusters 
are assigned to AP, BP and N lineages using a marker score-based approach 
incorporating a BP specificity index. Top 40 marker genes per cluster ranked 
by average log2 fold change are saved as outputs for the Monocle trajectory 
script.

---

## Monocle 2 Trajectory — Original GEO Data  `recreate_fig2A_from_fig1C_genes.R`

Performs pseudotemporal trajectory analysis on the original GEO dataset, 
replicating the approach of Camp et al. Fig 2A. Takes the cluster map and 
marker genes from recreate_fig1C_genes.R as inputs. A Monocle 2 CellDataSet 
object is constructed using the negative binomial size distribution. ICA with 
2 components is used for dimensionality reduction. Cells are ordered along the 
trajectory using the minimum spanning tree algorithm, with the root set to the 
Monocle state most enriched for AP1 cells. Programmatic compatibility patches 
for dplyr and igraph are applied at runtime. Generates:
- `figures/Fig2A_final_original_data.png`

---

## Fig1C-style Preprocessing — Recreated Pre-2016 Dataset  `recreated_fig1c_genes.R`

Performs cell clustering and biological classification of the recreated 
pre-2016 analysis dataset (camp_organoid_matrix_V2.csv) using Seurat v5. 
Ensembl gene IDs are converted to HGNC symbols using org.Hs.eg.db, retaining 
the highest-expressing Ensembl ID where duplicates exist. Fetal cells are 
identified by filtering to SRR accession numbers SRR2967101-SRR2967326. The 
same Seurat clustering pipeline as the GEO preprocessing script is applied. 
Automatic BP cluster detection failed due to weak BP marker expression relative 
to AP markers across all clusters, requiring manual assignment of the cluster 
with the highest BP specificity index to the BP lineage.

---

## Monocle 2 Trajectory — Recreated Pre-2016 Dataset  `recreated_fig2a_monocle.R`

Performs pseudotemporal trajectory analysis on the recreated pre-2016 analysis 
dataset. Takes the cluster map and marker genes from recreated_fig1c_genes.R 
as inputs. Gaussian jitter is applied to the expression matrix prior to ICA 
to resolve a kmeans initialisation error arising from near-identical cell 
expression profiles specific to this dataset. Programmatic compatibility 
patches for dplyr and igraph are applied at runtime. R must be launched with 
an increased stack size to avoid stack overflow errors during orderCells:

```bash
ulimit -s unlimited
R --max-ppsize=500000
```

Generates:
- `figures/Fig2_recreated_seven_classes.png`

---

## Fig1C-style Preprocessing — STAR Alternative Post-2016 Pipeline  `star_fig1c_genes.R`

Performs cell clustering and biological classification of the STAR-aligned 
count matrix (deseq2_counts_matrix.tsv) using Seurat v5. Ensembl gene IDs 
are converted to HGNC symbols using org.Hs.eg.db. Fetal cells are identified 
by filtering to SRR accession numbers SRR2967101-SRR2967326. The same Seurat 
clustering pipeline is applied. The STAR-aligned data contains raw integer 
counts rather than FPKM values, reflecting differences in upstream alignment 
and quantification compared to the GEO and recreated datasets.

---

## Monocle 2 Trajectory — STAR Alternative Post-2016 Pipeline  `star_fig2a_monocle.R`

Performs pseudotemporal trajectory analysis on the STAR-aligned alternative 
pipeline dataset. Takes the cluster map and marker genes from 
star_fig1c_genes.R as inputs. The same Monocle 2 pipeline as the GEO 
trajectory script is applied. The small AP1 cluster (n=6) produced by the 
STAR data prevented reliable pseudotime rooting at the cycling progenitor 
state. Programmatic compatibility patches for dplyr and igraph are applied 
at runtime. Generates:
- `figures/Fig3_STAR_seven_classes.png`

---

## Trial Script  `monocle_trail.R`

An exploratory trial script produced during early development of the Monocle 
2 pipeline. This script was not used in the final analysis and is included 
in the repository for reference only.

---

## Input Data

| File | Source | Used by |
|------|--------|---------|
| GSE75140_hOrg.fetal.master.data.frame.txt.gz | GEO accession GSE75140 | recreate_fig1C_genes.R, recreate_fig2A_from_fig1C_genes.R |
| camp_organoid_matrix_V2.csv | Group GitHub repository | recreated_fig1c_genes.R, recreated_fig2a_monocle.R |
| deseq2_counts_matrix.tsv | SRA project SRP066834 | star_fig1c_genes.R, star_fig2a_monocle.R |

---

## Output Figures

| Figure | Generated by |
|--------|-------------|
| Fig2A_final_original_data.png | recreate_fig2A_from_fig1C_genes.R |
| Fig2_recreated_seven_classes.png | recreated_fig2a_monocle.R |
| Fig3_STAR_seven_classes.png | star_fig2a_monocle.R |

---

## Required R Packages

```r
library(Seurat)        # v5 — Hao et al., 2023
library(monocle)       # v2 — Trapnell et al., 2014
library(ggplot2)       # Wickham, 2016
library(org.Hs.eg.db)  # Carlson, 2019
library(Matrix)
library(igraph)
library(dplyr)
library(rlang)
```

---

## Software Compatibility Patches

Monocle 2 was developed for older versions of dplyr and igraph and is 
incompatible with current versions without modification. All Monocle scripts 
include programmatic patches applied at runtime to resolve the following 
incompatibilities:

- dplyr underscore functions (select_(), group_by_() etc.) deprecated in dplyr 0.7.0
- igraph::dfs() argument names changed in newer igraph versions
- igraph::nei() removed in igraph 2.1.0 and replaced with igraph::neighbors()

These patches are applied automatically when each script is sourced and do 
not require manual intervention.

---

## References

Camp, J.G. et al. (2015) 'Human cerebral organoids recapitulate gene expression 
programs of fetal neocortex development', Proceedings of the National Academy of 
Sciences, 112(51), pp. 15672-15677.

Trapnell, C. et al. (2014) 'The dynamics and regulators of cell fate decisions 
are revealed by pseudotemporal ordering of single cells', Nature Biotechnology, 
32(4), pp. 381-386.

Hao, Y. et al. (2023) 'Dictionary learning for integrative, multimodal and 
scalable single-cell analysis', Nature Methods, 20(4), pp. 539-541.

Wickham, H. (2016) ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag.

Carlson, M. (2019) org.Hs.eg.db: Genome wide annotation for Human. R package.

Cao, J. et al. (2019) 'The single-cell transcriptional landscape of mammalian 
organogenesis', Nature, 566(7745), pp. 496-502.

Bergen, V. et al. (2020) 'Generalizing RNA velocity to transient cell states 
through dynamical modeling', Nature Biotechnology, 38(12), pp. 1408-1414.
