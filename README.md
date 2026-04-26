Overview of full post 2016 analysis pipeline: 

This repository implements an RNA-seq analysis pipeline using a modern (post-2016) workflow, aiming to process large SRR data from GEO and perform clustering and visualisation of gene expression patterns. 


The workflow consists of 5 main stages: 

Pre-processing and Quality Control 

Alignment of raw FASTQ reads using STAR 

Gene quantification using STAR (GeneCounts) 

Construction of a gene-by-sample count matrix 

Downstream clustering and visualisation using Seurat 

This stage retrieves raw sequencing data from the Sequence Read Archive (SRA) using SRR accession numbers obtained from GEO. The accession list is used as the input for automated batch downloading and conversion of sequencing files. 

The workflow uses the SRA Toolkit to download .sra files and converted them into paired-end FASTQ files. Prefetch was used to retrieve the raw SRR data, while fasterq-dump converted each SRR file into readable FASTQ format. This produced paired FASTQ files for each sample, which were compressed (Fastqz) 

Due to the large number of samples in the post-2016 dataset, this step was executed on ALICE HPC using job arrays. This allowed multiple SRR accessions to be processed in parallel while staying within HPC job limits (16). Temporary and output directories were assigned within the scratch file system for trouble shooting and avoiding overwriting.  

After conversion, the generated FASTQ files formed the main input for STAR alignment. This stage was essential because STAR requires FASTQ files rather than SRR archive files as input. 

 
(STAR) – STAR_align.sh 

This script performs the alignment step of the pipeline, mapping raw FASTQ reads to the human reference genome (GRCh38) to generate sorted BAM files. The workflow uses STAR, a highly efficient splice-aware aligner. 

STAR executed on a SLURM-managed HPC cluster (ALICE), enabling parallel processing of a large dataset (~730 samples). A pre-generated genome index based on the GRCh38 reference genome and corresponding GTF annotation is used to improve alignment speed and accuracy. 

Alignment is performed in batch mode using SLURM job arrays, allowing multiple samples to be processed simultaneously while respecting HPC resource limits. The output for each sample includes sorted BAM files and alignment logs. This stage is computationally intensive but significantly faster than legacy aligners, with the full dataset processed efficiently through parallelisation. 

 
Gene Quantification (STAR Gene Counts) 

This step performs gene-level quantification directly during alignment using STAR’s built-in quantMode GeneCounts functionality. This approach eliminates the need for external tools and improves overall efficiency. For each sample, STAR produced a ReadsPerGene.out.tab file containing raw read counts assigned to each gene. 

Count Matrix Construction 

This step aggregates gene count data from individual samples into a single count matrix for downstream analysis. All ReadsPerGene.out.tab files are merged into a unified matrix (post-2016_matrix), where rows represent genes and columns represent samples. This process ensures consistent gene ordering across samples and handles missing values appropriately. 

The resulting matrix serves as the primary input for downstream analysis in R. 

 
Downstream Analysis & Visualisation (Seurat) 

This stage performs clustering and visualisation using Seurat (v5) in R. The count matrix is loaded into Seurat to create a Seurat object, followed by standard preprocessing steps including normalisation, feature selection, and scaling. Dimensionality reduction is performed using Principal Component Analysis (PCA), followed by UMAP for visualisation of sample neighbouring relationships in reduced-dimensional space. Clustering is applied to identify groups of samples with similar gene expression patterns. Marker genes are analysed to interpret cluster identities and assess biological relevance. Although the dataset lacks true single-cell resolution, this approach enables exploration of transcriptional structure and comparison of expression patterns across samples. The outputs include UMAP plots and cluster assignments, providing a visual representation of the dataset structure and enabling comparison with expected biological trends. 

 
Notes and trouble shooting:
This pipeline was developed as part of an academic project to implement and evaluate a modern RNA-seq workflow. It highlights key differences between the original and post-2016 approaches, particularly in alignment efficiency.  
Some SRR accessions produced incomplete outputs, indicating silent processing failures. These were identified through log inspection and resolved by re-running affected samples, ensuring completeness of the dataset prior to downstream analysis. 
