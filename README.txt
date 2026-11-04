# TCGA-BRCA DESeq2 Differential Expression Analaysis

## Description
R scripts used for differential gene expression analysis of TCGA breast cancer RNA-seq data.

## Dataset used
The analysis was done using publicaly available TCGA-BRCA data.

Comparison:
 
- Luminal-A tumor subtype samples vs Normal
- Basal like tumor subtype samples vs Normal
- Basal like tumor subtype samples vs Luminal-A tumor subtype samples 

## Differential Gene Expression Analaysis
Differential gene expression analysis was conducted using DESeq2 Bioconductor package.

workflow:
- RNA-seq STAR counts were obtained for TCGA-BRCA project using the TCGAbiolinks R/Bioconductor package
- Basal-like, Luminal-A and Normal samples were selected
- Differential expression analysis was performed using DESeq2
- Differentially expressed genes for HSPGs were identified based on statistical significance p-adj < 0.05.
- Results for HSPGs for each comparison was plotted using ggplot2.

## Gene expression Correlation Analysis:

Expression correlation analysis was performed between target genes in Normal, Luminal-A and Basal-like subtype.

workflow:
- Count data were transformed using the variance stabilizing transformation (VST) presented in DESeq2
- VST-transformed expression values for our target genes were used for correlation analysis.
- Spearman correlation coefficients were calculated between target genes using "stat-cor()" function from ggpubr package.
- Correlation plots with p-values  were generated using ggplot2.

## Files
- "TCGA_BRCA_HSPGs_Expression_Analysis.R"
Differential and Correlation Expression analysis of genes of HSPGs workflow in Breast Cancer.

- "sessionInfo.txt"
R version and Packages information used for thee analysis.

## Software and Packages: 

Analysis was performed in R using Bioconductor packages:
- TCGAbiolinks
-SummarizedExperiment
-DESeq2
-ggplot2
-ggpubr

## Data Availability
TCGA_BRCA data are publicaly available through the Genomic Data Commoms (GDC) portal. 