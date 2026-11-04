

#LOADING THE DATASET TCGA BRCA AND ADDING COUNTS DATA AND CLINICAL DATA
library(TCGAbiolinks)
library(SummarizedExperiment)

query_BRCA <- GDCquery( project = "TCGA-BRCA", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")

options(timeout = 1800)
GDCdownload(query_BRCA, method = "api", files.per.chunk = 50)

brca_data <- GDCprepare(query_BRCA)

clinical_data <- colData(brca_data)
clinical_data_df <- as.data.frame(clinical_data)
counts_matrix <- SummarizedExperiment :: assay(brca_data, "unstranded")
counts_df <- as.data.frame(counts_matrix)

#TO DISCOVER NAMES AND NUMBERS OF SAMPLES PRESENT IN THE DATA SET 
table(clinical_data_df$paper_BRCA_Subtype_PAM50)
table(clinical_data_df$sample_type)

#creating column for DESEQ CONTAINING TUMOR-BASAL, TUMOR-LUMA AND NORMAL SUBTYPES
clinical_data_df <- clinical_data_df %>% mutate(AnalysisGroup = case_when(sample_type =="Solid Tissue Normal" ~ "Healthy_Control",paper_BRCA_Subtype_PAM50 == "Basal" & sample_type == "Primary Tumor" ~ "Basal_Tumor" , paper_BRCA_Subtype_PAM50 == "LumA" & sample_type == "Primary Tumor" ~ "LumA_Tumor", TRUE ~ NA_character_))



#Ceating counts and clinical subsets to be as input for dds (DGE) 
analysis_clinical <- clinical_data_df %>% filter(!is.na(AnalysisGroup))%>% mutate(AnalysisGroup= factor(AnalysisGroup, levels = c("Healthy_Control", "Basal_Tumor", "LumA_Tumor")))
target_samples <- rownames(analysis_clinical)
final_counts <- counts_df[target_samples]
counts_with_names <- as.data.frame(final_counts)
all(colnames(final_counts) ==rownames(analysis_clinical))



#Creating DESEQ object

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = final_counts, colData = analysis_clinical, design = ~ AnalysisGroup)
dds <- DESeq(dds)


#Extracting results for each condition 
gene_metadata <- as.data.frame(SummarizedExperiment:: rowData(brca_data))
res_Normal_vs_Basal <- results(dds, contrast = c("AnalysisGroup", "Basal_Tumor", "Healthy_Control"))
res_Normal_vs_LumA <- results(dds, contrast = c("AnalysisGroup", "LumA_Tumor", "Healthy_Control"))
res_Basal_vs_LumA <- results(dds, contrast = c("AnalysisGroup", "Basal_Tumor", "LumA_Tumor"))
res_Normal_vs_LumA_df <- as.data.frame (res_Normal_vs_LumA)
res_Normal_vs_Basal_df <- as.data.frame (res_Normal_vs_Basal)                             
res_Basal_vs_LumA_df <- as.data.frame (res_Basal_vs_LumA)                            

                             
#ADDING GENE SYMBOLS FOR OUR RESULTS
res_Normal_vs_LumA_df <- as.data.frame(res_Normal_vs_LumA_df) %>% tibble::rownames_to_column(var = "Ensembl_ID")
res_Normal_vs_LumA_df$gene_symbol <- gene_metadata$gene_name[match(res_Normal_vs_LumA_df$Ensembl_ID,gene_metadata$gene_id)]
res_Normal_vs_Basal_df <- res_Normal_vs_Basal_df %>% tibble::rownames_to_column(var = "Esembl_ID")
res_Normal_vs_Basal_df$gene_symbol <- gene_metadata$gene_name[match(res_Normal_vs_Basal_df$Esembl_ID, gene_metadata$gene_id)]
res_Basal_vs_LumA_df <- as.data.frame(res_Basal_vs_LumA) %>% tibble::rownames_to_column(var = "Esembl_ID")
res_Basal_vs_LumA_df$gene_symbol <- gene_metadata$gene_name [match(res_Basal_vs_LumA_df$Esembl_ID, gene_metadata$gene_id)]







#correlation between EXT1 and SDC genes
library(ggplot2)
library(ggpubr)
install.packages(c("tidyr", "ggpubr"))
library(tidyr)
library(dplyr)
library(DESeq2)
#NORMALIZING THE COUNTS
vsd <- vst(dds, blind = FALSE)

vst_matrix <- assay(vsd)
# PCA Plot 1: Standard DESeq2 PCA by Subtype (Healthy_Control, Basal_Tumor, LumA_Tumor)
plotPCA(vsd, intgroup = "AnalysisGroup") +
  scale_color_manual(values = c(
    "Healthy_Control" = "#F8766D",
    "Basal_Tumor"     = "#00BA38",
    "LumA_Tumor"      = "#619CFF"
  )) +
  theme_bw()

# PCA Plot 2: PCA grouped by Tissue Site / Center Code
# (Extracting tissue site prefix from TCGA barcodes, e.g., '3C', 'A1', 'A2', etc.)
vsd$tissue_site <- factor(substr(colnames(vsd), 6, 7))
plotPCA(vsd, intgroup = "tissue_site") +
  theme_bw()

# PCA Plot 3: PCA grouped by Sample Type & Tissue Preservation Condition
# (Categorizing sample preservation: FFPE vs OCT / Normal vs Tumor)
vsd$sample_preservation <- paste0(
  ifelse(grepl("FFPE", vsd$sample_type, ignore_case = TRUE), "FFPE",
         ifelse(grepl("OCT", vsd$sample_type, ignore_case = TRUE), "OCT", "Unknown")),
  ":",
  ifelse(vsd$AnalysisGroup == "Healthy_Control", "Normal", "Tumor")
)
plotPCA(vsd, intgroup = "sample_preservation") +
  theme_bw()
mapping_table <- res_Basal_vs_LumA_df %>% filter(gene_symbol %in% target_genes_coexpression) %>% select( Esembl_ID, gene_symbol)

target_ensg <- mapping_table$Esembl_ID

subset_vst <- vst_matrix [target_ensg, , drop = FALSE]

rownames(subset_vst) <- mapping_table$gene_symbol[match(rownames(subset_vst), mapping_table$Esembl_ID)]
exp_df <- as.data.frame (t(subset_vst))
exp_df$barcode <- rownames(exp_df)
plot_df <- inner_join(exp_df, analysis_clinical, by = "barcode")
analysis_clinical$SampleID <- rownames(analysis_clinical)
table(plot_df$AnalysisGroup)

plot_long <- plot_df %>% pivot_longer(cols = c( "SDC1","SDC2","SDC3","SDC4"), names_to = "SDC_Gene", values_to = "SDC_Expression")
plot_long$AnalysisGroup <- factor(plot_long$AnalysisGroup, levels = c("Healthy_Control", "LumA_Tumor", "Basal_Tumor"))

#SPEARMAN CORRELATION BETWEEN EXT1 AND SDC IN NORMAL, LUMA AND BASAL SUBTYPES
library(dplyr)
library(tidyr)
library(ggplot2)

#CORRELATION DATA FOR TARGET GENES
cor_input_df <- plot_df %>%
  select(barcode, AnalysisGroup, EXT1, SDC1, SDC2, SDC3, SDC4)

#CALCULATION OF EACH CORRELATION IN EVERY SUBTYPE WITH P VALUES
subtypes <- unique(cor_input_df$AnalysisGroup)
sdc_genes <- c("SDC1", "SDC2", "SDC3", "SDC4")
cor_results_list <- list()

for (grp in subtypes) {
  sub_data <- cor_input_df %>% filter(AnalysisGroup == grp)
  
  for (gene in sdc_genes) {
    test <- cor.test(sub_data$EXT1, sub_data[[gene]], method = "spearman", exact = FALSE)
    
    cor_results_list[[length(cor_results_list) + 1]] <- data.frame(
      Subtype = grp,
      Gene1 = "EXT1",
      Gene2 = gene,
      Correlation = test$estimate,
      P_value = test$p.value,
      N = sum(!is.na(sub_data$EXT1) & !is.na(sub_data[[gene]]))
    )
  }
}

correlation_results <- bind_rows(cor_results_list)

#COMPARISON OF CORRELATIONS BETWEEN EVERY SUBTYPE USING FISHER'S Z-TEST

compare_correlations <- function(cor_df) {
  wide_cor <- cor_df %>%
    select(Gene1, Gene2, Subtype, Correlation, N) %>%
    pivot_wider(names_from = Subtype, values_from = c(Correlation, N))
  
  comparisons <- list(
    c("LumA_Tumor", "Healthy_Control"), # LumA - Normal
    c("Basal_Tumor", "LumA_Tumor"), # Basal - LumA
    c("Basal_Tumor", "Healthy_Control") # Basal - Normal
  )
  
  diff_results <- list()
  
  for (pair in comparisons) {
    gA <- pair[1]
    gB <- pair[2]
    
    col_a <- paste0("Correlation_", gA)
    col_b <- paste0("Correlation_", gB)
    
    if (col_a %in% colnames(wide_cor) & col_b %in% colnames(wide_cor)) {
      res <- wide_cor %>% rowwise() %>% mutate(
        r1 = .data[[col_a]],
        n1 = .data[[paste0("N_", gA)]],
        r2 = .data[[col_b]],
        n2 = .data[[paste0("N_", gB)]],
        
        # Fisher z-transformation
        z1 = 0.5 * log((1 + r1) / (1 - r1)),
        z2 = 0.5 * log((1 + r2) / (1 - r2)),
        
        # Standard error and Z-score for difference (n - 3 correction)
        se_diff = sqrt(1 / (n1 - 3) + 1 / (n2 - 3)),
        z_score = (z1 - z2) / se_diff,
        p_value_diff = 2 * (1 - pnorm(abs(z_score))),
        
        Comparison = paste0(gA, " vs ", gB)
      ) %>% ungroup() %>% select(Gene1, Gene2, Comparison, r1, r2, z_score, p_value_diff)
      
      diff_results[[paste0(gA, "_vs_", gB)]] <- res
    }
  }
  
  bind_rows(diff_results)
}

diff_table <- compare_correlations(correlation_results)

#Generating plots to show the difference in correlations across the subtypes

#Annotations preparation for plots:
diff_table <- diff_table %>%
  mutate(
    Pair_Label = paste(Gene1, "-", Gene2),
    Sig_Stars = case_when(
      p_value_diff < 0.001 ~ "***",
      p_value_diff >= 0.001 & p_value_diff < 0.01 ~ "**",
      p_value_diff >= 0.01 & p_value_diff <= 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    Color_Group = case_when(
      p_value_diff < 0.05 & z_score > 0 ~ "Positive Significant",
      p_value_diff < 0.05 & z_score < 0 ~ "Negative Significant",
      TRUE ~ "Non-Significant"
    ),
    label_y = ifelse(z_score >= 0, z_score + 0.35, z_score - 0.45)
  )

custom_colors <- c(
  "Positive Significant" = "#2b8cbe", # Blue
  "Negative Significant" = "#e31a1c", # Red
  "Non-Significant" = "black" # Black
)




# Plot 1: LumA Tumor vs Healthy Control
plot_lum_norm <- diff_table %>%
  filter(Comparison == "LumA_Tumor vs Healthy_Control") %>%
  ggplot(aes(x = reorder(Pair_Label, z_score), y = z_score, fill = Color_Group)) +
  geom_col(width = 0.55, show.legend = TRUE) +
  geom_text(aes(y = label_y, label = Sig_Stars), size = 5, fontface = "bold", color = "black") +
  coord_flip() +
  scale_fill_manual(values = custom_colors) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  labs(title = "LumA Tumor vs Healthy Control",
       subtitle = "Change in EXT1-Syndecan Correlation",
       x = "Gene Pairs", y = "Z-score of Correlation Difference", fill = "Significance")

# Plot 2: Basal Tumor vs LumA Tumor
plot_basal_lum <- diff_table %>%
  filter(Comparison == "Basal_Tumor vs LumA_Tumor") %>%
  ggplot(aes(x = reorder(Pair_Label, z_score), y = z_score, fill = Color_Group)) +
  geom_col(width = 0.55, show.legend = TRUE) +
  geom_text(aes(y = label_y, label = Sig_Stars), size = 5, fontface = "bold", color = "black") +
  coord_flip() +
  scale_fill_manual(values = custom_colors) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  labs(title = "Basal Tumor vs LumA Tumor",
       subtitle = "Change in EXT1-Syndecan Correlation",
       x = "Gene Pairs", y = "Z-score of Correlation Difference", fill = "Significance")

# Plot 3: Basal Tumor vs Healthy Control
plot_basal_norm <- diff_table %>%
  filter(Comparison == "Basal_Tumor vs Healthy_Control") %>%
  ggplot(aes(x = reorder(Pair_Label, z_score), y = z_score, fill = Color_Group)) +
  geom_col(width = 0.55, show.legend = TRUE) +
  geom_text(aes(y = label_y, label = Sig_Stars), size = 5, fontface = "bold", color = "black") +
  coord_flip() +
  scale_fill_manual(values = custom_colors) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  labs(title = "Basal Tumor vs Healthy Control",
       subtitle = "Change in EXT1-Syndecan Correlation",
       x = "Gene Pairs", y = "Z-score of Correlation Difference", fill = "Significance")


print(plot_lum_norm)
print(plot_basal_lum)
print(plot_basal_norm)
target_genes_coexpression <- c("EXT1", "SDC1","SDC2", "SDC3", "SDC4")


#PLOT TO SHOW CORRELATIONS BETWEEN EXT1 AND SDC GENES IN EACH SUBTYPE
library(dplyr)
library(tidyr)
library(ggplot2)


cor_results <- plot_long %>%
  group_by(SDC_Gene, AnalysisGroup) %>%
  summarise(
    cor_test = list(cor.test(EXT1, SDC_Expression, method = "spearman", exact = FALSE)),
    .groups = "drop"
  ) %>%
  mutate(
    rho = sapply(cor_test, function(x) x$estimate),
    p_val = sapply(cor_test, function(x) x$p.value),
    # Format label text
    label = sprintf("r = %.3f\np = %.3g", rho, p_val),
    # Assign color: Red if positive & p < 0.05, Blue if negative & p < 0.05, Black otherwise
    text_color = case_when(
      p_val < 0.05 & rho > 0 ~ "#e41a1c", # Red
      p_val < 0.05 & rho < 0 ~ "#377eb8", # Blue
      TRUE ~ "black" # Black
    )
  )

ggplot(plot_long, aes(x = EXT1, y = SDC_Expression, color = AnalysisGroup)) +
 
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", color = "black", se = FALSE, linewidth = 0.8) +
  

  geom_text(
    data = cor_results,
    aes(x = -Inf, y = Inf, label = label, color = NULL, color = text_color),
    # Note: To force manual text colors per facet in geom_text, we pass the identity vector:
    colour = cor_results$text_color,
    hjust = -0.1,
    vjust = 1.2,
    size = 5.5, # Much bigger and bolded font size
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  
  facet_grid(SDC_Gene ~ AnalysisGroup, scales = "free") +
  
  scale_color_manual(values = c(
    "Healthy_Control" = "#2ca25f",
    "LumA_Tumor" = "#2b8cbe",
    "Basal_Tumor" = "#dd1c77"
  )) +
  
  theme_bw(base_size = 14) +
  labs(
    title = "Co-expression Analysis: EXT1 vs SYNDECANS in Different Breast Cancer Subtypes",
    subtitle = "Spearman Rank Correlation with Significance Highlighting (VST Normalized)",
    x = "EXT1 Expression",
    y = "Syndecan Expression"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    panel.spacing = unit(1, "lines")
  )








citation("DESeq2")
#Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)

citation("ggplot2")
# H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

citation("dplyr")
#Wickham H, François R, Henry L, Müller K, Vaughan D (2026). _dplyr: A Grammar of Data Manipulation_. doi:10.32614/CRAN.package.dplyr https://doi.org/10.32614/CRAN.package.dplyr>, R package version 1.2.0, <https://CRAN.R-project.org/package=dplyr>.

citation("TCGAbiolinks")
#Antonio Colaprico, Tiago Chedraoui Silva, Catharina Olsen, Luciano Garofano, Claudia Cava, Davide Garolini, Thais Sabedot, Tathiane Malta, Stefano M. Pagnotta, Isabella
#Castiglioni,Michele Ceccarelli, Gianluca Bontempi Houtan Noushmehr. TCGAbiolinks: An R/Bioconductor package for integrative analysis of TCGA data Nucleic Acids Research (05 May 2016) 44 (8): e71. (doi:10.1093/nar/gkv1507)

#Silva, Tiago C., et al. "TCGA Workflow: Analyze cancer genomics and epigenomics data using Bioconductor packages." F1000Research 5 (2016).

#Mounir, Mohamed, et al. "New functionalities in the TCGAbiolinks package for the study and integration of cancer data from GDC and GTEx." PLoS computational biology 15.3 (2019): e1006701.

citation("ggpubr")
#Kassambara A (2026). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.6.3.999. https://CRAN.R-project.org/package=ggpubr

citation("SummarizedExperiment")
#Morgan M, Obenchain V, Hester J, Pagès H (2025). _SummarizedExperiment: A container (S4 class) for matrix-like assays_. R package version 1.40.0,<https://bioconductor.org/packages/SummarizedExperiment>.



