"Statistical analysis and visualization of Saccharomyces cerevisiae RNA-Seq data"
# This analysis is based on the BINF 6110 - Genomic Methods Tutorials 7 & 8

# Install and load libraries ----
# if (!require("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")
# install.package("tidyverse")
# BiocManager::install("tximport")
# BiocManager::install("tximportData")
# BiocManager::install("rtracklayer")
# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")
# BiocManager::install("org.Sc.sgd.db")
# BiocManager::install("clusterProfiler")
# install.packages("pheatmap")

library(tidyverse)
library(tximport)
library(rtracklayer)
library(DESeq2)
library(apeglm)
library(pheatmap)
library(org.Sc.sgd.db)
library(clusterProfiler)

# Import data ----
# Set working directory
# Get the name of the current folder
current_dir <- basename(getwd())

# Only change directory if we aren't already in the 'scripts' folder
if (current_dir != "scripts") {
  setwd("./scripts/")
}

## Create objects for Salmon output ----
# Import sample metadata with stages of biofilm
samples <- read.csv("../data/metadata/metadata.csv")

# Make stage a factor with three levels - early, thin and mature biofilm
samples$stage <- factor(samples$stage, levels = c("early", "thin", "mature"))

# Make sra_accession the sample names 
rownames(samples) <- samples$sra_accession

# Point to directory with Salmon quant.sf files 
dir <- "../results/salmon_quant"

# Create a filepath for each sample and name each sample's filepath 
files <- file.path(dir, samples$sra_accession, "quant.sf")
names(files) <- samples$sra_accession

# Check that all files exist in the filepath location
all(file.exists(files))

## Create table for genes and their counts from Salmon ----
# Get list of names that Salmon used for quantifying and remove suffixes
all_salmon_names <- read.table(files[1], header=TRUE)$Name
clean_names <- gsub("\\.[0-9]+$", "", all_salmon_names)

# Map transcript names to gene symbols
gene_map <- bitr(clean_names, 
                 fromType = "REFSEQ", 
                 toType = c("ENTREZID", "GENENAME", "ORF"),
                 OrgDb = org.Sc.sgd.db)

# Build the final tx2gene table
tx2gene <- data.frame(TXNAME = all_salmon_names, 
                      GENEID = gene_map$ORF[match(clean_names, gene_map$REFSEQ)],
                      GENENAME = gene_map$GENENAME[match(clean_names, gene_map$REFSEQ)],
                      stringsAsFactors = FALSE)

## Import salmon data counts in
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)
head(txi$counts)

# Analyzing RNA-seq data with DESeq2 ----
# Create DESeqDataSet object
dds_txi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ stage)

# Run DESeq analysis with early as the baseline / intercept
dds_txi$stage <- factor(dds_txi$stage, levels = c("early", "thin", "mature"))
dds <- DESeq(dds_txi)

# Shrink "early to thin" and "early_to_mature" using apeglm
res_early_to_thin <- lfcShrink(dds, coef = "stage_thin_vs_early", type = "apeglm")
res_early_to_mature <- lfcShrink(dds, coef = "stage_mature_vs_early", type = "apeglm")

# Run DESeq analysis with thin as the baseline / intercept (releveling)
dds_thin <- dds
dds_thin$stage <- relevel(dds_thin$stage, ref = "thin")
dds_thin <- DESeq(dds_thin)

# Shrink "thin to mature" using apeglm
res_thin_to_mature <- lfcShrink(dds_thin, coef = "stage_mature_vs_thin", type = "apeglm")

# View coefficients
resultsNames(dds)
resultsNames(dds_thin)
samples$stage

# Plotting Differential Gene Expression Analysis Results ----

## Volcano Plot ----
# Mark genes as upregulated, downregulated, or not significant and exclude genes with under 2-fold change (log2FoldChange < 1)
# Create a function that labels results output from DESeq as significant or not
get_significance <- function(res_object, transition_label){
  df <- as.data.frame(res_object)
  df$gene <- rownames(df)
  df$significant <- ifelse(df$padj < 0.05 & abs(df$log2FoldChange) > 1, 
                           ifelse(df$log2FoldChange > 0, "Upregulated", "Downregulated"), "Not Significant")
  df$transition <- transition_label
  return(na.omit(df))
}

# Create data frames that store the results for each stage transition
df_res_early_to_thin <- get_significance(res_early_to_thin, "1. Early to Thin")
df_res_thin_to_mature <- get_significance(res_thin_to_mature, "2. Thin to Mature")
df_res_early_to_mature <- get_significance(res_early_to_mature, "3. Early to Mature")

# Create one combined data frame with all results from all transitions
df_all_res <- rbind(df_res_early_to_thin, df_res_thin_to_mature, df_res_early_to_mature)

# Ensure transition column is a factor with three levels
df_all_res$transition <- factor(df_all_res$transition)

# Plot log2foldchange against -log10(adjusted p value)
volcano_plot_transitions <- ggplot(df_all_res, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.3) +
  # Add "2-fold" and "p-value" threshold lines for what is considered relevant and significant
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  scale_color_manual(values = c("Downregulated" = "#0072B2", "Not Significant" = "#999999", "Upregulated" = "#E69F00")) +
  labs(x = "Log2 Fold Change", y = "-Log10 padj", 
       #title = "Volcano Plot: Gene Expression Changes Across Biofilm Stages of Development"
       color = NULL) +
  theme(legend.position = "right") + 
  facet_wrap(~ transition) + 
  theme(legend.position = "bottom") + 
  theme_minimal()

ggsave("../figs/03_volcano_plot_transitions.png", plot = volcano_plot_transitions, width = 10, height = 6, dpi = 600)

## Heatmap of Top 20 Genes (Early to Mature) ----
# Remove NA values
res_clean_early_to_mature <- na.omit(res_early_to_mature)

# Select top 20 genes by p values from non-NA genes
top_genes <- head(order(abs(res_clean_early_to_mature$padj), decreasing = FALSE), 20)
gene_names <- rownames(res_clean_early_to_mature)[top_genes]

# Extract transformed & normalized counts with a variance stabilizing transformation
vsd <- vst(dds)

# Store counts in a matrix for the heatmap
mat_vsd <- assay(vsd)[gene_names, ]

# Explicitly create heatmap labels to have gene names instead of ORFs if gene names exist
heatmap_labels <- tx2gene$GENENAME[match(rownames(mat_vsd), tx2gene$GENEID)]
rownames(mat_vsd) <- ifelse(is.na(heatmap_labels) | heatmap_labels == "",
                            paste0(rownames(mat_vsd), " (UNNAMED)"), 
                            heatmap_labels)

# Create an annotation data frame and assign the sample names as rownames
annotation_df <- data.frame(Stage = samples$stage)
rownames(annotation_df) <- colnames(mat_vsd)
colnames(mat_vsd) <- samples$sample_id

# Create heatmap
png("../figs/04_heatmap_top20.png", width = 10, height = 6, units = "in", res = 600)
pheatmap(mat_vsd, 
         scale = "row",
         cluster_cols = FALSE, 
         # Modify data frame to display formatted headers
         annotation_col = data.frame(
           "Biofilm Stage" = factor(samples$stage, 
                                    labels = c("Early", "Thin", "Mature")),
           row.names = colnames(mat_vsd)),
         show_colnames = TRUE,
         annotation_names_col = FALSE, 
         color = colorRampPalette(c("#0072B2", "white", "#E69F00"))(100))
dev.off()

## PCA Plot of Biofilm Stages ----
# Get the coordinates using plotPCA from DESeq2
pca_data <- plotPCA(vsd, intgroup = "stage", returnData = TRUE)

# Get percent variance explained by the top two principal components
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Display biofilm stage by shape
pca_plot_stages <- ggplot(pca_data, aes(PC1, PC2, shape = stage, color = stage)) +
  geom_point(size = 4, alpha = 0.75) + 
  scale_color_manual(name = "Biofilm Stage", 
                     values = c("early" = "#F8766D",
                                "thin" = "#00BA38",
                                "mature" = "#619CFF"),
                     labels = c("early" = "Early",
                                "thin" = "Thin",
                                "mature" = "Mature")) +
  scale_shape_manual(name = "Biofilm Stage",
    values = c("early" = 16,
               "thin" = 17,
               "mature" = 15),
    labels = c("early" = "Early",
               "thin" = "Thin",
               "mature" = "Mature")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  coord_fixed()

ggsave("../figs/05_pca_plot_stages.png", plot = pca_plot_stages, width = 10, height = 6, dpi = 600)

# Extract all Significant Genes from each stage ----
# Create a function that extracts significant genes and whether they were up- or down-regulated
get_sig_genes <- function(df, label) {
  up <- df %>% filter(significant == "Upregulated") %>% pull(gene)
  down <- df %>% filter(significant == "Downregulated") %>% pull(gene)
  
  res <- list()
  res[[paste0(label, " (Up)")]] <- up
  res[[paste0(label, " (Down)")]] <- down
  return(res)
}

# Create a list of all significant genes and their transitions
sig_list <- c(
  get_sig_genes(df_res_early_to_thin, "Early to Thin"),
  get_sig_genes(df_res_thin_to_mature, "Thin to Mature"),
  get_sig_genes(df_res_early_to_mature, "Early to Mature")
)

# Save results for functional enrichment analysis
saveRDS(sig_list, "../results/deseq2/sig_gene_lists.rds")
