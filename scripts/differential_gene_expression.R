'''
Statistical analysis and visualization of Saccharomyces cerevisiae RNA-Seq data 
'''
# Install and load libraries ----
# if (!require("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")
# BiocManager::install("tximport")
# BiocManager::install("tximportData")
# BiocManager::install("rtracklayer")
# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")

library(tidyverse)
library(tximport)
library(rtracklayer)
library(DESeq2)
library(apeglm)

# Import data ----

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
# Get list of names that Salmon used for quantifying
all_salmon_names <- read.table(files[1], header=TRUE)$Name

# Create a mapping table by taking the gene ID from the names used by Salmon - looks for "C4S56_" followed by digits (gene ID)
gene_ids <- regmatches(all_salmon_names, regexpr("C4S56_[0-9]+", all_salmon_names))

# Build the final tx2gene table
tx2gene <- data.frame(TXNAME = all_salmon_names,
                      GENEID = gene_ids,
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

# CONTINUE FROM HERE: SAVE RESULTS
# VISUALIZE WITH PLOTS: 
  # HEATMAP OF DIFFERENTIALLY EXPRESSED GENES
  # PCA PLOT
  # VOLCANO PLOT OR SOME OTHER INTERESTING VISUALIZATION
