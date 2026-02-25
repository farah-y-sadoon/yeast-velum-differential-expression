"Functional enrichment analysis of Saccharomyces cerevisiae RNA-Seq data between early and mature biofilm stages"
# This analysis is based on the BINF 6110 - Genomic Methods Tutorial 8

# Install and load libraries ----
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Sc.sgd.db")
# BiocManager::install("enrichplot")
# install.packages("tidyverse")

library(clusterProfiler)
library(org.Sc.sgd.db)
library(enrichplot)
library(tidyverse)


# Import data ----
# Set working directory
# Get the name of the current folder
current_dir <- basename(getwd())

# Only change directory if we aren't already in the 'scripts' folder
if (current_dir != "scripts") {
  setwd("./scripts/")
}

# Load significant gene lists with up- and down-regulated groups
sig_list <- readRDS("../results/deseq2/sig_gene_lists.rds")

# Convert all 6 lists from ORF to ENTREZID
sig_entrez_list <- lapply(sig_list, function(x) {
  bitr(x, fromType = "ORF", toType = "ENTREZID", OrgDb = org.Sc.sgd.db)$ENTREZID
})

# GO Enrichment: Over-representation Analysis (ORA) ----
# Conduct ORA across three stages of biofilm development and include up- and down-regulated distinctions
compare_go_results <- compareCluster(geneCluster = sig_entrez_list, fun = "enrichGO", 
                                     OrgDb = org.Sc.sgd.db, ont = "BP",
                                     pvalueCutoff = 0.05)

## Visualize Results ----
### Overall Results ----
dotplot(compare_go_results, showCategory = 5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Simplified Results ----
# Create a simplified results object to remove redundant GO terms
simplified_results <- clusterProfiler::simplify(compare_go_results, cutoff = 0.7)

# Clean up labels for plotting
levels(simplified_results@compareClusterResult$Cluster) <- 
  gsub("_", " ", levels(simplified_results@compareClusterResult$Cluster))

# Wrap the GO Term descriptions to prevent overlap
simplified_results@compareClusterResult$Description <- 
  str_wrap(simplified_results@compareClusterResult$Description, width = 100)

# Create final GO ORA plot
go_ora_plot <- dotplot(simplified_results, showCategory = 5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8)) + 
  labs(title = "A", x = "")

# KEGG Enrichment: Over-representation Analysis (ORA) ----
compare_kegg_results <- compareCluster(geneCluster = sig_list, 
                                       fun = "enrichKEGG", 
                                       organism = "sce", # Saccharomyces cerevisiae code in database = "sce"
                                       keyType = "kegg",
                                       pvalueCutoff = 0.05)

# Clean up labels for plotting
levels(compare_kegg_results@compareClusterResult$Cluster) <- 
  gsub("_", " ", levels(compare_kegg_results@compareClusterResult$Cluster))

## Visualize Results ----
kegg_ora_plot <- dotplot(compare_kegg_results, showCategory = 5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8)) + 
  labs(title = "B", x = "")

# Save Results and Figures ----
combined_ora_plot <- go_ora_plot + kegg_ora_plot
ggsave("../figs/06_ora_go_kegg_plot.png", plot = combined_ora_plot, width = 20, height = 12, dpi = 600)
go_final_results_table <- as.data.frame(simplified_results)
write.csv(go_final_results_table, "../results/enrichment/go_enrichment_results.csv", row.names = FALSE)
kegg_final_results_table <- as.data.frame(compare_kegg_results)
write.csv(kegg_final_results_table, "../results/enrichment/kegg_enrichment_results.csv", row.names = FALSE)
