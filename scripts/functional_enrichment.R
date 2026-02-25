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
sig_list <- readRDS("../results/deseq2_results/sig_gene_lists.rds")

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

ggsave("../figs/06_go_ora_plot.png", plot = go_ora_plot, width = 12, height = 10, dpi = 600)

# Save results 
go_final_results_table <- as.data.frame(simplified_results)
write.csv(go_final_results_table, "../results/deseq2_results/go_enrichment_results.csv", row.names = FALSE)
