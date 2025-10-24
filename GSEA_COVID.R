# Gene Set Enrichment Analysis (GSEA) Script
# ------------------------------------------
# This script performs Gene Set Enrichment Analysis (GSEA) using the `fgsea` package.
# It takes as input:
# 1. A CSV file containing significant gene results from our script DGEAnalysis (e.g., log2FoldChange values).
# 2. A downloaded GMT file containing gene sets (e.g., GO biological processes).
#     For example a downloaded c5.go.bp.v2023.2.Hs.symbols.gmt.txt file from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/
#
# Required Input Files:
# - `gene_results_file`: Path to the CSV file with significant gene results.
# - `GO_file`: Path to the GMT file with gene sets.
#
# Outputs:
# - A data frame of associated pathways after GSEA.
# - A plot of the top 10 up-regulated and down-regulated pathways.
#
# The script is modular and can be reused for different datasets and gene sets.
# ------------------------------------------

# Load required libraries
library(fgsea)
library(dplyr)
library(ggplot2)
library(data.table)

# Define input files (replace these with your own file paths)
gene_results_file <- "path/to/your_gene_results.csv"  # CSV file with significant gene results, these are preferred to be those prefiltered to a significant resPadj filename from DGEAnalysis script 
GO_file <- "path/to/your_gene_sets.gmt"              # GMT file with gene sets
pval_cutoff <- 0.05                                  # Adjusted p-value cutoff for significant pathways

# Load the gene results
gene_results <- read.csv(gene_results_file)

# Often data will contain rownames in a column X instead of as row.names, this adjusts if it is the case
if ("X" %in% colnames(gene_results)) {
  rownames(gene_results) <- gene_results$X
  gene_results$X <- NULL
}

# Prepare the gene list for GSEA using logfold change results, you can use other columns if needed but should change below
gene_list <- gene_results$log2FoldChange
names(gene_list) <- rownames(gene_results)
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!duplicated(names(gene_list))]  # Remove duplicates

# Define the GSEA function
GSEA <- function(gene_list, GO_file, pval) {
  # Check for duplicates and sort the gene list
  if (any(duplicated(names(gene_list)))) {
    warning("Duplicates in gene names")
    gene_list <- gene_list[!duplicated(names(gene_list))]
  }
  if (!all(order(gene_list, decreasing = TRUE) == 1:length(gene_list))) {
    warning("Gene list not sorted")
    gene_list <- sort(gene_list, decreasing = TRUE)
  }
  
  # Load the gene sets
  myGO <- fgsea::gmtPathways(GO_file)
  
  # Run GSEA
  fgRes <- fgsea::fgsea(
    pathways = myGO,
    stats = gene_list,
    minSize = 10,  # Minimum gene set size
    maxSize = 400, # Maximum gene set size
    nperm = 10000  # Number of permutations
  ) %>%
    as.data.frame() %>%
    dplyr::filter(padj < !!pval) %>%
    arrange(desc(NES))
  
  # Print the number of significant gene sets
  message(paste("Number of significant gene sets =", nrow(fgRes)))
  
  # Collapse pathways to remove redundancy
  message("Collapsing Pathways -----")
  concise_pathways <- collapsePathways(
    data.table::as.data.table(fgRes),
    pathways = myGO,
    stats = gene_list
  )
  fgRes <- fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  
  # Add enrichment direction
  fgRes$Enrichment <- ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  
  # Select the top 10 up-regulated and down-regulated pathways
  filtRes <- rbind(
    head(fgRes, n = 10),
    tail(fgRes, n = 10)
  )
  
  # Create a summary plot
  total_up <- sum(fgRes$Enrichment == "Up-regulated")
  total_down <- sum(fgRes$Enrichment == "Down-regulated")
  header <- paste0("Top 10 (Total pathways: Up=", total_up, ", Down=", total_down, ")")
  colos <- setNames(c("firebrick2", "dodgerblue2"), c("Up-regulated", "Down-regulated"))
  
  g1 <- ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point(aes(fill = Enrichment, size = size), shape = 21) +
    scale_fill_manual(values = colos) +
    scale_size_continuous(range = c(2, 10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(
      x = "Pathway",
      y = "Normalized Enrichment Score",
      title = header
    ) +
    theme_classic()
  
  # Return results and plot
  output <- list("Results" = fgRes, "Plot" = g1)
  return(output)
}

# Run GSEA
suppressWarnings({
  gsea_results <- GSEA(gene_list, GO_file, pval_cutoff)
})

# View the results
print(gsea_results$Results)

# Save the plot
ggsave("GSEA_Top_Pathways.png", gsea_results$Plot, width = 10, height = 6)

# Save the results to a CSV file
write.csv(gsea_results$Results, "GSEA_Results.csv", row.names = FALSE)
