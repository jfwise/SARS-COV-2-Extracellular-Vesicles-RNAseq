# Function for Differential Gene Expression (DGE) on datasets prefiltered on median of medians
# Specificallt, it performs Differential Gene Expression (DGE) analysis using DESeq2
# Requires
# Requires counts_filtered_medians a matrix of pre filtered counts table generated from FilteringMedian.R 
# Requires sample_info a metadata dataframe must include a Sample_ID column that has sample ids that match the counts_table columns.
# Requires Type a subtype or another column in sample_info for division 
# Requires comparison_variable the variable of interest for comparison commonly supplied as a column in sample_info
# Requires the output_dir directory where results and plots will be saved
# Optional is adjustment to give a variable to control the results for, here we use the example pool (referencing our multiple sequencing runs)
# The script outputs the DGE results and the filtered results (adjusted p-value < 0.05) as a CSV file. along wit a volcano plot highlighting significant genes.


# Load required libraries
library(data.table)
library(dplyr)
library(DESeq2)
library(EnhancedVolcano)

# Function to perform Differential Gene Expression (DGE) analysis
perform_dge_analysis <- function(counts_filtered_medians, sample_info, Type, comparison_variable, output_dir, adjustment = NULL) {
  
  # Set the working directory for saving results
  setwd(output_dir)
  
  # Filter metadata for the specified EV_Type
  sample_subset <- filter(sample_info, sample_metadata$Type == Type)
  
  # Match column order in counts table to sample metadata
  counts_filtered_medians <- setcolorder(counts_filtered_medians, sample_info$Sample_ID)
  
  # Subset counts table for the selected samples
  counts_subset <- counts_filtered_medians[, sample_subset$Sample_ID]
  
  # Create DESeq2 dataset
# Create the design formula
  if (!is.null(control_variable)) {
    design_formula <- as.formula(paste("~", control_variable, "+", comparison_variable))
  } else {
    design_formula <- as.formula(paste("~", comparison_variable))
  }
 dds <- DESeqDataSetFromMatrix(
    countData = counts_subset, 
    colData = sample_subset, 
    design = design_formula
  )

# Perform DGE analysis
  dds <- DESeq(dds, parallel = FALSE)
  
  # Calculate normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  
  # Extract results
  dds_results <- results(dds, alpha = 0.05)
  summary(dds_results)
  
  # Save full DGE results
  results_df <- as.data.frame(dds_results)
  filename <- paste(Sys.Date(), Type, comparison_variable, "dds_res_all", "alpha", 0.05, ".csv", sep = "_")
  write.csv(results_df, file = filename, row.names = TRUE)
  
  # Filter results by adjusted p-value (padj)
  dds_results <- dds_results[!is.na(dds_results$padj), ]
  dds_results_ordered <- dds_results[order(dds_results$padj), ]
  res_padj <- dds_results_ordered[dds_results_ordered$padj < 0.05, ]
  
  # Save filtered results
  filename_padj <- paste(Sys.Date(), Type, comparison_variable, "resPadj", "alpha", 0.05, ".csv", sep = "_")
  write.csv(res_padj, file = filename_padj, row.names = TRUE)

 # Generate volcano plot
  generate_volcano_plot(dds, res_padj, Type, comparison_variable)
}

# Function to generate a volcano plot
generate_volcano_plot <- function(dds, res_padj, Type, comparison_variable) {
  # Extract results for volcano plot
  res_volcano <- results(dds, alpha = 0.05)
  
  # Determine cutoffs for volcano plot
  volcano_pval_cutoff <- volcano_pval_cutoff <- max(resSig@listData[["pvalue"]])
  volcano_abs_lfc_cutoff <- min(abs(resSig@listData[["log2FoldChange"]]))
  
  # Define custom key-value pairs for coloring
  keyvals <- ifelse(
    res_volcano$padj < 0.05 & res_volcano$log2FoldChange >= volcano_abs_lfc_cutoff, 'red2',
    ifelse(
      res_volcano$padj < 0.05 & res_volcano$log2FoldChange <= -volcano_abs_lfc_cutoff, 'royalblue',
      'grey30'
    )
  )

keyvals[is.na(keyvals)] <- 'grey30'
  names(keyvals)[keyvals == 'red2'] <- 'upreg'
  names(keyvals)[keyvals == 'royalblue'] <- 'downreg'
  names(keyvals)[keyvals == 'grey30'] <- 'nonsig'
  
  # Generate and save volcano plot
  filename_volcano <- paste(Sys.Date(), Type, comparison_variable, "EnhancedVolcano", "alpha", 0.05, ".pdf", sep = "_")
  pdf(filename_volcano, height = 11, width = 8.5)
  
  EnhancedVolcano(
    res_volcano,
    title = paste("DESeq2 -", Type, comparison_variable),
    subtitle = paste("padj < 0.05; abs(log2FC) >=", volcano_abs_lfc_cutoff),
    lab = rownames(res_volcano),
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = rownames(res_volcano)[which(names(keyvals) %in% c('upreg', 'downreg'))],
    pCutoff = volcano_pval_cutoff,
    FCcutoff = volcano_abs_lfc_cutoff,
    colCustom = keyvals,
    colAlpha = 1,
    labSize = 1.5
  )
  
  dev.off()
}
