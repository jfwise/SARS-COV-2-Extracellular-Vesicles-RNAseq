#Function to iteratively test different filtering conditions on a dataset. 
#Specifically, it aims to identify the set of conditions that minimizes the standard deviation of median gene counts across samples, thereby reducing variability. 
#At the same time, the function ensures that within each defined subgroup (based on sample characteristics provided in the metadata table), the median of the medians of the filtered raw counts remains greater than zero.
#Requires a counts_table from gene expression
#Requires sample_info a metadata dataframe must include a Sample_ID column that has sample ids that match the counts_table columns.

#Function to calculate the median of medians across subgroups
#Required are the counts table and the sample data table
#May define subtype or another column division of samples with Type

calculate_median_of_medians <- function(
    counts_table, 
    sample_info, 
    Type) {
  # Ensure the column exists in the metadata
  if (!Type %in% colnames(sample_info)) {
    stop(paste("Column", Type, "not found in sample_info"))
  }
  
  # Dynamically reference the specified column
  medians <- apply(counts_table, 2, median)
  group_medians <- tapply(medians, sample_info[[Type]], median)
  
  return(group_medians)
}

# Placeholder for the best conditions
best_conditions <- NULL
min_stddev <- Inf

# Iterate through the number of samples allowed to have zero counts (from 0 to total number of samples)
for (samples_with_zero in 0:ncol(counts_table)) {
  
  # print(paste("samples_with_zero =", samples_with_zero))
  
  # Filter out rows based on the number of zeros in each row
  ##For each gene, this line counts the # of samples w 0 counts
  zero_counts <- rowSums(counts_table == 0)
  ##The filtered counts table include genes that meet criteria
  filtered_counts <- counts_table[zero_counts <= samples_with_zero, ]
  
  # Check if the filtered data is non-empty
  ##https://www.geeksforgeeks.org/r-next-statement/
  if (nrow(filtered_counts) == 0) next
  
  # Match sample_info with filtered_counts columns
  sample_info_subset <- sample_info[match(
    colnames(filtered_counts), 
    sample_info$Sample_ID), 
    ]
  
  # Calculate the median of the medians across subgroups
  group_medians <- calculate_median_of_medians(filtered_counts, sample_info_subset)
  
  # print(paste("group_medians =", group_medians))
  
  # Ensure the median of medians is more than zero for each group
  if (all(group_medians > 0)) {
    
    # print(paste("all(group_medians > 0"))
    
    # Calculate the standard deviation of the medians
    sample_medians <- apply(filtered_counts, 2, median)
    stddev_medians <- sd(sample_medians)
    
    # print(paste("stddev_medians =", stddev_medians))


    # Check if this is the best set of conditions found so far
    if (stddev_medians < min_stddev) {
      
      # print(paste("stddev_medians < min_stddev"))


      min_stddev <- stddev_medians
      best_conditions <- list(samples_with_zero = samples_with_zero, filtered_counts = filtered_counts)
    }
  }
}


# Display the optimal conditions found
if (!is.null(best_conditions)) {
  print("OPTIMAL CONDITIONS FOUND")
  print(paste("The number of samples allowed to have zero counts for a particular gene is:", best_conditions$samples_with_zero))
  print(paste("The associated standard deviation of the median filtered gene counts across samples was:",min_stddev))
  print(paste("For the three Type subgroups, the associated medians of median filtered gene counts across samples were:"))
  print(calculate_median_of_medians(best_conditions$filtered_counts,sample_info_subset))
  } else {
  print("No suitable conditions found.")
  }

counts_filtered_medians <- best_conditions$filtered_counts

return(counts_filtered_medians)
