# SARS-COV-2-Extracellular-Vesicles-RNAseq
Determining differential gene expression, correlations, and associations of RNA levels of three immune cell type's extracellular vesicles extracted through a microfluidics device

#FilteringMedian
Function to iteratively test different filtering conditions on a dataset. 
Specifically, it aims to identify the set of conditions that minimizes the standard deviation of median gene counts across samples, thereby reducing variability. 
At the same time, the function ensures that within each defined subgroup (based on sample characteristics provided in the metadata table), the median of the medians of the filtered raw counts remains greater than zero.

