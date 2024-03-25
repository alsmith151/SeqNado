library(edgeR)
library(tidyverse)

# Load the data
counts <- read_table(snakemake@input[[1]])
metadata <- read_table(snakemake@input[[2]])

# Create a DGEList object
dge_list <- DGEList(counts=counts, group=metadata)
dge_list <- calcNormFactors(dge_list)

# Save the scaling factors
library_sizes <- dge_list$samples

# Save the scaling factors
write_tsv(library_sizes, snakemake@output[[1]])

