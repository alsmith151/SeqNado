
library(tidyverse)

# Load the data
counts <- read_table(snakemake@input[[1]])
metadata <- read_table(snakemake@input[[2]])

spikein_rows <- which(rownames(counts) %in% c("AmpR_seq", "Cas9_5p_seq", "Cas9_3p_seq"))

spikein_counts <- counts[spikein_rows, ]
size_factors <- colSums(spikein_counts) / median(colSums(spikein_counts))
normalized_counts <- sweep(counts, 2, size_factors, FUN="/")

size_factors <- colData(dds)[, "sizeFactor"]
  names(size_factors) <- colData(dds)[, "sample"]
  sf <- toJSON(size_factors)
  writeLines(sf, "seqnado_output/resources/all_normalisation_factors.json")