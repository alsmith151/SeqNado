
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(DESeq2))

bam_dir <- normalizePath(snakemake@params$bam_dir)

# Load the data
metadata <- read_csv(snakemake@input[[2]]) %>%
  mutate(deseq2 = factor(deseq2))  # Convert to factor with 0 as reference
counts <- read.delim(snakemake@input[[1]], comment.char = "#") %>%
  rename_with(
  ~ tools::file_path_sans_ext(basename(.)),
  everything()
  ) %>%
  column_to_rownames(var = "Geneid")

# Select only numeric columns and match them to sample names
numeric_cols <- counts %>%
  dplyr::select(where(is.numeric)) %>%
  colnames()

# Match columns to sample IDs by checking if sample ID pattern is in column name
sample_cols <- sapply(metadata$sample_id, function(sid) {
  pattern <- gsub("-", ".", sid, fixed = TRUE)
  matching_col <- numeric_cols[grep(pattern, numeric_cols, fixed = TRUE)]
  if (length(matching_col) == 0) stop(paste("No matching column found for sample:", sid))
  return(matching_col[1])
})

counts <- counts %>%
  dplyr::select(all_of(sample_cols)) %>%
  setNames(metadata$sample_id) %>%
  filter(rowSums(.) > 0)


# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~deseq2)
spikein_genes <- which(rownames(counts) %in% snakemake@params$spikein_genes)
dds <- estimateSizeFactors(dds, controlGenes=spikein_genes)
counts(dds) <- counts(dds)[!rownames(counts(dds)) %in% spikein_genes, ]
dds <- DESeq(dds, quiet = T)

# Output size factors
size_factors <- colData(dds)[, "sizeFactor"]
names(size_factors) <- colData(dds)[, "sample_id"]

# Write JSON
sf_json <- toJSON(size_factors)
writeLines(sf_json, snakemake@output$normalisation_factors)

# Write TSV
sf_df <- data.frame(
  sample = names(size_factors),
  scale_factor = as.numeric(size_factors)
)
write_tsv(sf_df, snakemake@output$normalisation_table)
