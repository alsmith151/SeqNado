suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(edgeR))


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

# Create DGEList
dge <- DGEList(counts = counts)

# Identify spike-in genes
spikein <- rownames(dge) %in% snakemake@params$spikein_genes

# Estimate normalization factors using spike-ins only
dge <- calcNormFactors(dge, method = "TMM", subset = spikein)

# Remove spike-ins
dge <- dge[!spikein, , keep.lib.sizes = FALSE]

# edgeR model fitting
design <- model.matrix(~ deseq2, data = metadata)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

# Extract normalization factors
norm_factors <- dge$samples$norm.factors
names(norm_factors) <- metadata$sample_id

# Write JSON
sf_json <- toJSON(norm_factors)
writeLines(sf_json, snakemake@output$normalisation_factors)

# Write TSV
nf_df <- data.frame(
  sample = names(norm_factors),
  scale_factor = as.numeric(norm_factors)
)
write_tsv(nf_df, snakemake@output$normalisation_table)