
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(DESeq2))

# Load the data
metadata <- read_csv(snakemake@input[[2]])
counts <- read.delim(snakemake@input[[1]], comment.char = "#") %>%
  rename_with(~ gsub("(seqnado_output.aligned.)(.+)(.bam)", "\\2", .), everything()) %>%
  column_to_rownames(var = "Geneid") %>%
  dplyr::select(one_of(make.names(metadata$sample_name))) %>%
  setNames(metadata$sample_name) %>%
  filter(rowSums(.) > 0)


# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~deseq2)
spikein_genes <- which(rownames(counts) %in% snakemake@params$spikein_genes)
dds <- estimateSizeFactors(dds, controlGenes=spikein_genes)
counts(dds) <- counts(dds)[!rownames(counts(dds)) %in% spikein_genes, ]
dds <- DESeq(dds, quiet = T)

# Output size factors
size_factors <- colData(dds)[, "sizeFactor"]
names(size_factors) <- colData(dds)[, "sample_name"]
sf <- toJSON(size_factors)
writeLines(sf, snakemake@output[[1]])
