# Pairwise metabolic transformation analysis: R-side input preparation
# Default workflow: raw RNA-seq counts, source-vs-target DESeq2 contrast,
# target-state VST discretization, and GEM-compatible exports for MATLAB/rMTA.

library(DESeq2)
library(dplyr)
library(tibble)

# -----------------------------------------------------------------------------
# User-defined parameters
# -----------------------------------------------------------------------------
expression_file <- "data/expression_counts.tsv"
metadata_file   <- "data/metadata.csv"
gem_gene_file   <- "data/gem_genes.csv"

use_mrna_filter  <- TRUE
mrna_gene_file    <- "data/mrna_genes.csv"

use_id_conversion <- TRUE
id_mapping_file   <- "data/id_mapping.csv"

source_label <- "source"
target_label <- "target"

pval_cutoff <- 0.05
sd_multiplier <- 0.6
consensus_threshold <- 0.8

out_dir <- "results/R_outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
discretize_sample <- function(x, sd_multiplier = 0.6) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  ifelse(x <= m - sd_multiplier * s, -1,
         ifelse(x >= m + sd_multiplier * s, 1, 0))
}

build_consensus <- function(discrete_mat, consensus_threshold = 0.8) {
  apply(discrete_mat, 1, function(x) {
    x <- x[!is.na(x)]
    tab <- table(x)
    winner <- names(which.max(tab))
    if (max(tab) / length(x) >= consensus_threshold) as.integer(winner) else 0
  })
}

# -----------------------------------------------------------------------------
# Load expression matrix and metadata
# -----------------------------------------------------------------------------
expr <- read.delim(expression_file, check.names = FALSE)
rownames(expr) <- expr[[1]]
expr <- expr[, -1, drop = FALSE]
expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"

meta <- read.csv(metadata_file, check.names = FALSE)
meta$sample_id <- as.character(meta$sample_id)
meta$state <- as.character(meta$state)
rownames(meta) <- meta$sample_id

expr <- expr[, meta$sample_id, drop = FALSE]
meta <- meta[colnames(expr), , drop = FALSE]
meta$state <- factor(meta$state, levels = c(target_label, source_label))

# -----------------------------------------------------------------------------
# (Optional) Keep only mRNA / protein-coding genes
# -----------------------------------------------------------------------------
if (use_mrna_filter) {
  mrna_genes <- as.character(read.csv(mrna_gene_file, check.names = FALSE)$Gene_ID)
  expr <- expr[rownames(expr) %in% mrna_genes, , drop = FALSE]
}

# -----------------------------------------------------------------------------
# (Optional) Convert gene identifiers to GEM-compatible IDs
# -----------------------------------------------------------------------------
if (use_id_conversion) {
  id_map <- read.csv(id_mapping_file, check.names = FALSE)
  id_map$Input_ID <- as.character(id_map$Input_ID)
  id_map$Gene_ID <- as.character(id_map$Gene_ID)

  new_ids <- id_map$Gene_ID[match(rownames(expr), id_map$Input_ID)]
  expr <- expr[!is.na(new_ids), , drop = FALSE]
  rownames(expr) <- new_ids[!is.na(new_ids)]
}

# -----------------------------------------------------------------------------
# (Optional) Aggregate duplicate gene IDs after conversion
# -----------------------------------------------------------------------------
if (use_id_conversion) {
  expr <- rowsum(expr, group = rownames(expr))
}

# -----------------------------------------------------------------------------
# Define source and target states
# -----------------------------------------------------------------------------
source_samples <- rownames(meta)[meta$state == source_label]
target_samples <- rownames(meta)[meta$state == target_label]

# -----------------------------------------------------------------------------
# Differential expression analysis: source vs target
# -----------------------------------------------------------------------------
count_mat <- round(expr)
storage.mode(count_mat) <- "integer"

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta,
  design = ~ state
)

dds <- DESeq(dds)
res <- results(dds, contrast = c("state", source_label, target_label))

deg_all <- as.data.frame(res) %>%
  rownames_to_column("Gene_ID") %>%
  rename(logFC = log2FoldChange, pval = pvalue) %>%
  arrange(pval)

write.csv(
  deg_all,
  file.path(out_dir, "DEG_source_vs_target_all_genes.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Normalize expression data for target reference construction
# -----------------------------------------------------------------------------
vsd <- vst(dds, blind = FALSE)
norm_expr <- assay(vsd)

# -----------------------------------------------------------------------------
# Filter DEG table and normalized expression matrix to the GEM gene universe
# -----------------------------------------------------------------------------
gem_genes <- as.character(read.csv(gem_gene_file, check.names = FALSE)$Gene_ID)

deg_gem_all <- deg_all %>%
  filter(Gene_ID %in% gem_genes)

deg_gem_pval <- deg_gem_all %>%
  filter(!is.na(pval), pval <= pval_cutoff)

norm_expr_gem <- norm_expr[rownames(norm_expr) %in% gem_genes, , drop = FALSE]

# -----------------------------------------------------------------------------
# Discretize target-state expression
# -----------------------------------------------------------------------------
target_norm <- norm_expr_gem[, target_samples, drop = FALSE]

target_discrete <- apply(
  target_norm,
  2,
  discretize_sample,
  sd_multiplier = sd_multiplier
)

target_discrete <- as.data.frame(target_discrete)
rownames(target_discrete) <- rownames(target_norm)

# -----------------------------------------------------------------------------
# Build target-state consensus reference
# -----------------------------------------------------------------------------
target_consensus <- build_consensus(target_discrete, consensus_threshold)

target_reference <- data.frame(
  Gene_ID = rownames(target_discrete),
  Expression = target_consensus,
  row.names = NULL
)

write.table(
  target_reference,
  file = file.path(out_dir, "target_state_consensus_reference.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Export summary table
# -----------------------------------------------------------------------------
summary_table <- data.frame(
  Item = c(
    "Source samples",
    "Target samples",
    "Genes used in DESeq2",
    "Genes in GEM universe",
    "GEM genes with DE result",
    "GEM genes with p-value-filtered DE result",
    "Consensus low genes",
    "Consensus moderate genes",
    "Consensus high genes"
  ),
  Value = c(
    length(source_samples),
    length(target_samples),
    nrow(expr),
    length(gem_genes),
    nrow(deg_gem_all),
    nrow(deg_gem_pval),
    sum(target_reference$Expression == -1),
    sum(target_reference$Expression == 0),
    sum(target_reference$Expression == 1)
  )
)

write.csv(summary_table, file.path(out_dir, "R_preparation_summary.csv"), row.names = FALSE)
print(summary_table)
