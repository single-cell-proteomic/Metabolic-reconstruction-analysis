# rMTA result aggregation, filtering, enrichment analysis, and plotting
# Run from the repository root:
# Rscript src/analyze_rMTA_results.R

library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(openxlsx)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(ggplot2)
library(grid)

# If this script is run while the working directory is src/, move to repo root.
if (basename(getwd()) == "src" && dir.exists("../results")) {
  setwd("..")
}

# -----------------------------
# User-defined parameters
# -----------------------------
comparison_name <- "source_vs_target"

rmta_dir <- "results/MATLAB_outputs"
analysis_dir <- "results/analysis_results"
plot_dir <- "plots"

rmta_sheet <- "alpha = 0.66"
rmta_file_pattern <- paste0("^", comparison_name, "_rep\\d+_padj005_FC0_CHRR_rMTA\\.xlsx$")

rank_variation_threshold <- 10
min_valid_experiments <- 8

require_positive_rTS <- TRUE
invalid_rTS_value <- 65535

partitions <- c(1, 0.75, 0.50, 0.25)

gene_symbol_column <- "gene"

org_db <- org.Mm.eg.db
org_keytype <- "SYMBOL"
kegg_organism <- "mmu"

top_n_terms <- 50

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helper functions
# -----------------------------
parse_rts <- function(x) {
  readr::parse_number(
    as.character(x),
    locale = readr::locale(decimal_mark = ".")
  )
}

convert_gene_ratio <- function(x) {
  if (is.numeric(x)) {
    return(x)
  }

  vapply(
    x,
    function(z) {
      if (is.na(z)) {
        return(NA_real_)
      }

      sp <- strsplit(as.character(z), "/", fixed = TRUE)[[1]]

      if (length(sp) != 2) {
        return(NA_real_)
      }

      as.numeric(sp[1]) / as.numeric(sp[2])
    },
    numeric(1)
  )
}

normalize_partition <- function(x) {
  as.character(as.numeric(as.character(x)))
}

get_partition_n_labels <- function(ora_genes, partitions) {
  partition_levels <- normalize_partition(partitions)

  partition_n <- ora_genes %>%
    mutate(
      partition = normalize_partition(partition),
      GeneSymbol = as.character(GeneSymbol)
    ) %>%
    filter(
      !is.na(partition),
      !is.na(GeneSymbol),
      GeneSymbol != ""
    ) %>%
    group_by(partition) %>%
    summarise(
      n_genes = n_distinct(GeneSymbol),
      .groups = "drop"
    )

  label_map <- setNames(
    paste0(partition_n$partition, "\n(n = ", partition_n$n_genes, ")"),
    partition_n$partition
  )

  partition_labels <- setNames(partition_levels, partition_levels)

  for (p in partition_levels) {
    if (p %in% names(label_map)) {
      partition_labels[p] <- label_map[p]
    } else {
      partition_labels[p] <- paste0(p, "\n(n = 0)")
    }
  }

  partition_labels
}

read_and_tidy_enrichment <- function(df, top_n = 50) {
  if (nrow(df) == 0) {
    return(df)
  }

  df$GeneRatio <- convert_gene_ratio(df$GeneRatio)

  df <- df %>%
    mutate(
      Cluster = normalize_partition(Cluster),
      Description = as.character(Description)
    ) %>%
    filter(
      is.finite(GeneRatio),
      is.finite(p.adjust),
      !is.na(Cluster),
      !is.na(Description)
    )

  if (nrow(df) == 0) {
    return(df)
  }

  df$Description <- ifelse(
    nchar(df$Description) > 60,
    paste0(substr(df$Description, 1, 60), "..."),
    df$Description
  )

  df <- df %>% arrange(p.adjust)

  keep_terms <- unique(df$Description)[
    1:min(top_n, length(unique(df$Description)))
  ]

  df %>% filter(Description %in% keep_terms)
}

get_term_order <- function(df) {
  df %>%
    group_by(Description) %>%
    summarise(
      best_padj = min(p.adjust, na.rm = TRUE),
      max_gene_ratio = max(GeneRatio, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(best_padj, desc(max_gene_ratio)) %>%
    pull(Description)
}

get_horizontal_dims <- function(df) {
  n_terms <- length(unique(df$Description))
  list(
    width = max(11, min(36, 5 + 0.70 * n_terms)),
    height = 6.5
  )
}

get_vertical_dims <- function(df) {
  n_terms <- length(unique(df$Description))
  list(
    width = 12,
    height = max(7, min(34, 3.5 + 0.45 * n_terms))
  )
}

plot_horizontal <- function(df, plot_title, partition_labels, partitions) {
  term_order <- get_term_order(df)

  df_plot <- df %>%
    mutate(
      x_term = factor(Description, levels = term_order),
      Cluster = factor(Cluster, levels = normalize_partition(partitions))
    )

  ggplot(
    df_plot,
    aes(
      x = x_term,
      y = Cluster,
      size = GeneRatio,
      color = p.adjust
    )
  ) +
    geom_point(alpha = 0.95) +
    scale_size(range = c(3, 10)) +
    scale_color_gradient(low = "red", high = "blue", trans = "reverse") +
    scale_y_discrete(labels = partition_labels) +
    labs(
      title = plot_title,
      x = "",
      y = "Partition",
      size = "GeneRatio",
      color = "p.adjust"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0),
      axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 12),
      legend.position = "top",
      legend.direction = "horizontal",
      plot.margin = margin(t = 8, r = 14, b = 8, l = 8)
    ) +
    guides(
      color = guide_colorbar(
        title.position = "top",
        barwidth = grid::unit(5, "cm"),
        barheight = grid::unit(0.5, "cm")
      ),
      size = guide_legend(
        title.position = "top",
        nrow = 1
      )
    )
}

plot_vertical <- function(df, plot_title, partition_labels, partitions) {
  term_order <- get_term_order(df)

  df_plot <- df %>%
    mutate(
      Cluster = factor(Cluster, levels = normalize_partition(partitions)),
      y_term = factor(Description, levels = rev(term_order))
    )

  ggplot(
    df_plot,
    aes(
      x = Cluster,
      y = y_term,
      size = GeneRatio,
      color = p.adjust
    )
  ) +
    geom_point(alpha = 0.95) +
    scale_size(range = c(3, 10)) +
    scale_color_gradient(low = "red", high = "blue", trans = "reverse") +
    scale_x_discrete(labels = partition_labels) +
    labs(
      title = plot_title,
      x = "Partition",
      y = "",
      size = "GeneRatio",
      color = "p.adjust"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.position = "top",
      legend.direction = "horizontal",
      plot.margin = margin(t = 8, r = 18, b = 8, l = 8)
    ) +
    guides(
      color = guide_colorbar(
        title.position = "top",
        barwidth = grid::unit(5, "cm"),
        barheight = grid::unit(0.5, "cm")
      ),
      size = guide_legend(
        title.position = "top",
        nrow = 1
      )
    )
}

save_enrichment_plots <- function(enrichment_table, analysis_name, ora_genes) {
  plot_table <- read_and_tidy_enrichment(
    enrichment_table,
    top_n = top_n_terms
  )

  if (nrow(plot_table) == 0) {
    message("No ", analysis_name, " enrichment rows to plot.")
    return(invisible(NULL))
  }

  partition_labels <- get_partition_n_labels(ora_genes, partitions)

  plot_title <- paste0(
    comparison_name,
    " ",
    analysis_name,
    " enrichment | var=",
    rank_variation_threshold,
    ", minValid=",
    min_valid_experiments
  )

  p_h <- plot_horizontal(
    plot_table,
    plot_title,
    partition_labels,
    partitions
  )

  p_v <- plot_vertical(
    plot_table,
    plot_title,
    partition_labels,
    partitions
  )

  dims_h <- get_horizontal_dims(plot_table)
  dims_v <- get_vertical_dims(plot_table)

  out_h <- file.path(
    plot_dir,
    paste0(
      comparison_name,
      "_",
      analysis_name,
      "_dotplot_horizontal_var",
      rank_variation_threshold,
      "_minValid",
      min_valid_experiments,
      ".png"
    )
  )

  out_v <- file.path(
    plot_dir,
    paste0(
      comparison_name,
      "_",
      analysis_name,
      "_dotplot_vertical_var",
      rank_variation_threshold,
      "_minValid",
      min_valid_experiments,
      ".png"
    )
  )

  ggsave(
    filename = out_h,
    plot = p_h,
    width = dims_h$width,
    height = dims_h$height,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  ggsave(
    filename = out_v,
    plot = p_v,
    width = dims_v$width,
    height = dims_v$height,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  message("Saved plot: ", out_h)
  message("Saved plot: ", out_v)

  invisible(NULL)
}

# -----------------------------
# 1. Merge rMTA replicate results
# -----------------------------
rmta_files <- list.files(
  rmta_dir,
  pattern = rmta_file_pattern,
  full.names = TRUE
)

rmta_files <- sort(rmta_files)

if (length(rmta_files) == 0) {
  stop("No rMTA result files found in: ", rmta_dir)
}

message("rMTA files found: ", length(rmta_files))

all_runs <- list()

for (i in seq_along(rmta_files)) {
  file_i <- rmta_files[i]
  run_id <- sprintf("%02d", i)

  message("Reading rMTA file: ", file_i)

  rmta_res <- readxl::read_excel(
    file_i,
    sheet = rmta_sheet
  )

  required_cols <- c("gene", "rTS")
  missing_cols <- setdiff(required_cols, colnames(rmta_res))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in file: ",
      file_i,
      "\nMissing: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  rmta_res <- rmta_res %>%
    mutate(
      gene = trimws(as.character(gene)),
      rTS = parse_rts(rTS)
    ) %>%
    filter(
      !is.na(gene),
      gene != "",
      !is.na(rTS)
    )

  if (gene_symbol_column %in% colnames(rmta_res)) {
    rmta_res$GeneSymbol <- trimws(as.character(rmta_res[[gene_symbol_column]]))
  } else {
    rmta_res$GeneSymbol <- rmta_res$gene
  }

  valid_rows <- rmta_res %>%
    filter(rTS != invalid_rTS_value) %>%
    arrange(desc(rTS))

  invalid_rows <- rmta_res %>%
    filter(rTS == invalid_rTS_value)

  rmta_res <- bind_rows(valid_rows, invalid_rows) %>%
    mutate(
      run_id = run_id,
      rank = row_number()
    )

  all_runs[[paste0("run", run_id)]] <- rmta_res
}

all_genes <- Reduce(
  union,
  lapply(all_runs, function(x) x$gene)
)

merged_ranks <- data.frame(
  gene = all_genes,
  stringsAsFactors = FALSE
)

first_symbol_table <- bind_rows(all_runs) %>%
  select(gene, GeneSymbol) %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
  distinct(gene, .keep_all = TRUE)

merged_ranks <- merged_ranks %>%
  left_join(first_symbol_table, by = "gene")

for (i in seq_along(all_runs)) {
  run_df <- all_runs[[i]]

  merged_ranks[[paste0("rTS", i)]] <- run_df$rTS[
    match(merged_ranks$gene, run_df$gene)
  ]

  merged_ranks[[paste0("rank", i)]] <- run_df$rank[
    match(merged_ranks$gene, run_df$gene)
  ]
}

rts_cols <- grep("^rTS", colnames(merged_ranks), value = TRUE)
rank_cols <- grep("^rank", colnames(merged_ranks), value = TRUE)

rts_matrix <- as.matrix(merged_ranks[, rts_cols, drop = FALSE])
rank_matrix <- as.matrix(merged_ranks[, rank_cols, drop = FALSE])

valid_matrix <- !is.na(rts_matrix) & rts_matrix != invalid_rTS_value

merged_ranks$valid_experiment_count <- rowSums(valid_matrix, na.rm = TRUE)

rank_matrix_valid <- rank_matrix
rank_matrix_valid[!valid_matrix] <- NA

merged_ranks$valid_min <- apply(
  rank_matrix_valid,
  1,
  function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
)

merged_ranks$valid_max <- apply(
  rank_matrix_valid,
  1,
  function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
)

merged_ranks$valid_median <- apply(
  rank_matrix_valid,
  1,
  function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
)

merged_ranks$valid_mean <- apply(
  rank_matrix_valid,
  1,
  function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
)

merged_ranks$valid_std <- apply(
  rank_matrix_valid,
  1,
  function(x) {
    if (sum(!is.na(x)) <= 1) NA_real_ else sd(x, na.rm = TRUE)
  }
)

merged_ranks$valid_rank_range <- merged_ranks$valid_max - merged_ranks$valid_min

merged_file <- file.path(
  analysis_dir,
  paste0(comparison_name, "_rMTA_merged_rank_statistics.csv")
)

write.csv(
  merged_ranks,
  merged_file,
  row.names = FALSE
)

message("Saved merged rMTA statistics: ", merged_file)

# -----------------------------
# 2. Filter retained rMTA genes
# -----------------------------
positive_valid_rTS <- apply(
  rts_matrix,
  1,
  function(x) {
    x_valid <- x[!is.na(x) & x != invalid_rTS_value]

    if (length(x_valid) == 0) {
      return(FALSE)
    }

    all(x_valid > 0)
  }
)

keep_rows <- merged_ranks$valid_experiment_count >= min_valid_experiments &
  !is.na(merged_ranks$valid_rank_range) &
  merged_ranks$valid_rank_range <= rank_variation_threshold

if (require_positive_rTS) {
  keep_rows <- keep_rows & positive_valid_rTS
}

retained_genes <- merged_ranks %>%
  filter(keep_rows) %>%
  arrange(valid_median)

if (nrow(retained_genes) == 0) {
  stop("No genes retained after rMTA filtering. Relax filtering parameters.")
}

retained_file <- file.path(
  analysis_dir,
  paste0(
    comparison_name,
    "_rMTA_retained_genes_var",
    rank_variation_threshold,
    "_minValid",
    min_valid_experiments,
    ".csv"
  )
)

write.csv(
  retained_genes,
  retained_file,
  row.names = FALSE
)

message("Saved retained rMTA genes: ", retained_file)

# -----------------------------
# 3. Build retained-gene partitions
# -----------------------------
partition_tables <- list()

for (partition in partitions) {
  n_retained <- nrow(retained_genes)
  cutoff <- floor(partition * n_retained)

  if (cutoff <= 0) {
    next
  }

  tmp_partition <- retained_genes %>%
    arrange(valid_median) %>%
    slice(1:cutoff) %>%
    mutate(
      partition = normalize_partition(partition)
    )

  partition_tables[[normalize_partition(partition)]] <- tmp_partition
}

retained_partitions <- bind_rows(partition_tables)

partition_file <- file.path(
  analysis_dir,
  paste0(
    comparison_name,
    "_rMTA_retained_gene_partitions_var",
    rank_variation_threshold,
    "_minValid",
    min_valid_experiments,
    ".csv"
  )
)

write.csv(
  retained_partitions,
  partition_file,
  row.names = FALSE
)

message("Saved retained gene partitions: ", partition_file)

# -----------------------------
# 4. GO-BP and KEGG enrichment
# -----------------------------
ora_genes <- retained_partitions %>%
  select(GeneSymbol, partition, valid_experiment_count) %>%
  mutate(
    GeneSymbol = as.character(GeneSymbol),
    partition = as.character(partition)
  ) %>%
  filter(
    !is.na(GeneSymbol),
    GeneSymbol != "",
    !is.na(partition)
  ) %>%
  distinct()

symbol_keys <- unique(ora_genes$GeneSymbol)

symbol_to_entrez <- AnnotationDbi::select(
  org_db,
  keys = symbol_keys,
  keytype = org_keytype,
  columns = c(org_keytype, "ENTREZID")
)

colnames(symbol_to_entrez)[colnames(symbol_to_entrez) == org_keytype] <- "GeneSymbol"

symbol_to_entrez <- symbol_to_entrez %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(GeneSymbol, ENTREZID)

ora_genes_mapped <- ora_genes %>%
  left_join(symbol_to_entrez, by = "GeneSymbol") %>%
  filter(!is.na(ENTREZID)) %>%
  distinct()

if (nrow(ora_genes_mapped) == 0) {
  stop("No retained genes could be mapped to ENTREZID.")
}

partition_gene_lists <- ora_genes_mapped %>%
  split(.$partition) %>%
  lapply(function(x) unique(as.character(x$ENTREZID)))

partition_gene_lists <- partition_gene_lists[vapply(partition_gene_lists, length, numeric(1)) > 0]

if (length(partition_gene_lists) == 0) {
  stop("No non-empty partition gene lists available for enrichment analysis.")
}

message("Running GO-BP enrichment...")

go_results <- compareCluster(
  geneCluster = partition_gene_lists,
  fun = "enrichGO",
  OrgDb = org_db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

message("Running KEGG enrichment...")

kegg_results <- compareCluster(
  geneCluster = partition_gene_lists,
  fun = "enrichKEGG",
  organism = kegg_organism,
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

go_table <- as.data.frame(go_results)
kegg_table <- as.data.frame(kegg_results)

go_file <- file.path(
  analysis_dir,
  paste0(
    comparison_name,
    "_GO_BP_var",
    rank_variation_threshold,
    "_minValid",
    min_valid_experiments,
    ".csv"
  )
)

kegg_file <- file.path(
  analysis_dir,
  paste0(
    comparison_name,
    "_KEGG_var",
    rank_variation_threshold,
    "_minValid",
    min_valid_experiments,
    ".csv"
  )
)

ora_input_file <- file.path(
  analysis_dir,
  paste0(
    comparison_name,
    "_ORA_input_genes_var",
    rank_variation_threshold,
    "_minValid",
    min_valid_experiments,
    ".csv"
  )
)

ora_mapped_file <- file.path(
  analysis_dir,
  paste0(
    comparison_name,
    "_ORA_mapped_entrez_var",
    rank_variation_threshold,
    "_minValid",
    min_valid_experiments,
    ".csv"
  )
)

write.csv(go_table, go_file, row.names = FALSE)
write.csv(kegg_table, kegg_file, row.names = FALSE)
write.csv(ora_genes, ora_input_file, row.names = FALSE)
write.csv(ora_genes_mapped, ora_mapped_file, row.names = FALSE)

message("Saved GO-BP enrichment results: ", go_file)
message("Saved KEGG enrichment results: ", kegg_file)
message("Saved ORA input genes: ", ora_input_file)
message("Saved mapped ORA genes: ", ora_mapped_file)

# -----------------------------
# 5. Plot enrichment results
# -----------------------------
save_enrichment_plots(
  enrichment_table = go_table,
  analysis_name = "GO_BP",
  ora_genes = ora_genes
)

save_enrichment_plots(
  enrichment_table = kegg_table,
  analysis_name = "KEGG",
  ora_genes = ora_genes
)

message("rMTA result analysis completed.")
