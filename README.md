# Pairwise Metabolic Transformation Analysis Pipeline

This repository contains a two-step pipeline for pairwise metabolic transformation analysis between a **source** state and a **target** state.

The pipeline is intentionally simple and handover-friendly:

1. The **R section** prepares differential expression and target-state expression-reference files.
2. The **MATLAB section** uses those files to build a target-state tissue-specific model and run rMTA.

The default workflow assumes raw RNA-seq gene counts, but the data-format requirements and points that need adaptation for other expression platforms are described below.

## Repository layout

```text
.
├── README.md
├── data/
│   ├── expression_counts.tsv
│   ├── metadata.csv
│   ├── gem_genes.csv
│   ├── mrna_genes.csv              # (Optional)
│   ├── id_mapping.csv              # (Optional)
│   └── Recon3D.mat                 # or another GEM model
├── results/
│   ├── R_outputs/
│   └── MATLAB_outputs/
└── src/
    ├── pairwise_prepare_inputs.R
    └── pairwise_rMTA_general.m
```

The `README.md` explains the workflow step by step. The `src/` files contain the same logic as directly executable scripts.

---

# Part I — R-side preparation

## 1. Purpose

The R side prepares MATLAB/rMTA-ready input files from an expression matrix and metadata table.

It performs the following tasks:

1. loads the expression matrix and metadata,
2. (Optional) removes non-mRNA / non-protein-coding genes,
3. (Optional) converts gene identifiers to the GEM-compatible ID type,
4. computes differential expression as **source versus target**,
5. normalizes raw RNA-seq counts using VST for target-reference construction,
6. filters exported files to the GEM gene universe,
7. discretizes target-state expression,
8. builds a target-state consensus expression reference,
9. exports MATLAB/rMTA-ready files.

The R code is written for raw RNA-seq counts. For microarray or already normalized data, see [(Optional) Notes for non-RNA-seq datasets](#optional-notes-for-non-rna-seq-datasets).

## 2. Required input files

### 2.1 Expression matrix

The expression matrix must contain genes in rows and samples in columns.

| Gene_ID | Sample_01 | Sample_02 | Sample_03 |
|---|---:|---:|---:|
| GeneA | 120 | 95 | 210 |
| GeneB | 0 | 4 | 1 |
| GeneC | 51 | 43 | 38 |

For the default RNA-seq workflow, values must be raw, non-normalized gene counts.

**Note on example data**

The `data/expression_counts.tsv` file included in this repository is a synthetic RNA-seq count matrix generated only for testing and demonstration purposes. It does not represent real biological samples. Users should replace this file with their own raw RNA-seq count matrix before running a real analysis.

### 2.2 Metadata table

The metadata table must contain at least two columns.

| sample_id | state |
|---|---|
| Sample_01 | target |
| Sample_02 | source |
| Sample_03 | source |

Column meanings:

| Column | Meaning |
|---|---|
| `sample_id` | Sample identifier matching the expression matrix column names. |
| `state` | Biological state of each sample: `source` or `target`. |
| `source` | Initial / old / patient / disease state. |
| `target` | Desired / young / healthy state. |

Differential expression is always computed as **source state versus target state**.

> **Direction rule**
>
> Positive `logFC` = higher expression in the **source** state.  
> Negative `logFC` = higher expression in the **target** state.

For example, if `source = old / patient / disease` and `target = young / healthy / desired state`, positive `logFC` values represent genes increased in the source state relative to the target state.

### 2.3 (Optional) mRNA / protein-coding gene list

Use this file only if the expression matrix contains ncRNAs, pseudogenes, antisense transcripts, or other features that should not be used in the GEM/rMTA analysis.

| Gene_ID |
|---|
| GeneA |
| GeneB |
| GeneC |

The IDs in this file must match the current row identifiers of the expression matrix at the moment this filter is applied.

### 2.4 (Optional) Gene ID mapping table

Use this file if the expression matrix is indexed by a different gene identifier type than the GEM model.

| Input_ID | Gene_ID |
|---|---|
| ENSG00000141510 | 7157 |
| ENSG00000171862 | 2064 |

Column meanings:

| Column | Meaning |
|---|---|
| `Input_ID` | Gene ID currently used in the expression matrix. |
| `Gene_ID` | Gene ID required by the GEM model. |

The final `Gene_ID` type must match the gene IDs used by the GEM model. For example, if Recon3D is used with Entrez IDs, Ensembl-indexed RNA-seq data should be converted to Entrez IDs before exporting MATLAB inputs.

### 2.5 GEM gene universe

The GEM gene universe is the list of genes present in the metabolic model.

| Gene_ID |
|---|
| 7157 |
| 2064 |
| ... |

DE analysis is **not** performed only on GEM genes. DE analysis is performed on all usable expression genes, and GEM filtering is applied later to the exported DEG table and target reference files.

## 3. R packages

```r
library(DESeq2)
library(dplyr)
library(tibble)
```

## 4. User-defined parameters

Edit this block before running the pipeline on a new dataset.

```r
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
```

## 5. Load expression matrix and metadata

The expression matrix and metadata are aligned by sample ID. The metadata `state` column is encoded so that the DESeq2 contrast can be explicitly defined as **source versus target**.

```r
expr <- read.delim(expression_file, check.names = FALSE)
rownames(expr) <- expr[[1]]
expr <- expr[, -1]
expr <- as.matrix(expr)

meta <- read.csv(metadata_file, check.names = FALSE)
meta$sample_id <- as.character(meta$sample_id)
meta$state <- as.character(meta$state)
rownames(meta) <- meta$sample_id

expr <- expr[, meta$sample_id]
meta <- meta[colnames(expr), ]
meta$state <- factor(meta$state, levels = c(target_label, source_label))
```

## 6. (Optional) Keep only mRNA / protein-coding genes

Run this step if the expression matrix contains ncRNAs or other non-mRNA features that should be excluded before DE analysis and GEM mapping.

```r
if (use_mrna_filter) {
  mrna_genes <- as.character(read.csv(mrna_gene_file, check.names = FALSE)$Gene_ID)
  expr <- expr[rownames(expr) %in% mrna_genes, ]
}
```

## 7. (Optional) Convert gene identifiers to GEM-compatible IDs

Run this step if the expression matrix uses a different gene ID type than the GEM model.

After this step, row names of the expression matrix become model-compatible `Gene_ID` values.

```r
if (use_id_conversion) {
  id_map <- read.csv(id_mapping_file, check.names = FALSE)
  id_map$Input_ID <- as.character(id_map$Input_ID)
  id_map$Gene_ID <- as.character(id_map$Gene_ID)

  new_ids <- id_map$Gene_ID[match(rownames(expr), id_map$Input_ID)]

  expr <- expr[!is.na(new_ids), ]
  rownames(expr) <- new_ids[!is.na(new_ids)]
}
```

## 8. (Optional) Aggregate duplicate gene IDs after conversion

If multiple input genes map to the same GEM-compatible `Gene_ID`, duplicate rows must be aggregated before DESeq2.

For raw RNA-seq counts, duplicates are summed. For microarray or already normalized expression values, mean or median aggregation may be more appropriate.

```r
if (use_id_conversion) {
  expr <- rowsum(expr, group = rownames(expr))
}
```

## 9. Define source and target states

This pipeline expects exactly two biological states: one source state and one target state.

```r
source_samples <- rownames(meta)[meta$state == source_label]
target_samples <- rownames(meta)[meta$state == target_label]
```

## 10. Differential expression analysis: source vs target

This step is written for raw RNA-seq counts.

DE analysis is performed before filtering to the GEM gene universe.

```r
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
```

## 11. Normalize expression data for target reference construction

For raw RNA-seq counts, VST-normalized expression values are used for target-state discretization and consensus reference construction.

```r
vsd <- vst(dds, blind = FALSE)
norm_expr <- assay(vsd)
```

## 12. Filter DEG table and normalized expression matrix to the GEM gene universe

GEM filtering is applied after DE analysis and normalization.

The R-side p-value filter is a practical pre-filter. The MATLAB/rMTA side still uses `padj` through the `padj_threshold` parameter.

```r
gem_genes <- as.character(read.csv(gem_gene_file, check.names = FALSE)$Gene_ID)

deg_gem_all <- deg_all %>%
  filter(Gene_ID %in% gem_genes)

deg_gem_pval <- deg_gem_all %>%
  filter(!is.na(pval), pval <= pval_cutoff)

norm_expr_gem <- norm_expr[rownames(norm_expr) %in% gem_genes, ]

write.csv(
  deg_gem_all,
  file.path(out_dir, "DEG_source_vs_target_GEM_genes_all.csv"),
  row.names = FALSE
)

write.csv(
  deg_gem_pval,
  file.path(out_dir, paste0("DEG_source_vs_target_GEM_genes_pval_", pval_cutoff, ".csv")),
  row.names = FALSE
)

write.table(
  norm_expr_gem,
  file = file.path(out_dir, "normalized_expression_GEM_genes_VST.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
```

### DEG filtering note for rMTA

The exported p-value-filtered DEG table can be further filtered inside the MATLAB/rMTA workflow by adjusted p-value, fold change, or top-ranked DEG count.

In the original rMTA-style workflow, the selected DEG set is often adjusted so that the number of active reactions remains approximately between 100 and 200. However, this criterion is not strongly justified as a universal rule. It should be optimized according to the dataset, biological contrast, and desired transformation behavior.

If too many DEGs are passed to rMTA, too many reactions may become active. In that case, knock-out simulations may fail to produce negative `wTS` values, which can prevent `wTS` and `bTS` from contributing meaningfully to the final `rTS` score. Therefore, if `wTS`/`bTS`-based scoring is desired, it can be useful to keep the DEG set relatively compact by using stricter adjusted p-value and/or fold-change thresholds.

## 13. Discretize target-state expression

Only target-state samples are used to generate the desired expression reference.

Each target sample is discretized gene-wise using the expression distribution of that sample:

- low expression = `-1` if expression `<= mean - 0.6 × SD`
- moderate expression = `0` if expression is between thresholds
- high expression = `1` if expression `>= mean + 0.6 × SD`

The default `0.6 × SD` threshold follows commonly used iMAT-style expression discretization settings. This value can be optimized depending on the expression distribution and expected biological structure of the target state.

```r
discretize_sample <- function(x, sd_multiplier = 0.6) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  ifelse(x <= m - sd_multiplier * s, -1,
         ifelse(x >= m + sd_multiplier * s, 1, 0))
}

target_norm <- norm_expr_gem[, target_samples, drop = FALSE]

target_discrete <- apply(
  target_norm,
  2,
  discretize_sample,
  sd_multiplier = sd_multiplier
)

target_discrete <- as.data.frame(target_discrete)
rownames(target_discrete) <- rownames(target_norm)
```

## 14. Build target-state consensus reference

The consensus reference summarizes the desired target state across target samples.

For each gene:

- if at least `consensus_threshold` fraction of target samples agree on the same discrete expression state, that value is assigned;
- otherwise, the consensus value is set to `0`.

A higher threshold gives a stricter and more homogeneous target reference. A lower threshold may be useful for heterogeneous cohorts or tissues.

```r
build_consensus <- function(discrete_mat, consensus_threshold = 0.8) {
  apply(discrete_mat, 1, function(x) {
    x <- x[!is.na(x)]
    tab <- table(x)
    winner <- names(which.max(tab))
    if (max(tab) / length(x) >= consensus_threshold) as.integer(winner) else 0
  })
}

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
```

## 15. Export R-side summary table

```r
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
summary_table
```

## 16. R-side output files

The most important MATLAB/rMTA input files are:

| File | Purpose |
|---|---|
| `results/R_outputs/DEG_source_vs_target_GEM_genes_pval_0.05.csv` | GEM-filtered source-vs-target DEG table. |
| `results/R_outputs/target_state_consensus_reference.tsv` | Target-state consensus expression reference. |

The DEG file uses this format:

| Column | Meaning |
|---|---|
| `Gene_ID` | GEM-compatible gene identifier. |
| `baseMean` | Mean normalized count from DESeq2. |
| `logFC` | Source-vs-target log2 fold change. |
| `lfcSE` | Standard error of log fold change. |
| `stat` | DESeq2 test statistic. |
| `pval` | Raw p-value. |
| `padj` | Adjusted p-value. |

The target reference file uses this format:

| Column | Meaning |
|---|---|
| `Gene_ID` | GEM-compatible gene identifier. |
| `Expression` | Discrete target-state consensus expression value. |

`Expression` values:

| Value | Meaning |
|---:|---|
| `-1` | Low expression in target state. |
| `0` | Moderate / uncertain / no consensus. |
| `1` | High expression in target state. |

## (Optional) Notes for non-RNA-seq datasets

This R workflow is written for raw RNA-seq gene count data.

For microarray data or already normalized expression matrices:

- Do not use DESeq2.
- Do not use VST normalization.
- Use an appropriate differential expression method, such as limma.
- If the expression values are already normalized and approximately continuous, pass the normalized matrix directly to the discretization and consensus steps.
- For microarray-like expression values that are already normalized and approximately normally distributed, additional VST-like transformation is usually unnecessary.

The key requirements remain the same regardless of data type:

1. DEG direction must be source versus target.
2. Positive `logFC` must mean higher expression in the source state.
3. The final `Gene_ID` type must match the GEM model.
4. The target consensus reference must be generated only from target-state samples.


---

# ======================================================================
# Part II — MATLAB/rMTA transformation analysis
# ======================================================================

## 1. Purpose

The MATLAB side continues from the R-side outputs. It expects:

1. one target-state consensus expression file,
2. one source-vs-target DEG file,
3. one GEM model compatible with the `Gene_ID` values.

It then:

1. loads the GEM model,
2. aligns target expression with model genes,
3. maps target gene expression to reaction expression,
4. builds a target-state tissue-specific model using iMAT,
5. samples the target-state model to estimate `Vref`,
6. maps source-vs-target DEG information to reaction-level `rxnFBS`,
7. runs rMTA for multiple sampling replicates,
8. exports replicate-level Excel and MAT outputs.

The target-state model represents the desired / healthy / young metabolic reference. rMTA evaluates source-state perturbation signals relative to this target reference state.

## 2. Required MATLAB / COBRA dependencies

Required software and functions:

- MATLAB
- COBRA Toolbox
- Gurobi or another compatible MILP solver
- rMTA helper functions:
  - `diffexprs2rxnFBS`
  - `calculateEPSILON`
  - `rMTA`
  - `rMTAsaveInExcel`

## 3. Required MATLAB input files

The default paths below match the R-side outputs in this README.

| File | Purpose |
|---|---|
| `results/R_outputs/target_state_consensus_reference.tsv` | Target-state consensus expression reference. |
| `results/R_outputs/DEG_source_vs_target_GEM_genes_pval_0.05.csv` | Source-vs-target DEG table. |
| `data/Recon3D.mat` | GEM model used by COBRA/rMTA. |

The target expression file must contain:

| Column | Meaning |
|---|---|
| `Gene_ID` | Gene identifier compatible with the cleaned model gene IDs. |
| `Expression` | Discrete target-state value: `-1`, `0`, or `1`. |

The DEG file must contain at least:

| Column | Meaning |
|---|---|
| `Gene_ID` | Gene identifier compatible with the cleaned model gene IDs. |
| `logFC` | Source-vs-target log2 fold change. |
| `padj` | Adjusted p-value used by the MATLAB/rMTA filter. |

The rMTA helper function expects a column named `pval`. In this pipeline, adjusted p-values are used conceptually and computationally. Therefore, `padj` is copied into a temporary `pval` column only for compatibility with the existing rMTA helper code.

> **Direction rule**
>
> The DEG file must represent **source versus target**.  
> Positive `logFC` = higher expression in source.  
> Negative `logFC` = higher expression in target.

## 4. User-defined MATLAB parameters

```matlab
clear; clc;

model_file = 'data/Recon3D.mat';
target_expression_file = 'results/R_outputs/target_state_consensus_reference.tsv';
deg_file = 'results/R_outputs/DEG_source_vs_target_GEM_genes_pval_0.05.csv';

comparison_name = 'source_vs_target';
out_dir = 'results/MATLAB_outputs';
log_dir = fullfile(out_dir, 'Logs');
mkdir(out_dir);
mkdir(log_dir);

target_gene_id_column = 'Gene_ID';
target_expression_column = 'Expression';
deg_gene_id_column = 'Gene_ID';
deg_logfc_column = 'logFC';
deg_padj_column = 'padj';

fc_threshold = 0;
padj_threshold = 0.05;
n_replicates = 10;

solver_name = 'gurobi';
transcript_separator = '.';
sampling_method = 'CHRR';
alpha_values = 0.66;

imat_time_limit = 60;
rmta_time_limit = 10;
n_workers = 4;
n_threads = 4;
print_level = 1;

sampling_options = struct();
sampling_options.nPointsReturned = 2000;
sampling_options.nStepsPerPoint = 200;
sampling_options.nThreads = n_threads;
```

## 5. Load COBRA model and prepare model gene IDs

Some GEMs store genes with transcript-like suffixes, for example `1234.1`, while the R outputs may contain the parent gene ID, for example `1234`.

This step automatically creates two columns:

| Column | Meaning |
|---|---|
| `Gene_ID` | Cleaned ID used to match R output files. |
| `Model_Gene_ID` | Original ID used internally by COBRA and `model.genes`. |

If the model genes do not contain transcript-like suffixes, both IDs remain identical.

```matlab
changeCobraSolver(solver_name, 'all');
model = readCbModel(model_file);
model_gene_table = makeModelGeneTable(model, transcript_separator);
```

## 6. Load and map target-state consensus expression

The target-state expression file is the desired reference state. Genes not present in the target file are assigned `Expression = 0`, meaning no expression evidence / moderate expression.

The discrete expression values are multiplied by 2 before iMAT reaction mapping so that the `-1 / 0 / 1` input matches the `+/-1` thresholds used in this COBRA/iMAT workflow.

```matlab
[target_data, target_rxn_expression, parsedGPR] = prepareTargetExpression( ...
    model, model_gene_table, target_expression_file, ...
    target_gene_id_column, target_expression_column);

fprintf('\nMapped target expression from gene level to reaction level.\n');
fprintf('Target genes matched to model: %d / %d\n', ...
    sum(target_data.Expression ~= 0), height(target_data));
```

## 7. Build target-state tissue-specific model with iMAT

The iMAT model represents the target/healthy/desired metabolic state.

```matlab
imat_options = struct();
imat_options.solver = 'iMAT';
imat_options.threshold_lb = +1;
imat_options.threshold_ub = -1;
imat_options.expressionRxns = target_rxn_expression;
imat_options.timelimit = imat_time_limit;
imat_options.printLevel = print_level;
imat_options.numWorkers = n_workers;
imat_options.numThreads = n_threads;

tic;
tissueModel_ref = createTissueSpecificModel(model, imat_options, 1);
TIME.iMAT = toc;

save(fullfile(log_dir, [comparison_name '_target_tissueModel_ref.mat']), ...
    'tissueModel_ref', 'imat_options', 'TIME');
```

## 8. Load and prepare source-vs-target DEG table

This pipeline uses adjusted p-values. The rMTA helper function expects a column named `pval`, so `padj` is copied into a temporary `pval` column only for compatibility.

Required internal variables for `diffexprs2rxnFBS`:

| Internal variable | Meaning |
|---|---|
| `gene` | Gene IDs compatible with the cleaned model `Gene_ID`. |
| `logFC` | Source-vs-target log fold-change. |
| `pval` | Adjusted p-value alias copied from `padj`. |

```matlab
differ_genes = prepareDegTable( ...
    deg_file, deg_gene_id_column, deg_logfc_column, deg_padj_column);

fprintf('\nDEG table loaded: %d genes\n', height(differ_genes));
fprintf('Using adjusted p-value threshold: %.4g\n', padj_threshold);
fprintf('Using absolute logFC threshold: %.4g\n', fc_threshold);
```

## 9. Run rMTA replicates

Each replicate repeats flux sampling on the same target tissue-specific model. The sampled mean flux vector `Vref` is used as the target metabolic reference state for rMTA.

Important QC checkpoint:

- After DEG-to-reaction mapping, check `sum(rxnFBS ~= 0)`.
- Very large DEG sets may activate too many reactions.
- If too many reactions are active, knockout simulations may fail to produce negative `wTS` values.
- In that case, `wTS` and `bTS` may not contribute meaningfully to final `rTS`.
- If this happens, stricter `padj` and/or `logFC` thresholds can be used to decrease number of active reactions.

```matlab
for rep = 1:n_replicates
    rep_label = sprintf('rep%02d', rep);
    fprintf('\n=== %s | %s ===\n', comparison_name, rep_label);

    tic;
    [modelSampling, samples_raw] = sampleCbModel( ...
        tissueModel_ref, 'sampleFiles', sampling_method, sampling_options);
    TIME.sampling = toc;

    [samples, sampleStats, Vref, rxnInactive] = expandSamplingToFullModel( ...
        model, tissueModel_ref, samples_raw);

    rxnFBS = diffexprs2rxnFBS(model, differ_genes, Vref, ...
        'SeparateTranscript', transcript_separator, ...
        'logFC', fc_threshold, ...
        'pval', padj_threshold);

    rxnFBS(rxnInactive) = 0;
    active_reaction_count = sum(rxnFBS ~= 0);
    fprintf('Active differential reactions after curation: %d\n', active_reaction_count);

    rxnFBS_table = table(model.rxns(:), rxnFBS(:), ...
        'VariableNames', {'Reaction_ID', 'rxnFBS'});
    save(fullfile(log_dir, [comparison_name '_' rep_label '_rxnFBS.mat']), ...
        'rxnFBS_table', 'active_reaction_count');

    epsilon = calculateEPSILON(samples, rxnFBS);

    changeCobraSolver(solver_name, 'all');
    tic;
    [TSscore, deletedGenes, Vres] = rMTA( ...
        model, rxnFBS, Vref, alpha_values, epsilon, ...
        'timelimit', rmta_time_limit, ...
        'SeparateTranscript', transcript_separator, ...
        'printLevel', print_level, ...
        'numWorkers', n_workers);
    TIME.rMTA = toc;

    gene_info = table(cellstr(string(deletedGenes(:))), 'VariableNames', {'gene'});

    excel_file = fullfile(out_dir, ...
        [comparison_name '_' rep_label '_padj005_FC0_' sampling_method '_rMTA.xlsx']);

    rMTAsaveInExcel(excel_file, TSscore, deletedGenes, alpha_values, ...
        'differ_genes', differ_genes, 'gene_info', gene_info);

    save(fullfile(log_dir, [comparison_name '_' rep_label '_rMTA_workspace.mat']), ...
        'TSscore', 'deletedGenes', 'Vres', 'Vref', 'epsilon', 'rxnFBS', ...
        'active_reaction_count', 'TIME', 'sampleStats', 'modelSampling');
end
```

## 10. Threshold tuning notes

Default thresholds are:

```matlab
fc_threshold = 0;
padj_threshold = 0.05;
```

These are intentionally permissive. rMTA can become less informative if too many DEGs create too many active reaction-level signals. If the number of active differential reactions is too high, consider increasing `fc_threshold` or decreasing `padj_threshold`.

The original rMTA-style workflow often aims to keep the number of active reactions approximately between 100 and 200, but this should not be treated as a universal rule. It should be optimized according to dataset size, biological contrast, and the desired transformation behavior.

## 11. MATLAB local functions

These functions can stay at the bottom of the MATLAB script or be moved into separate files under `src/`.

```matlab
function model_gene_table = makeModelGeneTable(model, transcript_separator)
    model_gene_raw = string(model.genes(:));
    sep = regexptranslate('escape', transcript_separator);
    model_gene_clean = regexprep(cellstr(model_gene_raw), [sep '.*$'], '');
    model_gene_table = table(model_gene_clean, cellstr(model_gene_raw), ...
        'VariableNames', {'Gene_ID', 'Model_Gene_ID'});
end

function [target_data, rxn_expression, parsedGPR] = prepareTargetExpression( ...
    model, model_gene_table, expression_file, gene_col, expr_col)

    target_ref = readtable(expression_file, 'ReadVariableNames', true);
    expr_values = target_ref.(expr_col);
    if iscell(expr_values) || isstring(expr_values)
        expr_values = str2double(string(expr_values));
    end

    target_ref = table( ...
        cellstr(string(target_ref.(gene_col))), double(expr_values), ...
        'VariableNames', {'Gene_ID', 'Expression'});

    target_data = outerjoin(model_gene_table, target_ref, ...
        'Keys', 'Gene_ID', 'MergeKeys', true, 'Type', 'left');
    target_data.Expression(isnan(target_data.Expression)) = 0;

    expressionData = struct();
    expressionData.gene = target_data.Model_Gene_ID;
    expressionData.value = target_data.Expression * 2;

    [rxn_expression, parsedGPR] = mapExpressionToReactions(model, expressionData);
    rxn_expression(rxn_expression == -1) = 0;
end

function differ_genes = prepareDegTable(deg_file, gene_col, logfc_col, padj_col)
    deg_raw = readtable(deg_file, 'ReadVariableNames', true);

    padj_values = deg_raw.(padj_col);
    if iscell(padj_values) || isstring(padj_values)
        padj_values = str2double(string(padj_values));
    end

    logfc_values = deg_raw.(logfc_col);
    if iscell(logfc_values) || isstring(logfc_values)
        logfc_values = str2double(string(logfc_values));
    end

    differ_genes = table( ...
        cellstr(string(deg_raw.(gene_col))), double(logfc_values), double(padj_values), ...
        'VariableNames', {'gene', 'logFC', 'pval'});

    differ_genes = differ_genes(~isnan(differ_genes.logFC) & ~isnan(differ_genes.pval), :);
end

function [samples_full, sampleStats_full, Vref, rxnInactive] = expandSamplingToFullModel( ...
    model, tissueModel_ref, samples_raw)

    sampleStats_raw = calcSampleStats(samples_raw);
    n_rxns = numel(model.rxns);
    n_tissue_rxns = numel(tissueModel_ref.rxns);

    idx = zeros(n_tissue_rxns, 1);
    flip_sign = false(n_tissue_rxns, 1);

    for k = 1:n_tissue_rxns
        hit = find(strcmp(model.rxns, tissueModel_ref.rxns{k}), 1);
        if isempty(hit)
            hit = find(contains(model.rxns, tissueModel_ref.rxns{k}), 1);
            flip_sign(k) = true;
        end
        idx(k) = hit;
    end

    fields = fieldnames(sampleStats_raw);
    sampleStats_full = struct();
    for i = 1:numel(fields)
        aux = sampleStats_raw.(fields{i});
        aux = aux(:);
        aux(flip_sign) = -aux(flip_sign);
        tmp = zeros(n_rxns, 1);
        tmp(idx) = aux;
        sampleStats_full.(fields{i}) = tmp;
    end

    samples_adj = samples_raw;
    samples_adj(flip_sign, :) = -samples_adj(flip_sign, :);
    samples_full = zeros(n_rxns, size(samples_raw, 2));
    samples_full(idx, :) = samples_adj;

    Vref = sampleStats_full.mean;
    rxnInactive = setdiff((1:n_rxns)', idx);
end
```

## 12. MATLAB output files

For each replicate, the MATLAB pipeline writes:

| File pattern | Purpose |
|---|---|
| `results/MATLAB_outputs/source_vs_target_rep01_padj005_FC0_CHRR_rMTA.xlsx` | rMTA Excel output for replicate 1. |
| `results/MATLAB_outputs/source_vs_target_rep02_padj005_FC0_CHRR_rMTA.xlsx` | rMTA Excel output for replicate 2. |
| `results/MATLAB_outputs/Logs/source_vs_target_rep01_rxnFBS.mat` | Reaction-level DEG signal and active reaction count. |
| `results/MATLAB_outputs/Logs/source_vs_target_rep01_rMTA_workspace.mat` | Replicate-level MATLAB workspace/log variables. |
---

# ======================================================================
# Part III — rMTA result analysis
# ======================================================================

## 1. Purpose

The final R step analyzes the replicate-level rMTA outputs produced by the MATLAB pipeline.

It performs the following tasks:

1. reads and merges all rMTA replicate result files,
2. calculates rank stability statistics across replicates,
3. filters retained genes by rank variation and valid experiment count,
4. runs GO Biological Process and KEGG enrichment analysis with `clusterProfiler`,
5. saves analysis tables under `results/analysis_results/`,
6. saves enrichment plots under `plots/`.

The result-analysis script is:

```text
src/analyze_rMTA_results.R
```

## 2. Required R packages

```r
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
```

For mouse data, the default annotation package is:

```r
org.Mm.eg.db
```

For human data, replace it with:

```r
org.Hs.eg.db
```

and update the KEGG organism code from `mmu` to `hsa`.

## 3. User-defined parameters

Edit this block before running the result-analysis script.

```r
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
```

## 4. rMTA input files

This step reads the `n` rMTA result files generated by the MATLAB pipeline under:

```text
results/MATLAB_outputs/
```

Default file pattern:

```text
source_vs_target_rep01_padj005_FC0_CHRR_rMTA.xlsx
source_vs_target_rep02_padj005_FC0_CHRR_rMTA.xlsx
...
source_vs_target_rep10_padj005_FC0_CHRR_rMTA.xlsx
```

The script reads the rMTA result sheet:

```text
alpha = 0.66
```

The required columns are:

| Column | Meaning |
|---|---|
| `gene` | Gene ID or gene symbol returned by rMTA. |
| `rTS` | Final rMTA score used for ranking and filtering. |
| `GeneSymbol` | Optional gene symbol column, if already present in the rMTA output. |

If the rMTA `gene` column already contains gene symbols, keep:

```r
gene_symbol_column <- "gene"
```

If the rMTA output contains a separate `GeneSymbol` column, use:

```r
gene_symbol_column <- "GeneSymbol"
```

## 5. Merge rMTA replicate results

The script reads all matching rMTA Excel files, ranks genes in each replicate by descending `rTS`, and builds one merged table.

```r
rmta_files <- list.files(
  rmta_dir,
  pattern = rmta_file_pattern,
  full.names = TRUE
)

rmta_files <- sort(rmta_files)

if (length(rmta_files) == 0) {
  stop("No rMTA result files found in: ", rmta_dir)
}

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
      rTS = readr::parse_number(
        as.character(rTS),
        locale = readr::locale(decimal_mark = ".")
      )
    ) %>%
    filter(
      !is.na(gene),
      gene != "",
      !is.na(rTS)
    ) %>%
    arrange(desc(rTS))
  
  if (gene_symbol_column %in% colnames(rmta_res)) {
    rmta_res$GeneSymbol <- trimws(as.character(rmta_res[[gene_symbol_column]]))
  } else {
    rmta_res$GeneSymbol <- rmta_res$gene
  }
  
  rmta_res <- rmta_res %>%
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
```

The merged table is saved as:

```text
results/analysis_results/source_vs_target_rMTA_merged_rank_statistics.csv
```

The merged table contains:

| Column type | Meaning |
|---|---|
| `gene` | Gene ID from rMTA. |
| `GeneSymbol` | Gene symbol used for enrichment analysis. |
| `rTS1`, `rTS2`, ... | rTS values across rMTA replicates. |
| `rank1`, `rank2`, ... | Rank positions across rMTA replicates. |
| `valid_experiment_count` | Number of valid rMTA replicates for each gene. |
| `valid_median` | Median rank across valid replicates. |
| `valid_rank_range` | Difference between worst and best valid rank. |

## 6. Filter retained rMTA genes

The merged rMTA table is filtered using rank stability and valid experiment count.

```r
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
```

The filtered table is saved as:

```text
results/analysis_results/source_vs_target_rMTA_retained_genes_var10_minValid8.csv
```

A gene is retained if it satisfies:

```r
valid_experiment_count >= min_valid_experiments
valid_rank_range <= rank_variation_threshold
```

If `require_positive_rTS = TRUE`, all valid rTS values for that gene must also be positive.

## 7. Build retained-gene partitions

Retained genes are ordered by `valid_median`.

The ordered gene list is then split into cumulative partitions:

| Partition | Meaning |
|---:|---|
| `1` | All retained genes. |
| `0.75` | Top 75% of retained genes. |
| `0.50` | Top 50% of retained genes. |
| `0.25` | Top 25% of retained genes. |

```r
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
      partition = as.character(partition)
    )
  
  partition_tables[[as.character(partition)]] <- tmp_partition
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
```

The partitioned retained-gene table is saved as:

```text
results/analysis_results/source_vs_target_rMTA_retained_gene_partitions_var10_minValid8.csv
```

## 8. GO-BP and KEGG enrichment analysis

The script maps retained gene symbols to Entrez IDs and runs enrichment analysis with `clusterProfiler`.

```r
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

partition_gene_lists <- ora_genes_mapped %>%
  split(.$partition) %>%
  lapply(function(x) unique(as.character(x$ENTREZID)))

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
```

The enrichment results are saved as:

```text
results/analysis_results/source_vs_target_GO_BP_var10_minValid8.csv
results/analysis_results/source_vs_target_KEGG_var10_minValid8.csv
results/analysis_results/source_vs_target_ORA_input_genes_var10_minValid8.csv
results/analysis_results/source_vs_target_ORA_mapped_entrez_var10_minValid8.csv
```

## 9. Plot enrichment results

The script uses the included custom plotting functions to visualize GO-BP and KEGG enrichment results.

```r
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

get_partition_n_labels <- function(ora_genes) {
  
  partition_levels <- as.character(partitions)
  
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
    paste0(
      partition_n$partition,
      "\n(n = ",
      partition_n$n_genes,
      ")"
    ),
    partition_n$partition
  )
  
  partition_labels <- setNames(
    partition_levels,
    partition_levels
  )
  
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
  
  df <- df %>%
    arrange(p.adjust)
  
  keep_terms <- unique(df$Description)[
    1:min(top_n, length(unique(df$Description)))
  ]
  
  df %>%
    filter(Description %in% keep_terms)
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

plot_horizontal <- function(df, plot_title, partition_labels) {
  
  term_order <- get_term_order(df)
  
  df_plot <- df %>%
    mutate(
      x_term = factor(Description, levels = term_order),
      Cluster = factor(Cluster, levels = as.character(partitions))
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

plot_vertical <- function(df, plot_title, partition_labels) {
  
  term_order <- get_term_order(df)
  
  df_plot <- df %>%
    mutate(
      Cluster = factor(Cluster, levels = as.character(partitions)),
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
```

The GO-BP and KEGG plots are generated and saved as follows:

```r
partition_labels <- get_partition_n_labels(ora_genes)

go_plot_table <- read_and_tidy_enrichment(
  go_table,
  top_n = top_n_terms
)

kegg_plot_table <- read_and_tidy_enrichment(
  kegg_table,
  top_n = top_n_terms
)

if (nrow(go_plot_table) > 0) {
  
  p_go_h <- plot_horizontal(
    go_plot_table,
    paste0(
      comparison_name,
      " GO-BP enrichment | var=",
      rank_variation_threshold,
      ", minValid=",
      min_valid_experiments
    ),
    partition_labels
  )
  
  p_go_v <- plot_vertical(
    go_plot_table,
    paste0(
      comparison_name,
      " GO-BP enrichment | var=",
      rank_variation_threshold,
      ", minValid=",
      min_valid_experiments
    ),
    partition_labels
  )
  
  dims_go_h <- get_horizontal_dims(go_plot_table)
  dims_go_v <- get_vertical_dims(go_plot_table)
  
  ggsave(
    filename = file.path(
      plot_dir,
      paste0(
        comparison_name,
        "_GO_BP_dotplot_horizontal_var",
        rank_variation_threshold,
        "_minValid",
        min_valid_experiments,
        ".png"
      )
    ),
    plot = p_go_h,
    width = dims_go_h$width,
    height = dims_go_h$height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(
      plot_dir,
      paste0(
        comparison_name,
        "_GO_BP_dotplot_vertical_var",
        rank_variation_threshold,
        "_minValid",
        min_valid_experiments,
        ".png"
      )
    ),
    plot = p_go_v,
    width = dims_go_v$width,
    height = dims_go_v$height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}

if (nrow(kegg_plot_table) > 0) {
  
  p_kegg_h <- plot_horizontal(
    kegg_plot_table,
    paste0(
      comparison_name,
      " KEGG enrichment | var=",
      rank_variation_threshold,
      ", minValid=",
      min_valid_experiments
    ),
    partition_labels
  )
  
  p_kegg_v <- plot_vertical(
    kegg_plot_table,
    paste0(
      comparison_name,
      " KEGG enrichment | var=",
      rank_variation_threshold,
      ", minValid=",
      min_valid_experiments
    ),
    partition_labels
  )
  
  dims_kegg_h <- get_horizontal_dims(kegg_plot_table)
  dims_kegg_v <- get_vertical_dims(kegg_plot_table)
  
  ggsave(
    filename = file.path(
      plot_dir,
      paste0(
        comparison_name,
        "_KEGG_dotplot_horizontal_var",
        rank_variation_threshold,
        "_minValid",
        min_valid_experiments,
        ".png"
      )
    ),
    plot = p_kegg_h,
    width = dims_kegg_h$width,
    height = dims_kegg_h$height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(
      plot_dir,
      paste0(
        comparison_name,
        "_KEGG_dotplot_vertical_var",
        rank_variation_threshold,
        "_minValid",
        min_valid_experiments,
        ".png"
      )
    ),
    plot = p_kegg_v,
    width = dims_kegg_v$width,
    height = dims_kegg_v$height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}
```

Default plot outputs:

```text
plots/source_vs_target_GO_BP_dotplot_horizontal_var10_minValid8.png
plots/source_vs_target_GO_BP_dotplot_vertical_var10_minValid8.png
plots/source_vs_target_KEGG_dotplot_horizontal_var10_minValid8.png
plots/source_vs_target_KEGG_dotplot_vertical_var10_minValid8.png
```





