# Pairwise Metabolic Transformation Analysis Pipeline

This repository contains a three-step pipeline for pairwise metabolic transformation analysis between a **source** state and a **target** state.

The full detailed manual is available in:

```text
 detailed_manual.md
```

This README provides the minimal information needed to run the pipeline: required input files, parameters, commands, and expected outputs.

---

## Repository structure

All executable scripts are stored under `src/`. Run all commands from the repository root.

```text
.
├── README.md
├── detailed_manual.md
├── data/
│   ├── expression_counts.tsv
│   ├── metadata.csv
│   ├── gem_genes.csv
│   ├── mrna_genes.csv              # optional
│   ├── id_mapping.csv              # optional
│   └── Recon3D.mat                 # or another COBRA-readable GEM model
├── results/
│   ├── R_outputs/
│   ├── MATLAB_outputs/
│   │   └── Logs/
│   └── analysis_results/
├── plots/
└── src/
    ├── pairwise_prepare_inputs.R
    ├── pairwise_rMTA_general.m
    └── analyze_rMTA_results.R
```

---

## Dependencies

### R

The R preparation script requires:

```r
library(DESeq2)
library(dplyr)
library(tibble)
```

The R result-analysis script requires:

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

For mouse data, the default annotation package is `org.Mm.eg.db`. For human data, use `org.Hs.eg.db` and set the KEGG organism code to `hsa`.

### MATLAB / COBRA

The MATLAB/rMTA script requires:

- MATLAB,
- COBRA Toolbox,
- a compatible MILP solver such as Gurobi,
- a COBRA-readable GEM model,
- rMTA helper functions available on the MATLAB path:
  - `diffexprs2rxnFBS`,
  - `calculateEPSILON`,
  - `rMTA`,
  - `rMTAsaveInExcel`.

Before running the MATLAB script, initialize COBRA Toolbox and check the solver:

```matlab
initCobraToolbox(false)
changeCobraSolver('gurobi', 'all')
```

---

# Step 1 — Prepare R-side rMTA input files

## Script

```text
src/pairwise_prepare_inputs.R
```

## Run command

```bash
Rscript src/pairwise_prepare_inputs.R
```

## Input files

| File | Description |
|---|---|
| `data/expression_counts.tsv` | Raw RNA-seq count matrix with genes in rows and samples in columns. |
| `data/metadata.csv` | Sample metadata table with at least `sample_id` and `state` columns. |
| `data/gem_genes.csv` | Gene universe of the GEM model. The `Gene_ID` column must match the model gene IDs. |
| `data/mrna_genes.csv` | Optional mRNA/protein-coding gene list used when `use_mrna_filter = TRUE`. |
| `data/id_mapping.csv` | Optional ID conversion table used when `use_id_conversion = TRUE`. |

The metadata `state` column must contain the two states used as `source_label` and `target_label`. Differential expression is always computed as **source versus target**.

## Parameters

| Parameter | Default | Possible / recommended values | Meaning |
|---|---:|---|---|
| `expression_file` | `data/expression_counts.tsv` | File path | Raw count matrix. |
| `metadata_file` | `data/metadata.csv` | File path | Sample metadata. |
| `gem_gene_file` | `data/gem_genes.csv` | File path | GEM-compatible gene list. |
| `use_mrna_filter` | `TRUE` | `TRUE` / `FALSE` | Whether to keep only mRNA/protein-coding genes before analysis. |
| `mrna_gene_file` | `data/mrna_genes.csv` | File path | Used only if `use_mrna_filter = TRUE`. |
| `use_id_conversion` | `TRUE` | `TRUE` / `FALSE` | Whether to convert expression IDs to GEM-compatible IDs. |
| `id_mapping_file` | `data/id_mapping.csv` | File path | Used only if `use_id_conversion = TRUE`. |
| `source_label` | `source` | Any value present in `metadata$state` | Initial/disease/old/patient state. |
| `target_label` | `target` | Any value present in `metadata$state` | Desired/healthy/young/reference state. |
| `pval_cutoff` | `0.05` | `0–1` | Raw p-value cutoff for exporting the GEM-filtered DEG table. |
| `sd_multiplier` | `0.6` | Positive numeric; usually `0.4–1.0` | Controls low/high expression discretization for target samples. |
| `consensus_threshold` | `0.8` | `0.5–1.0` | Minimum fraction of target samples that must agree for consensus expression. |
| `out_dir` | `results/R_outputs` | Folder path | Output directory for R-side files. |

## Output files

| File | Description |
|---|---|
| `results/R_outputs/DEG_source_vs_target_all_genes.csv` | Full source-vs-target DESeq2 result table. |
| `results/R_outputs/DEG_source_vs_target_GEM_genes_all.csv` | DESeq2 result table restricted to GEM genes. |
| `results/R_outputs/DEG_source_vs_target_GEM_genes_pval_0.05.csv` | GEM-filtered DEG table used by MATLAB/rMTA. |
| `results/R_outputs/normalized_expression_all_genes_VST.tsv` | VST-normalized expression matrix for all retained genes. |
| `results/R_outputs/normalized_expression_GEM_genes_VST.tsv` | VST-normalized expression matrix restricted to GEM genes. |
| `results/R_outputs/target_state_discretized_expression.tsv` | Discretized target-sample expression values. |
| `results/R_outputs/target_state_consensus_reference.tsv` | Target-state consensus expression reference used by iMAT. |
| `results/R_outputs/R_preparation_summary.csv` | Summary of sample counts, gene counts, and consensus states. |

The two main files required by Step 2 are `DEG_source_vs_target_GEM_genes_pval_0.05.csv` and `target_state_consensus_reference.tsv`.

---

# Step 2 — Run MATLAB/rMTA transformation analysis

## Script

```text
src/pairwise_rMTA_general.m
```

## Run command

From MATLAB, with the repository root as the working directory:

```matlab
run('src/pairwise_rMTA_general.m')
```

## Input files

| File | Description |
|---|---|
| `data/Recon3D.mat` | COBRA-readable GEM model. This can be replaced by another compatible model file. |
| `results/R_outputs/target_state_consensus_reference.tsv` | Target-state expression reference produced by Step 1. |
| `results/R_outputs/DEG_source_vs_target_GEM_genes_pval_0.05.csv` | Source-vs-target DEG table produced by Step 1. |

The DEG file must represent **source versus target**. Positive `logFC` means higher expression in source; negative `logFC` means higher expression in target.

## Parameters

| Parameter | Default | Possible / recommended values | Meaning |
|---|---:|---|---|
| `model_file` | `data/Recon3D.mat` | COBRA-readable model path | GEM model used by COBRA/rMTA. |
| `target_expression_file` | `results/R_outputs/target_state_consensus_reference.tsv` | File path | Target-state consensus expression reference. |
| `deg_file` | `results/R_outputs/DEG_source_vs_target_GEM_genes_pval_0.05.csv` | File path | DEG table used for rMTA. |
| `comparison_name` | `source_vs_target` | Text label | Prefix used in output file names. |
| `out_dir` | `results/MATLAB_outputs` | Folder path | Main MATLAB output folder. |
| `log_dir` | `results/MATLAB_outputs/Logs` | Folder path | MATLAB log/workspace output folder. |
| `target_gene_id_column` | `Gene_ID` | Column name | Gene ID column in the target expression file. |
| `target_expression_column` | `Expression` | Column name | Discrete expression column in the target expression file. |
| `deg_gene_id_column` | `Gene_ID` | Column name | Gene ID column in the DEG file. |
| `deg_logfc_column` | `logFC` | Column name | Source-vs-target log fold-change column. |
| `deg_padj_column` | `padj` | Column name | Adjusted p-value column. Internally copied to `pval` for rMTA compatibility. |
| `fc_threshold` | `0` | Numeric `>= 0` | Absolute logFC cutoff used in DEG-to-reaction mapping. `0` means no fold-change filter. |
| `padj_threshold` | `0.05` | `0–1` | Adjusted p-value cutoff used in DEG-to-reaction mapping. |
| `n_replicates` | `10` | Positive integer | Number of rMTA sampling/rMTA replicates. |
| `solver_name` | `gurobi` | COBRA-compatible solver | Optimization solver. |
| `transcript_separator` | `.` | Character | Removes transcript-like suffixes from model gene IDs if needed. |
| `sampling_method` | `CHRR` | COBRA sampling method | Sampling method used by `sampleCbModel`. |
| `alpha_values` | `0.66` | Numeric scalar or vector | rMTA alpha value(s). |
| `imat_time_limit` | `60` | Positive numeric | Time limit for iMAT in seconds. |
| `rmta_time_limit` | `10` | Positive numeric | Time limit for rMTA optimization in seconds. |
| `n_workers` | `4` | Positive integer | Number of workers used by COBRA/rMTA where supported. |
| `n_threads` | `4` | Positive integer | Number of solver/sampling threads. |
| `print_level` | `1` | Usually `0–3` | Verbosity level. |
| `sampling_options.nPointsReturned` | `2000` | Positive integer | Number of sampled flux points. |
| `sampling_options.nStepsPerPoint` | `200` | Positive integer | CHRR steps between returned points. |
| `sampling_options.nThreads` | `n_threads` | Positive integer | Sampling thread count. |

## Output files

| File | Description |
|---|---|
| `results/MATLAB_outputs/source_vs_target_rep01_padj005_FC0_CHRR_rMTA.xlsx` | rMTA Excel result for replicate 1. Additional replicates follow the same naming pattern. |
| `results/MATLAB_outputs/Logs/source_vs_target_target_tissueModel_ref.mat` | Target-state tissue-specific model generated by iMAT. |
| `results/MATLAB_outputs/Logs/source_vs_target_rep01_rxnFBS.mat` | Reaction-level DEG signal and active reaction count for replicate 1. |
| `results/MATLAB_outputs/Logs/source_vs_target_rep01_rMTA_workspace.mat` | Replicate-level rMTA workspace variables. |

The Excel files in `results/MATLAB_outputs/` are the main input for Step 3.

---

# Step 3 — Analyze rMTA results and run enrichment analysis

## Script

```text
src/analyze_rMTA_results.R
```

## Run command

```bash
Rscript src/analyze_rMTA_results.R
```

## Input files

| File | Description |
|---|---|
| `results/MATLAB_outputs/source_vs_target_rep*_padj005_FC0_CHRR_rMTA.xlsx` | Replicate-level rMTA Excel files generated by Step 2. |

The script reads the rMTA result sheet named `alpha = 0.66`. Each rMTA result table must contain at least `gene` and `rTS` columns. If available, `GeneSymbol` can be used for enrichment analysis.

## Parameters

| Parameter | Default | Possible / recommended values | Meaning |
|---|---:|---|---|
| `comparison_name` | `source_vs_target` | Text label | Prefix used to identify input/output files. |
| `rmta_dir` | `results/MATLAB_outputs` | Folder path | Folder containing rMTA Excel files. |
| `analysis_dir` | `results/analysis_results` | Folder path | Folder for merged tables and enrichment results. |
| `plot_dir` | `plots` | Folder path | Folder for enrichment plots. |
| `rmta_sheet` | `alpha = 0.66` | Excel sheet name | Sheet to read from each rMTA Excel file. |
| `rmta_file_pattern` | Derived from `comparison_name` | Regular expression | Pattern used to select rMTA replicate files. |
| `rank_variation_threshold` | `10` | Positive integer | Maximum allowed rank range across valid replicates. |
| `min_valid_experiments` | `8` | Integer from `1` to `n_replicates` | Minimum number of valid rMTA replicates required for a gene. |
| `require_positive_rTS` | `TRUE` | `TRUE` / `FALSE` | Whether retained genes must have positive valid rTS values. |
| `invalid_rTS_value` | `65535` | Numeric sentinel value | rMTA invalid-score value excluded from valid statistics. |
| `partitions` | `c(1, 0.75, 0.50, 0.25)` | Fractions in `(0, 1]` | Cumulative retained-gene subsets used for enrichment analysis. |
| `gene_symbol_column` | `gene` | Column name | Column used as gene symbol input for annotation mapping. |
| `org_db` | `org.Mm.eg.db` | Annotation database object | Organism annotation database. Use `org.Hs.eg.db` for human. |
| `org_keytype` | `SYMBOL` | Valid key type in `org_db` | Key type used for mapping to Entrez IDs. |
| `kegg_organism` | `mmu` | KEGG organism code | Use `mmu` for mouse and `hsa` for human. |
| `top_n_terms` | `50` | Positive integer | Maximum number of enriched terms shown in plots. |

## Output files

| File | Description |
|---|---|
| `results/analysis_results/source_vs_target_rMTA_merged_rank_statistics.csv` | Merged rMTA table containing rTS values, ranks, and rank-stability statistics across replicates. |
| `results/analysis_results/source_vs_target_rMTA_retained_genes_var10_minValid8.csv` | Filtered candidate genes after rank variation and valid-experiment filtering. |
| `results/analysis_results/source_vs_target_rMTA_retained_gene_partitions_var10_minValid8.csv` | Retained candidate genes split into cumulative partitions for enrichment analysis. |
| `results/analysis_results/source_vs_target_GO_BP_var10_minValid8.csv` | GO Biological Process enrichment results. |
| `results/analysis_results/source_vs_target_KEGG_var10_minValid8.csv` | KEGG enrichment results. |
| `results/analysis_results/source_vs_target_ORA_input_genes_var10_minValid8.csv` | Gene symbols used as enrichment input. |
| `results/analysis_results/source_vs_target_ORA_mapped_entrez_var10_minValid8.csv` | Gene symbol to Entrez ID mapping used for enrichment. |
| `plots/source_vs_target_GO_BP_dotplot_horizontal_var10_minValid8.png` | GO-BP enrichment dotplot. |
| `plots/source_vs_target_GO_BP_dotplot_vertical_var10_minValid8.png` | GO-BP enrichment dotplot. |
| `plots/source_vs_target_KEGG_dotplot_horizontal_var10_minValid8.png` | KEGG enrichment dotplot. |
| `plots/source_vs_target_KEGG_dotplot_vertical_var10_minValid8.png` | KEGG enrichment dotplot. |

---

## Complete run order

Run the three scripts in this order from the repository root:

```bash
Rscript src/pairwise_prepare_inputs.R
```

```matlab
run('src/pairwise_rMTA_general.m')
```

```bash
Rscript src/analyze_rMTA_results.R
```
