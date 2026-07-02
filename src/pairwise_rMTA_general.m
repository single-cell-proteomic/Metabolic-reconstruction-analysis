% Pairwise metabolic transformation analysis: MATLAB/rMTA general script
% This script continues from src/pairwise_prepare_inputs.R.
% Run the complete file from the MATLAB Editor; local helper functions are at the end.

% Direction rule:
% DEG input must be source state versus target state.
% Positive logFC = higher expression in source.
% Negative logFC = higher expression in target.

%% User-defined parameters
clear; clc;

model_file = 'data/Recon3D.mat';
target_expression_file = 'results/R_outputs/target_state_consensus_reference.tsv';
deg_file = 'results/R_outputs/DEG_source_vs_target_all_genes.csv';

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

%% Load COBRA model and prepare model gene IDs
% Gene_ID is the cleaned ID used to match R outputs.
% Model_Gene_ID is the original ID stored inside model.genes and used by COBRA.
% If model genes do not contain transcript-like suffixes after '.', IDs are unchanged.

changeCobraSolver(solver_name, 'all');
model = readCbModel(model_file);
model_gene_table = makeModelGeneTable(model, transcript_separator);

%% Load and map target-state consensus expression
% Missing model genes are assigned Expression = 0, meaning no expression evidence / moderate.
% Discrete -1/0/1 values are multiplied by 2 to match the +/-1 iMAT thresholds.

[target_data, target_rxn_expression, parsedGPR] = prepareTargetExpression( ...
    model, model_gene_table, target_expression_file, ...
    target_gene_id_column, target_expression_column);

fprintf('\nMapped target expression from gene level to reaction level.\n');
fprintf('Target genes with non-zero expression evidence: %d / %d\n', ...
    sum(target_data.Expression ~= 0), height(target_data));

%% Build target-state tissue-specific model with iMAT
% This model represents the desired/healthy/young target metabolic state.

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

%% Load and prepare source-vs-target DEG table
% rMTA helper functions expect a column named pval.
% This script uses adjusted p-values, so padj is copied to the internal pval field.

differ_genes = prepareDegTable( ...
    deg_file, deg_gene_id_column, deg_logfc_column, deg_padj_column);

fprintf('\nDEG table loaded: %d genes\n', height(differ_genes));
fprintf('Using adjusted p-value threshold: %.4g\n', padj_threshold);
fprintf('Using absolute logFC threshold: %.4g\n', fc_threshold);

%% Run rMTA replicates
% Each replicate repeats flux sampling on the same target tissue-specific model.
% Vref is the sampled mean flux vector of the target metabolic state.
% Check active_reaction_count after DEG-to-reaction mapping; very high values
% may make wTS/bTS less informative and may require stricter padj/logFC thresholds.

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

%% Local functions
function model_gene_table = makeModelGeneTable(model, transcript_separator)
    model_gene_raw = string(model.genes(:));
    sep = regexptranslate('escape', transcript_separator);
    model_gene_clean = regexprep(cellstr(model_gene_raw), [sep '.*$'], '');
    model_gene_table = table(model_gene_clean, cellstr(model_gene_raw), ...
        'VariableNames', {'Gene_ID', 'Model_Gene_ID'});
end

function [target_data, rxn_expression, parsedGPR] = prepareTargetExpression( ...
    model, model_gene_table, expression_file, gene_col, expr_col)

    target_ref = readtable(expression_file, 'ReadVariableNames', true, 'FileType', 'text');
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
