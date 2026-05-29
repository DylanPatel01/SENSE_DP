%% analysis1_technique_THA.m
% Analysis 1: Compare RA-THA vs C-THA across brain regions and stage conditions
%
% Design:
%   - Independent groups: RA-THA vs C-THA
%   - Unit of analysis: surgeon (averaged across their procedures)
%   - Regions: 8 brain regions
%   - Conditions: shared_same, shared_diff, robotic_only, conventional_only
%   - Test: Mann-Whitney U (non-parametric, small n)
%   - Correction: Benjamini-Hochberg FDR across all region x condition tests
%   - Effect size: Hedges' g (bias-corrected Cohen's d, appropriate for small n)
%
% Hedges' g formula:
%   g = d * J   where J = 1 - 3 / (4*df - 1)   and   df = n1 + n2 - 2
%   d = (mean1 - mean2) / SD_pooled
%   SD_pooled = sqrt(((n1-1)*SD1^2 + (n2-1)*SD2^2) / df)
%
% Benchmarks (Cohen 1988): Trivial<0.2, Small 0.2-0.5, Medium 0.5-0.8, Large>=0.8
%
% Input:  master_grouped_region_HbDiff.csv (from synthesise_betas_grouped_SD.m)
% Output: analysis1_shared_conditions_HbDiff.csv
%         analysis1_unique_conditions_HbDiff.csv
%
% Dylan Patel | SENSE Project

clear all; close all; clc;

%% -------------------------------------------------------------------------
%  CONFIGURATION
%% -------------------------------------------------------------------------

signal      = 'HbDiff';
input_dir   = '/Volumes/Dylan SSD/DYLAN/Results/Synthesised_Grouped';
output_dir  = '/Volumes/Dylan SSD/DYLAN/Results/Analysis';
master_file = fullfile(input_dir, sprintf('master_grouped_region_%s.csv', signal));

if ~isfolder(output_dir), mkdir(output_dir); end

alpha             = 0.05;
shared_conditions = {'shared_same', 'shared_diff'};
unique_conditions = {'robotic_only', 'conventional_only'};

%% -------------------------------------------------------------------------
%  HEDGES' G FUNCTION
%% -------------------------------------------------------------------------
% Inputs:  x1, x2 — vectors of observations for each group
% Outputs: g       — Hedges' g effect size
%          d       — Cohen's d (uncorrected)
%          sd_pool — pooled SD

hedges_g = @(x1, x2) deal( ...
    ((mean(x1) - mean(x2)) / sqrt(((length(x1)-1)*std(x1)^2 + (length(x2)-1)*std(x2)^2) / (length(x1)+length(x2)-2))) ...
    * (1 - 3 / (4*(length(x1)+length(x2)-2) - 1)), ...
    (mean(x1) - mean(x2)) / sqrt(((length(x1)-1)*std(x1)^2 + (length(x2)-1)*std(x2)^2) / (length(x1)+length(x2)-2)), ...
    sqrt(((length(x1)-1)*std(x1)^2 + (length(x2)-1)*std(x2)^2) / (length(x1)+length(x2)-2)) ...
);

%% -------------------------------------------------------------------------
%  LOAD AND FILTER TO THA ONLY
%% -------------------------------------------------------------------------

T = readtable(master_file, 'TextType', 'string');

str_cols = {'Surgeon','Joint','Tech','Group_Key','Condition','Region','Region_Label'};
for s = 1:length(str_cols)
    if ismember(str_cols{s}, T.Properties.VariableNames)
        T.(str_cols{s}) = string(T.(str_cols{s}));
    end
end

T_tha = T(T.Joint == "H", :);
fprintf('THA rows: %d\n', height(T_tha));

surgeons_rob  = unique(T_tha.Surgeon(T_tha.Tech == "R"));
surgeons_conv = unique(T_tha.Surgeon(T_tha.Tech == "C"));
fprintf('\nRA-THA surgeons (%d): %s\n', length(surgeons_rob),  strjoin(surgeons_rob,  ', '));
fprintf('C-THA  surgeons (%d): %s\n', length(surgeons_conv), strjoin(surgeons_conv, ', '));

% Canonical region order (matches Table 3)
region_order = ["DLPFC_BA9_46", "OFC_BA11_47", "MPFC_BA9_10", "SFG_BA8_9", ...
                "PSC_BA1_2_3", "SPL_BA5_7", "SMG_BA40", "Angular_Gyrus_BA39"];

% Use canonical order, falling back to any regions present in data
available = unique(T_tha.Region);
regions   = region_order(ismember(region_order, available));
% Append any unexpected regions not in canonical list
extra     = available(~ismember(available, region_order));
regions   = [regions, extra'];
n_regions = length(regions);
fprintf('\nRegions: %d\n', n_regions);

%% -------------------------------------------------------------------------
%  MANN-WHITNEY U + HEDGES' G: SHARED CONDITIONS
%% -------------------------------------------------------------------------

fprintf('\n--- Running Mann-Whitney U tests with Hedges'' g (shared conditions) ---\n');

results_shared = table();

for c = 1:length(shared_conditions)
    cond = shared_conditions{c};

    for r = 1:n_regions
        reg = regions(r);

        mask_rob  = T_tha.Tech == "R" & T_tha.Condition == cond & T_tha.Region == reg;
        mask_conv = T_tha.Tech == "C" & T_tha.Condition == cond & T_tha.Region == reg;

        betas_rob  = T_tha.Beta(mask_rob);
        betas_conv = T_tha.Beta(mask_conv);
        n_rob      = length(betas_rob);
        n_conv     = length(betas_conv);
        reg_label  = char(T_tha.Region_Label(find(T_tha.Region == reg, 1)));

        if n_rob < 2 || n_conv < 2
            fprintf('  WARNING: Insufficient data for %s / %s (n_rob=%d, n_conv=%d) — skipping\n', ...
                cond, char(reg), n_rob, n_conv);
            new_row = {cond, char(reg), reg_label, ...
                mean(betas_rob,'omitnan'), mean(betas_conv,'omitnan'), ...
                mean(betas_rob,'omitnan') - mean(betas_conv,'omitnan'), ...
                NaN, NaN, NaN, NaN, NaN, NaN, NaN, "Insufficient data", n_rob, n_conv};
            results_shared = [results_shared; cell2table(new_row, 'VariableNames', ...
                {'Condition','Region','Region_Label', ...
                 'Mean_Robotic','Mean_Conventional','Difference', ...
                 'U_stat','p_uncorrected','p_fdr', ...
                 'Cohens_d','Hedges_g','Abs_g','J_correction', ...
                 'Interpretation','n_rob','n_conv'})];
            continue;
        end

        % Mann-Whitney U
        [p, ~, stats] = ranksum(betas_rob, betas_conv, 'alpha', alpha);
        U = stats.ranksum - n_rob*(n_rob+1)/2;

        % Hedges' g
        df     = n_rob + n_conv - 2;
        sd_p   = sqrt(((n_rob-1)*std(betas_rob)^2 + (n_conv-1)*std(betas_conv)^2) / df);
        d      = (mean(betas_rob) - mean(betas_conv)) / sd_p;
        J      = 1 - 3 / (4*df - 1);   % bias correction factor
        g      = d * J;

        difference = mean(betas_rob) - mean(betas_conv);

        new_row = {cond, char(reg), reg_label, ...
            mean(betas_rob), mean(betas_conv), difference, ...
            U, p, NaN, d, g, abs(g), J, ...
            "Pending FDR", n_rob, n_conv};
        results_shared = [results_shared; cell2table(new_row, 'VariableNames', ...
            {'Condition','Region','Region_Label', ...
             'Mean_Robotic','Mean_Conventional','Difference', ...
             'U_stat','p_uncorrected','p_fdr', ...
             'Cohens_d','Hedges_g','Abs_g','J_correction', ...
             'Interpretation','n_rob','n_conv'})];
    end
end

%% -------------------------------------------------------------------------
%  FDR CORRECTION (Benjamini-Hochberg)
%% -------------------------------------------------------------------------

valid_idx = ~isnan(results_shared.p_uncorrected);
p_raw     = results_shared.p_uncorrected(valid_idx);
n_valid   = sum(valid_idx);

[p_sorted, sort_idx] = sort(p_raw);
fdr_thresh = (1:n_valid)' / n_valid * alpha;
below      = p_sorted <= fdr_thresh;

if any(below)
    k_max = find(below, 1, 'last');
    p_fdr = min(p_sorted(k_max) * n_valid ./ (1:n_valid)', 1);
else
    p_fdr = ones(n_valid, 1);
end

p_fdr_reordered           = NaN(n_valid, 1);
p_fdr_reordered(sort_idx) = p_fdr;
results_shared.p_fdr(valid_idx) = p_fdr_reordered;

% Assign interpretations
for i = 1:height(results_shared)
    if isnan(results_shared.p_fdr(i)), continue; end

    p   = results_shared.p_fdr(i);
    g   = results_shared.Abs_g(i);
    dif = results_shared.Difference(i);

    if p < 0.001,      sig = "***";
    elseif p < 0.01,   sig = "**";
    elseif p < alpha,  sig = "*";
    else,              sig = "ns";
    end

    % Cohen benchmarks for g: Trivial<0.2, Small 0.2-0.5, Medium 0.5-0.8, Large>=0.8
    if g < 0.2,       eff = "Trivial";
    elseif g < 0.5,   eff = "Small";
    elseif g < 0.8,   eff = "Medium";
    else,             eff = "Large";
    end

    if dif > 0,       dir = "Robotic > Conventional";
    elseif dif < 0,   dir = "Conventional > Robotic";
    else,             dir = "No difference";
    end

    results_shared.Interpretation(i) = sig + " | " + eff + " | " + dir;
end

%% -------------------------------------------------------------------------
%  DESCRIPTIVE STATISTICS: UNIQUE CONDITIONS
%% -------------------------------------------------------------------------

fprintf('\n--- Descriptive stats for unique conditions ---\n');

results_unique = table();

for c = 1:length(unique_conditions)
    cond     = unique_conditions{c};
    tech_grp = "R";
    if cond == "conventional_only", tech_grp = "C"; end

    for r = 1:n_regions
        reg       = regions(r);
        mask      = T_tha.Tech == tech_grp & T_tha.Condition == cond & T_tha.Region == reg;
        betas     = T_tha.Beta(mask);
        reg_label = char(T_tha.Region_Label(find(T_tha.Region == reg, 1)));
        n_surg    = length(betas);

        if n_surg == 0, continue; end

        new_row = {cond, char(reg), reg_label, ...
            mean(betas,'omitnan'), std(betas,'omitnan'), ...
            median(betas,'omitnan'), min(betas), max(betas), n_surg};
        results_unique = [results_unique; cell2table(new_row, 'VariableNames', ...
            {'Condition','Region','Region_Label', ...
             'Mean','SD','Median','Min','Max','n_surgeons'})];
    end
end

%% -------------------------------------------------------------------------
%  SAVE RESULTS
%% -------------------------------------------------------------------------

% Sort results by canonical region order before saving
[~, reg_sort_idx] = ismember(results_shared.Region, cellstr(region_order));
reg_sort_idx(reg_sort_idx == 0) = length(region_order) + 1;
[~, final_sort]   = sortrows([reg_sort_idx, (1:height(results_shared))']);
results_shared    = results_shared(final_sort, :);

[~, reg_sort_idx] = ismember(results_unique.Region, cellstr(region_order));
reg_sort_idx(reg_sort_idx == 0) = length(region_order) + 1;
[~, final_sort]   = sortrows([reg_sort_idx, (1:height(results_unique))']);
results_unique    = results_unique(final_sort, :);

out_shared = fullfile(output_dir, sprintf('analysis1_shared_conditions_%s.csv', signal));
out_unique = fullfile(output_dir, sprintf('analysis1_unique_conditions_%s.csv', signal));
writetable(results_shared, out_shared);
writetable(results_unique, out_unique);
fprintf('\nSaved: %s\n', out_shared);
fprintf('Saved: %s\n', out_unique);

%% -------------------------------------------------------------------------
%  PRINT SUMMARY
%% -------------------------------------------------------------------------

fprintf('\n');
fprintf('========================================================================\n');
fprintf('  ANALYSIS 1: RA-THA vs C-THA — SHARED CONDITIONS\n');
fprintf('========================================================================\n');
fprintf('  FDR: Benjamini-Hochberg (alpha=%.2f)  |  Effect: Hedges'' g\n', alpha);
fprintf('  Benchmarks: Trivial<0.2, Small 0.2-0.5, Medium 0.5-0.8, Large>=0.8\n');
fprintf('  Significance: * p<.05  ** p<.01  *** p<.001  ns = not significant\n');
fprintf('  J correction factor applied (df = n_rob + n_conv - 2 = %d)\n', ...
    results_shared.n_rob(1) + results_shared.n_conv(1) - 2);
fprintf('------------------------------------------------------------------------\n');

for c = 1:length(shared_conditions)
    cond = shared_conditions{c};
    mask = strcmp(results_shared.Condition, cond);
    T_c  = results_shared(mask, :);

    % Sort by raw p-value
    [~, sidx] = sort(T_c.p_uncorrected);
    T_c = T_c(sidx, :);

    fprintf('\n  Condition: %s\n', upper(cond));
    fprintf('  %-40s %6s %6s %12s %12s %12s %8s %8s %7s %7s  %s\n', ...
        'Region', 'n_Rob', 'n_Con', 'Mean_Rob', 'Mean_Con', 'Diff', ...
        'p_raw', 'p_fdr', 'd', 'g', 'Result');
    fprintf('  %s\n', repmat('-', 1, 145));

    for i = 1:height(T_c)
        if T_c.p_uncorrected(i) < 0.05
            flag = ' <-- p<.05 uncorrected';
        elseif T_c.p_uncorrected(i) < 0.10
            flag = ' <-- trend (p<.10)';
        else
            flag = '';
        end
        fprintf('  %-40s %6d %6d %12.4e %12.4e %12.4e %8.4f %8.4f %7.3f %7.3f  %s%s\n', ...
            T_c.Region_Label{i}, T_c.n_rob(i), T_c.n_conv(i), ...
            T_c.Mean_Robotic(i), T_c.Mean_Conventional(i), T_c.Difference(i), ...
            T_c.p_uncorrected(i), T_c.p_fdr(i), ...
            T_c.Cohens_d(i), T_c.Hedges_g(i), ...
            char(T_c.Interpretation(i)), flag);
    end
end

fprintf('\n');
fprintf('========================================================================\n');
fprintf('  UNIQUE CONDITIONS (descriptive only)\n');
fprintf('========================================================================\n');
fprintf('  %-20s %-40s %6s %12s %12s %12s\n', ...
    'Condition', 'Region', 'n', 'Mean', 'SD', 'Median');
fprintf('  %s\n', repmat('-', 1, 100));
for i = 1:height(results_unique)
    fprintf('  %-20s %-40s %6d %12.4e %12.4e %12.4e\n', ...
        results_unique.Condition{i}, results_unique.Region_Label{i}, ...
        results_unique.n_surgeons(i), results_unique.Mean(i), ...
        results_unique.SD(i), results_unique.Median(i));
end

fprintf('\n');
fprintf('========================================================================\n');
fprintf('  SIGNIFICANT FINDINGS AFTER FDR CORRECTION\n');
fprintf('========================================================================\n');
sig_mask = results_shared.p_fdr < alpha;
if any(sig_mask)
    T_sig = results_shared(sig_mask, :);
    for i = 1:height(T_sig)
        fprintf('  [%s] %s: p_fdr=%.4f, g=%.3f (d=%.3f, J=%.3f) — %s\n', ...
            T_sig.Condition{i}, T_sig.Region_Label{i}, T_sig.p_fdr(i), ...
            T_sig.Hedges_g(i), T_sig.Cohens_d(i), T_sig.J_correction(i), ...
            char(T_sig.Interpretation(i)));
    end
else
    fprintf('  No significant differences after FDR correction.\n');
end

fprintf('\n========================================================================\n');
fprintf('  NOTE ON INTERPRETATION\n');
fprintf('========================================================================\n');
fprintf('  Unit of analysis: surgeon (betas averaged across procedures per surgeon).\n');
fprintf('  Hedges'' g used instead of Cohen''s d — bias-corrected for small n.\n');
fprintf('  J correction factor = 1 - 3/(4*df-1) where df = n1+n2-2.\n');
fprintf('  With n=4 per group, J = %.4f (%.1f%% downward correction on d).\n', ...
    1 - 3/(4*6-1), (1-(1-3/(4*6-1)))*100);
fprintf('  Mann-Whitney U used for p-values; Hedges'' g for effect magnitude.\n');
fprintf('\nDone.\n');
