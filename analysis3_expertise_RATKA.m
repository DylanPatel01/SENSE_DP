%% analysis3_expertise_RATKA.m
% Analysis 3: Compare Expert vs Novice surgeons in RA-TKA
%
% Design:
%   - Independent groups: Expert vs Novice (RA-TKA only)
%   - Unit of analysis: surgeon (averaged across their procedures)
%   - Regions: 8 brain regions (canonical order)
%   - Conditions: shared_same, shared_diff, robotic_only
%   - Primary output: Hedges' g effect sizes (descriptive)
%   - Secondary output: Mann-Whitney U for ALL conditions where n>=2 per group
%                       (reported with explicit power caveat)
%
% Surgeon classification (RA-TKA):
%   Expert : 1 (n=2 proc), 2 (n=3 proc), 3 (n=1 proc), 4 (n=1 proc) — 4 surgeons
%   Novice : 5 (n=1 proc), 7 (n=1 proc)                               — 2 surgeons
%
% IMPORTANT POWER NOTE:
%   With n=4 experts vs n=2 novices, the minimum achievable Mann-Whitney p-value
%   is 0.133. Statistical significance is therefore not achievable regardless of
%   the true effect. Hedges' g is the primary outcome; p-values are reported
%   for completeness only and should not be used to draw inferential conclusions.
%
% Input:  master_grouped_region_HbDiff.csv (from synthesise_betas_grouped_SD.m)
% Output: analysis3_shared_conditions_HbDiff.csv   (shared_same, shared_diff)
%         analysis3_unique_conditions_HbDiff.csv    (robotic_only)
%         analysis3_descriptive_HbDiff.csv
%         (all output files share the same column schema including p-values)
%
% Dylan Patel | SENSE Project

clear all; close all; clc;

%% -------------------------------------------------------------------------
%  CONFIGURATION
%% -------------------------------------------------------------------------

signal     = 'HbDiff';
input_dir  = '/Volumes/Dylan SSD/DYLAN/Results/Synthesised_Grouped';
output_dir = '/Volumes/Dylan SSD/DYLAN/Results/Analysis';
master_file = fullfile(input_dir, sprintf('master_grouped_region_%s.csv', signal));

if ~isfolder(output_dir), mkdir(output_dir); end

alpha = 0.05;

% ── Expertise classification (surgeon ID -> group) ────────────────────────
% Update these lists if surgeon assignments change
expert_ids = ["1", "2", "3", "4"];
novice_ids = ["5", "7"];

% ── Canonical region order ────────────────────────────────────────────────
region_order = ["DLPFC_BA9_46", "OFC_BA11_47", "MPFC_BA9_10", "SFG_BA8_9", ...
                "PSC_BA1_2_3", "SPL_BA5_7", "SMG_BA40", "Angular_Gyrus_BA39"];

% ── Conditions to analyse ─────────────────────────────────────────────────
% shared_same & shared_diff: both groups may have data → test + effect size
% robotic_only: present in RA-TKA, compare expert vs novice descriptively
shared_conditions = {'shared_same', 'shared_diff'};
unique_conditions = {'robotic_only'};  % conventional_only not expected in RA-TKA
all_conditions    = [shared_conditions, unique_conditions];

%% -------------------------------------------------------------------------
%  LOAD DATA AND FILTER TO RA-TKA
%% -------------------------------------------------------------------------

T = readtable(master_file, 'TextType', 'string');

str_cols = {'Surgeon','Joint','Tech','Group_Key','Condition','Region','Region_Label'};
for s = 1:length(str_cols)
    if ismember(str_cols{s}, T.Properties.VariableNames)
        T.(str_cols{s}) = string(T.(str_cols{s}));
    end
end

% Filter: Knee joint, Robotic technique
T_tka = T(T.Joint == "K" & T.Tech == "R", :);
fprintf('RA-TKA rows (before expertise labelling): %d\n', height(T_tka));

% Assign expertise groups
T_tka.Expertise = repmat("Unknown", height(T_tka), 1);
T_tka.Expertise(ismember(T_tka.Surgeon, expert_ids)) = "Expert";
T_tka.Expertise(ismember(T_tka.Surgeon, novice_ids)) = "Novice";

% Report
exp_surgs = unique(T_tka.Surgeon(T_tka.Expertise == "Expert"));
nov_surgs = unique(T_tka.Surgeon(T_tka.Expertise == "Novice"));
unk_surgs = unique(T_tka.Surgeon(T_tka.Expertise == "Unknown"));

fprintf('\nExpert surgeons (%d): %s\n', length(exp_surgs), strjoin(exp_surgs, ', '));
fprintf('Novice surgeons (%d): %s\n', length(nov_surgs), strjoin(nov_surgs, ', '));
if ~isempty(unk_surgs)
    fprintf('WARNING — Unclassified surgeons (%d): %s\n', length(unk_surgs), strjoin(unk_surgs, ', '));
end

% Remove unclassified
T_tka = T_tka(T_tka.Expertise ~= "Unknown", :);
fprintf('RA-TKA rows after filtering: %d\n', height(T_tka));

% Build ordered region list from available data
available = unique(T_tka.Region);
regions   = region_order(ismember(region_order, available));
extra     = available(~ismember(available, region_order));
regions   = [cellstr(regions), cellstr(extra)'];
n_regions = length(regions);
fprintf('Regions: %d\n', n_regions);

%% -------------------------------------------------------------------------
%  DESCRIPTIVE STATISTICS (all conditions)
%% -------------------------------------------------------------------------

fprintf('\n--- Computing descriptive statistics ---\n');

desc_results = table();

for c = 1:length(all_conditions)
    cond = all_conditions{c};

    for r = 1:n_regions
        reg = regions{r};

        for grp = ["Expert", "Novice"]
            mask      = T_tka.Expertise == grp & ...
                        T_tka.Condition == cond & ...
                        T_tka.Region    == reg;
            betas     = T_tka.Beta(mask);
            reg_label = char(T_tka.Region_Label(find(T_tka.Region == reg, 1)));
            n_surg    = length(betas);

            if n_surg == 0, continue; end

            if n_surg > 1
                sd_val  = std(betas);
                sem_val = sd_val / sqrt(n_surg);
            else
                sd_val  = NaN;
                sem_val = NaN;
            end

            new_row = {cond, char(reg), reg_label, char(grp), ...
                mean(betas,'omitnan'), sd_val, sem_val, ...
                median(betas,'omitnan'), min(betas), max(betas), n_surg};
            desc_results = [desc_results; cell2table(new_row, 'VariableNames', ...
                {'Condition','Region','Region_Label','Expertise', ...
                 'Mean','SD','SEM','Median','Min','Max','n_surgeons'})];
        end
    end
end

%% -------------------------------------------------------------------------
%  HEDGES' G + MANN-WHITNEY U — ALL CONDITIONS
%  Mirrors the structure of analysis3_physiology_RATKA.m:
%    - Shared conditions (shared_same, shared_diff): Hedges' g + Mann-Whitney U
%    - Unique conditions (robotic_only):             Hedges' g only (p = NaN)
%  All conditions stored in a single unified results table.
%% -------------------------------------------------------------------------

fprintf('\n--- Computing Hedges'' g and Mann-Whitney U (all conditions) ---\n');

results = table();

for c = 1:length(all_conditions)
    cond = all_conditions{c};

    for r = 1:n_regions
        reg = regions{r};

        mask_exp  = T_tka.Expertise == "Expert" & T_tka.Condition == cond & T_tka.Region == reg;
        mask_nov  = T_tka.Expertise == "Novice" & T_tka.Condition == cond & T_tka.Region == reg;

        betas_exp = T_tka.Beta(mask_exp);
        betas_nov = T_tka.Beta(mask_nov);
        n_exp     = length(betas_exp);
        n_nov     = length(betas_nov);
        reg_label = char(T_tka.Region_Label(find(T_tka.Region == reg, 1)));

        if n_exp == 0 && n_nov == 0, continue; end

        diff_val = mean(betas_exp,'omitnan') - mean(betas_nov,'omitnan');

        % ── Hedges' g (always computed where possible) ───────────────────
        if n_exp >= 2 && n_nov >= 2
            df   = n_exp + n_nov - 2;
            sd_p = sqrt(((n_exp-1)*std(betas_exp)^2 + (n_nov-1)*std(betas_nov)^2) / df);
            g_note = sprintf('n_exp=%d, n_nov=%d', n_exp, n_nov);
        elseif n_exp >= 2 && n_nov == 1
            df   = n_exp - 1;
            sd_p = std(betas_exp);   % use expert SD as proxy
            g_note = sprintf('n_exp=%d, n_nov=1 — expert SD proxy', n_exp);
        elseif n_exp == 1 && n_nov >= 2
            df   = n_nov - 1;
            sd_p = std(betas_nov);   % use novice SD as proxy
            g_note = sprintf('n_exp=1, n_nov=%d — novice SD proxy', n_nov);
        else
            df   = NaN;
            sd_p = NaN;
            g_note = 'n=1 both groups — g not calculable';
        end

        if ~isnan(sd_p) && sd_p > 0
            d = diff_val / sd_p;
            if ~isnan(df) && df >= 1
                J = 1 - 3 / (4*df - 1);
            else
                J = 1;   % no correction possible
            end
            g = d * J;
        else
            d = NaN; J = NaN; g = NaN;
        end

        % ── Mann-Whitney U (all conditions, n >= 2 per group) ────────────
        if n_exp >= 2 && n_nov >= 2
            [p, ~, stats] = ranksum(betas_exp, betas_nov, 'alpha', alpha);
            U = stats.ranksum - n_exp*(n_exp+1)/2;
            interp_str = "Pending FDR";
        else
            p = NaN; U = NaN;
            interp_str = "g only (n<2 in one group)";
        end

        new_row = {cond, char(reg), reg_label, ...
            mean(betas_exp,'omitnan'), mean(betas_nov,'omitnan'), diff_val, ...
            U, p, NaN, d, g, abs(g), J, interp_str, n_exp, n_nov, g_note};
        results = [results; cell2table(new_row, 'VariableNames', ...
            {'Condition','Region','Region_Label', ...
             'Mean_Expert','Mean_Novice','Difference', ...
             'U_stat','p_uncorrected','p_fdr', ...
             'Cohens_d','Hedges_g','Abs_g','J_correction', ...
             'Interpretation','n_exp','n_nov','g_note'})];
    end
end

%% -------------------------------------------------------------------------
%  FDR CORRECTION (Benjamini-Hochberg, all conditions where p was computed)
%% -------------------------------------------------------------------------

valid_idx = ~isnan(results.p_uncorrected);
p_raw     = results.p_uncorrected(valid_idx);
n_valid   = sum(valid_idx);

if n_valid > 0
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
    results.p_fdr(valid_idx)  = p_fdr_reordered;
end

% Assign interpretations
for i = 1:height(results)
    if isnan(results.p_fdr(i)) || ~strcmp(char(results.Interpretation(i)), 'Pending FDR')
        % For g-only rows (unique conditions or n<2), assign effect label without significance
        if ~isnan(results.Abs_g(i))
            g   = results.Abs_g(i);
            dif = results.Difference(i);
            if g < 0.2,      eff = "Trivial";
            elseif g < 0.5,  eff = "Small";
            elseif g < 0.8,  eff = "Medium";
            else,            eff = "Large";
            end
            if dif > 0,      dir = "Expert > Novice";
            elseif dif < 0,  dir = "Novice > Expert";
            else,            dir = "No difference";
            end
            if strcmp(char(results.Interpretation(i)), 'Pending FDR')
                results.Interpretation(i) = "ns | " + eff + " | " + dir;
            else
                results.Interpretation(i) = results.Interpretation(i) + ...
                    " | " + eff + " | " + dir;
            end
        end
        continue;
    end

    p   = results.p_fdr(i);
    g   = results.Abs_g(i);
    dif = results.Difference(i);

    if p < 0.001,     sig = "***";
    elseif p < 0.01,  sig = "**";
    elseif p < alpha, sig = "*";
    else,             sig = "ns";
    end
    if g < 0.2,      eff = "Trivial";
    elseif g < 0.5,  eff = "Small";
    elseif g < 0.8,  eff = "Medium";
    else,            eff = "Large";
    end
    if dif > 0,      dir = "Expert > Novice";
    elseif dif < 0,  dir = "Novice > Expert";
    else,            dir = "No difference";
    end
    results.Interpretation(i) = sig + " | " + eff + " | " + dir;
end

%% -------------------------------------------------------------------------
%  SORT ALL OUTPUTS BY CANONICAL REGION ORDER
%% -------------------------------------------------------------------------

[~, ridx] = ismember(results.Region, cellstr(region_order));
ridx(ridx==0) = length(region_order)+1;
results = sortrows([results, table(ridx)], {'Condition','ridx'});
results.ridx = [];

[~, ridx] = ismember(desc_results.Region, cellstr(region_order));
ridx(ridx==0) = length(region_order)+1;
desc_results = sortrows([desc_results, table(ridx)], {'Condition','Expertise','ridx'});
desc_results.ridx = [];

%% -------------------------------------------------------------------------
%  SAVE RESULTS
%  - Shared conditions (shared_same, shared_diff): results rows with p-values
%  - Unique conditions (robotic_only):              results rows, p = NaN
%  Both subsets share the same column schema for consistent downstream use.
%% -------------------------------------------------------------------------

results_shared = results(ismember(results.Condition, shared_conditions), :);
results_unique = results(ismember(results.Condition, unique_conditions), :);

out_shared = fullfile(output_dir, sprintf('analysis3_shared_conditions_%s.csv', signal));
out_unique = fullfile(output_dir, sprintf('analysis3_unique_conditions_%s.csv', signal));
out_desc   = fullfile(output_dir, sprintf('analysis3_descriptive_%s.csv',       signal));

writetable(results_shared, out_shared);
writetable(results_unique, out_unique);
writetable(desc_results,   out_desc);

fprintf('\nSaved: %s\n', out_shared);
fprintf('Saved: %s\n', out_unique);
fprintf('Saved: %s\n', out_desc);

%% -------------------------------------------------------------------------
%  PRINT SUMMARY
%% -------------------------------------------------------------------------

fprintf('\n');
fprintf('========================================================================\n');
fprintf('  ANALYSIS 3: Expert vs Novice — RA-TKA\n');
fprintf('========================================================================\n');
fprintf('  Expert surgeons: %s (n=%d)\n', strjoin(exp_surgs, ', '), length(exp_surgs));
fprintf('  Novice surgeons: %s (n=%d)\n', strjoin(nov_surgs, ', '), length(nov_surgs));
fprintf('\n  *** POWER WARNING ***\n');
fprintf('  Minimum achievable p-value with n=%d vs n=%d is %.3f.\n', ...
    length(exp_surgs), length(nov_surgs), ...
    2 * min(1, (factorial(length(exp_surgs)) * factorial(length(nov_surgs))) / ...
    factorial(length(exp_surgs) + length(nov_surgs))));
fprintf('  Statistical significance is not achievable with this sample size.\n');
fprintf('  Hedges'' g is the PRIMARY outcome. p-values reported for completeness only.\n');
fprintf('------------------------------------------------------------------------\n');

% ── Descriptive stats ────────────────────────────────────────────────────
fprintf('\n  DESCRIPTIVE STATISTICS (Mean ± SD)\n');
for c = 1:length(all_conditions)
    cond = all_conditions{c};
    fprintf('\n  Condition: %s\n', upper(cond));
    fprintf('  %-40s %12s %12s %12s %12s\n', ...
        'Region', 'Expert Mean', 'Expert SD', 'Novice Mean', 'Novice SD');
    fprintf('  %s\n', repmat('-', 1, 92));

    for r = 1:n_regions
        reg      = regions{r};
        exp_mask = strcmp(desc_results.Condition, cond) & ...
                   strcmp(desc_results.Region, reg) & ...
                   strcmp(desc_results.Expertise, 'Expert');
        nov_mask = strcmp(desc_results.Condition, cond) & ...
                   strcmp(desc_results.Region, reg) & ...
                   strcmp(desc_results.Expertise, 'Novice');

        if ~any(exp_mask) && ~any(nov_mask), continue; end

        exp_mean = NaN; exp_sd = NaN; nov_mean = NaN; nov_sd = NaN;
        if any(exp_mask)
            exp_mean = desc_results.Mean(exp_mask);
            exp_sd   = desc_results.SD(exp_mask);
        end
        if any(nov_mask)
            nov_mean = desc_results.Mean(nov_mask);
            nov_sd   = desc_results.SD(nov_mask);
        end

        reg_label = char(T_tka.Region_Label(find(T_tka.Region == reg, 1)));
        fprintf('  %-40s %12.4e %12.4e %12.4e %12.4e\n', ...
            reg_label, exp_mean, exp_sd, nov_mean, nov_sd);
    end
end

% ── All conditions: effect sizes + test ──────────────────────────────────
fprintf('\n');
fprintf('========================================================================\n');
fprintf('  HEDGES'' g + MANN-WHITNEY U — ALL CONDITIONS\n');
fprintf('========================================================================\n');
fprintf('  Benchmarks: Trivial<0.2, Small 0.2-0.5, Medium 0.5-0.8, Large>=0.8\n');
fprintf('  FDR: Benjamini-Hochberg across ALL condition tests (where n>=2 per group)\n');
fprintf('  p = N/A only where one group has n<2\n');
fprintf('------------------------------------------------------------------------\n');

for c = 1:length(all_conditions)
    cond      = all_conditions{c};
    mask      = strcmp(results.Condition, cond);
    T_c       = results(mask, :);

    % Display in canonical region order
    [~, ridx] = ismember(T_c.Region, cellstr(region_order));
    ridx(ridx==0) = length(region_order)+1;
    [~, sidx] = sort(ridx);
    T_c = T_c(sidx, :);

    fprintf('\n  Condition: %s\n', upper(cond));

    fprintf('  %-40s %6s %6s %12s %12s %12s %7s %7s %8s %8s  %-30s  %s\n', ...
        'Region', 'n_Exp', 'n_Nov', 'Mean_Exp', 'Mean_Nov', 'Diff', ...
        'd', 'g', 'p_raw', 'p_fdr', 'Result', 'Note');
    fprintf('  %s\n', repmat('-', 1, 175));

    for i = 1:height(T_c)
        if isnan(T_c.p_uncorrected(i))
            p_raw_str = '     N/A ';
            p_fdr_str = '     N/A ';
        else
            p_raw_str = sprintf('%8.4f', T_c.p_uncorrected(i));
            p_fdr_str = sprintf('%8.4f', T_c.p_fdr(i));
        end
        fprintf('  %-40s %6d %6d %12.4e %12.4e %12.4e %7.3f %7.3f %s %s  %-30s  %s\n', ...
            T_c.Region_Label{i}, T_c.n_exp(i), T_c.n_nov(i), ...
            T_c.Mean_Expert(i), T_c.Mean_Novice(i), T_c.Difference(i), ...
            T_c.Cohens_d(i), T_c.Hedges_g(i), ...
            p_raw_str, p_fdr_str, ...
            char(T_c.Interpretation(i)), char(T_c.g_note(i)));
    end
end

% ── Top effects summary (all conditions) ─────────────────────────────────
fprintf('\n');
fprintf('========================================================================\n');
fprintf('  LARGEST EFFECT SIZES (|g| >= 0.5, all conditions)\n');
fprintf('========================================================================\n');
all_g = results(results.Abs_g >= 0.5, ...
    {'Condition','Region_Label','Mean_Expert','Mean_Novice','Hedges_g','Abs_g','Interpretation'});
if ~isempty(all_g)
    [~, sidx] = sort(all_g.Abs_g, 'descend');
    all_g = all_g(sidx, :);
    for i = 1:height(all_g)
        fprintf('  [%s] %-40s  g=%.3f  %s\n', ...
            all_g.Condition{i}, all_g.Region_Label{i}, ...
            all_g.Hedges_g(i), char(all_g.Interpretation(i)));
    end
else
    fprintf('  No effects >= 0.5 found.\n');
end

fprintf('\n========================================================================\n');
fprintf('  INTERPRETATION NOTES\n');
fprintf('========================================================================\n');
fprintf('  Unit of analysis: surgeon (betas averaged across procedures per surgeon).\n');
fprintf('  Expert n=%d, Novice n=%d — underpowered for inferential testing.\n', ...
    length(exp_surgs), length(nov_surgs));
fprintf('  Hedges'' g uses surgeon-level SD; J = 1 - 3/(4*df-1).\n');
fprintf('  Positive g/difference = Expert > Novice activation.\n');
fprintf('  Negative g/difference = Novice > Expert activation.\n');
fprintf('  p-values included for completeness but should not be used inferentially.\n');
fprintf('  Mann-Whitney U computed for all conditions where n>=2 per group.\n');
fprintf('  FDR (Benjamini-Hochberg) applied across all valid p-values.\n');
fprintf('\nDone.\n');