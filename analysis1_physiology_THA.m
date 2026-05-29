%% analysis1_physiology_THA.m
% Analysis 1 (Physiology): Compare RA-THA vs C-THA physiological stress markers
%
% Metrics: HR (heart rate), BR (breathing rate), RMSSD (HRV proxy)
% Note: RVT excluded — it reflects respiratory mechanics rather than stress directly.
%
% Design:
%   - Independent groups: RA-THA vs C-THA
%   - Unit of analysis: surgeon mean per condition
%     (physiology_by_group.m computes mean across all samples in that stage group
%      for one acquisition — we average across acquisitions per surgeon x tech first)
%   - Test: Mann-Whitney U (non-parametric, small n)
%   - Effect size: Hedges' g (bias-corrected, primary outcome)
%
% Input:
%   Concatenated physiology_by_group_long.csv files from all acquisitions.
%   Script first concatenates all per-acquisition long CSVs from the output
%   directory, then averages to surgeon level before testing.
%
%   Expected columns: Acquisition, Surgeon, Joint, Tech, Group, Metric, Mean, SD, Min, Max
%
% Output:
%   analysis1_physiology_shared_HbDiff.csv   — full results
%   analysis1_physiology_unique_HbDiff.csv   — descriptive only
%   analysis1_physiology_descriptive.csv     — mean ± SD per group per metric
%
% Dylan Patel | SENSE Project

clear all; close all; clc;

%% -------------------------------------------------------------------------
%  CONFIGURATION
%% -------------------------------------------------------------------------

% Base directory containing per-acquisition physiology folders
% Each acquisition folder should contain physiology_by_group_long.csv
data_base  = '/Volumes/Dylan SSD/DYLAN/Data';
output_dir = '/Volumes/Dylan SSD/DYLAN/Results/Analysis';

if ~isfolder(output_dir), mkdir(output_dir); end

alpha = 0.05;

% Metrics to analyse
metrics           = {'HR', 'BR', 'RMSSD'};
shared_conditions = {'shared_same', 'shared_diff'};
unique_conditions = {'robotic_only', 'conventional_only'};

% Expert/novice not needed here — Analysis 1 pools all surgeons by technique
% Surgeon IDs with RA-THA or C-THA Hip procedures
% (used only for reporting — filtering done by Joint/Tech columns)

%% -------------------------------------------------------------------------
%  STEP 1: FIND AND CONCATENATE ALL PER-ACQUISITION LONG CSVs
%% -------------------------------------------------------------------------

fprintf('Searching for physiology_by_group_long.csv files...\n');

% Recursively find all matching files
csv_files = dir(fullfile(data_base, '**', 'physiology_by_group_long.csv'));

if isempty(csv_files)
    error(['No physiology_by_group_long.csv files found under:\n  %s\n' ...
           'Run physiology_by_group.m for each acquisition first.'], data_base);
end

fprintf('Found %d acquisition files.\n', length(csv_files));

master = table();
for f = 1:length(csv_files)
    fpath = fullfile(csv_files(f).folder, csv_files(f).name);
    T = readtable(fpath, 'TextType', 'string');
    master = [master; T];
    fprintf('  Loaded: %s (%d rows)\n', csv_files(f).name, height(T));
end

fprintf('Total rows: %d\n', height(master));

% Force string columns
str_cols = {'Acquisition','Surgeon','Joint','Tech','Group','Metric'};
for s = 1:length(str_cols)
    if ismember(str_cols{s}, master.Properties.VariableNames)
        master.(str_cols{s}) = string(master.(str_cols{s}));
    end
end

%% -------------------------------------------------------------------------
%  STEP 2: FILTER TO THA (Hip) AND AVERAGE TO SURGEON LEVEL
%% -------------------------------------------------------------------------
% Each surgeon may have multiple acquisitions for the same tech.
% Average their Mean values across acquisitions to get one value per
% surgeon x tech x condition x metric — this is the unit of analysis.

T_tha = master(master.Joint == "H", :);
fprintf('\nTHA rows: %d\n', height(T_tha));

% Get unique surgeons per tech
surgs_rob  = unique(T_tha.Surgeon(T_tha.Tech == "R"));
surgs_conv = unique(T_tha.Surgeon(T_tha.Tech == "C"));
fprintf('RA-THA surgeons (%d): %s\n', length(surgs_rob),  strjoin(surgs_rob,  ', '));
fprintf('C-THA  surgeons (%d): %s\n', length(surgs_conv), strjoin(surgs_conv, ', '));

% Build surgeon-level table
surg_level = table();

all_surgs = unique(T_tha.Surgeon);
all_techs = ["R", "C"];
all_conds = [shared_conditions, unique_conditions];

for s = 1:length(all_surgs)
    surg = all_surgs(s);
    for tk = 1:length(all_techs)
        tech = all_techs(tk);
        for c = 1:length(all_conds)
            cond = string(all_conds{c});
            for m = 1:length(metrics)
                met = string(metrics{m});

                mask = T_tha.Surgeon == surg & T_tha.Tech == tech & ...
                       T_tha.Group == cond & T_tha.Metric == met;

                vals = T_tha.Mean(mask);
                if isempty(vals), continue; end

                % Average across acquisitions for this surgeon x tech
                surg_mean = mean(vals, 'omitnan');
                n_acq_val = sum(~isnan(vals));

                new_row = {char(surg), char(tech), char(cond), char(met), ...
                           surg_mean, n_acq_val};
                surg_level = [surg_level; cell2table(new_row, 'VariableNames', ...
                    {'Surgeon','Tech','Condition','Metric','Mean','N_Acquisitions'})];
            end
        end
    end
end

fprintf('\nSurgeon-level rows: %d\n', height(surg_level));

% Force string columns
str_cols2 = {'Surgeon','Tech','Condition','Metric'};
for s = 1:length(str_cols2)
    surg_level.(str_cols2{s}) = string(surg_level.(str_cols2{s}));
end

%% -------------------------------------------------------------------------
%  STEP 3: DESCRIPTIVE STATISTICS
%% -------------------------------------------------------------------------

fprintf('\n--- Computing descriptive statistics ---\n');

desc_results = table();

for m = 1:length(metrics)
    met = string(metrics{m});
    for c = 1:length(all_conds)
        cond = string(all_conds{c});

        for grp = ["R", "C"]
            mask  = surg_level.Tech == grp & ...
                    surg_level.Condition == cond & ...
                    surg_level.Metric == met;
            vals  = surg_level.Mean(mask);
            n_s   = length(vals);

            if n_s == 0, continue; end

            grp_label = "RA-THA";
            if grp == "C", grp_label = "C-THA"; end

            if n_s > 1
                sd_val  = std(vals);
                sem_val = sd_val / sqrt(n_s);
            else
                sd_val = NaN; sem_val = NaN;
            end

            new_row = {char(met), char(cond), char(grp_label), ...
                mean(vals,'omitnan'), sd_val, sem_val, ...
                median(vals,'omitnan'), min(vals), max(vals), n_s};
            desc_results = [desc_results; cell2table(new_row, 'VariableNames', ...
                {'Metric','Condition','Group', ...
                 'Mean','SD','SEM','Median','Min','Max','n_surgeons'})];
        end
    end
end

%% -------------------------------------------------------------------------
%  STEP 4: MANN-WHITNEY U + HEDGES' G — SHARED CONDITIONS
%% -------------------------------------------------------------------------

fprintf('\n--- Running Mann-Whitney U + Hedges'' g (shared conditions) ---\n');

results_shared = table();

for m = 1:length(metrics)
    met = string(metrics{m});
    for c = 1:length(shared_conditions)
        cond = string(shared_conditions{c});

        mask_rob  = surg_level.Tech == "R" & surg_level.Condition == cond & surg_level.Metric == met;
        mask_conv = surg_level.Tech == "C" & surg_level.Condition == cond & surg_level.Metric == met;

        vals_rob  = surg_level.Mean(mask_rob);
        vals_conv = surg_level.Mean(mask_conv);
        n_rob     = length(vals_rob);
        n_conv    = length(vals_conv);
        diff_val  = mean(vals_rob,'omitnan') - mean(vals_conv,'omitnan');

        % ── Hedges' g (always computed) ──────────────────────────────────
        if n_rob >= 2 && n_conv >= 2
            df   = n_rob + n_conv - 2;
            sd_p = sqrt(((n_rob-1)*std(vals_rob)^2 + (n_conv-1)*std(vals_conv)^2) / df);
            g_note = sprintf('n_rob=%d, n_conv=%d', n_rob, n_conv);
        elseif n_rob >= 2 && n_conv == 1
            df = n_rob - 1; sd_p = std(vals_rob);
            g_note = 'RA-THA SD proxy (n_conv=1)';
        elseif n_rob == 1 && n_conv >= 2
            df = n_conv - 1; sd_p = std(vals_conv);
            g_note = 'C-THA SD proxy (n_rob=1)';
        else
            df = NaN; sd_p = NaN; g_note = 'n=1 both — not calculable';
        end

        if ~isnan(sd_p) && sd_p > 0
            d = diff_val / sd_p;
            J = 1 - 3 / (4*df - 1);
            g = d * J;
        else
            d = NaN; J = NaN; g = NaN;
        end

        % ── Mann-Whitney U (only when n >= 2 per group) ──────────────────
        if n_rob >= 2 && n_conv >= 2
            [p, ~, stats] = ranksum(vals_rob, vals_conv, 'alpha', alpha);
            U = stats.ranksum - n_rob*(n_rob+1)/2;
            interp_str = "Pending";
        else
            p = NaN; U = NaN;
            interp_str = "g only (n<2 in one group)";
        end

        new_row = {char(met), char(cond), ...
            mean(vals_rob,'omitnan'), mean(vals_conv,'omitnan'), diff_val, ...
            U, p, d, g, abs(g), J, interp_str, n_rob, n_conv, g_note};
        results_shared = [results_shared; cell2table(new_row, 'VariableNames', ...
            {'Metric','Condition', ...
             'Mean_Robotic','Mean_Conventional','Difference', ...
             'U_stat','p_value', ...
             'Cohens_d','Hedges_g','Abs_g','J_correction', ...
             'Interpretation','n_rob','n_conv','g_note'})];
    end
end

% Assign interpretations
for i = 1:height(results_shared)
    if strcmp(char(results_shared.Interpretation(i)), 'g only (n<2 in one group)'), continue; end
    if isnan(results_shared.p_value(i)), continue; end

    p   = results_shared.p_value(i);
    g   = results_shared.Abs_g(i);
    dif = results_shared.Difference(i);

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

    if dif > 0,      dir = "Robotic > Conventional";
    elseif dif < 0,  dir = "Conventional > Robotic";
    else,            dir = "No difference";
    end

    results_shared.Interpretation(i) = sig + " | " + eff + " | " + dir;
end

%% -------------------------------------------------------------------------
%  STEP 6: UNIQUE CONDITIONS — DESCRIPTIVE ONLY
%% -------------------------------------------------------------------------

results_unique = table();

for m = 1:length(metrics)
    met = string(metrics{m});
    for c = 1:length(unique_conditions)
        cond     = string(unique_conditions{c});
        tech_grp = "R";
        if cond == "conventional_only", tech_grp = "C"; end

        mask = surg_level.Tech == tech_grp & ...
               surg_level.Condition == cond & surg_level.Metric == met;
        vals = surg_level.Mean(mask);
        n_s  = length(vals);
        if n_s == 0, continue; end

        new_row = {char(met), char(cond), ...
            mean(vals,'omitnan'), std(vals,'omitnan'), ...
            median(vals,'omitnan'), min(vals), max(vals), n_s};
        results_unique = [results_unique; cell2table(new_row, 'VariableNames', ...
            {'Metric','Condition','Mean','SD','Median','Min','Max','n_surgeons'})];
    end
end

%% -------------------------------------------------------------------------
%  SAVE RESULTS
%% -------------------------------------------------------------------------

out_shared = fullfile(output_dir, 'analysis1_physiology_shared.csv');
out_unique = fullfile(output_dir, 'analysis1_physiology_unique.csv');
out_desc   = fullfile(output_dir, 'analysis1_physiology_descriptive.csv');

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
fprintf('  ANALYSIS 1 (PHYSIOLOGY): RA-THA vs C-THA\n');
fprintf('========================================================================\n');
fprintf('  Metrics: HR (bpm), BR (breaths/min), RMSSD (ms)\n');
fprintf('  Unit of analysis: surgeon (averaged across acquisitions)\n');
fprintf('  Effect size: Hedges'' g  |  Trivial<0.2, Small 0.2-0.5, Medium 0.5-0.8, Large>=0.8\n');
fprintf('------------------------------------------------------------------------\n');

% ── Descriptive stats ────────────────────────────────────────────────────
fprintf('\n  DESCRIPTIVE STATISTICS\n');
for m = 1:length(metrics)
    met = string(metrics{m});
    fprintf('\n  Metric: %s\n', char(met));
    fprintf('  %-22s %-12s %12s %12s %12s %12s\n', ...
        'Condition', 'Group', 'Mean', 'SD', 'Median', 'n');
    fprintf('  %s\n', repmat('-', 1, 80));

    for c = 1:length(all_conds)
        cond = string(all_conds{c});
        for grp = ["RA-THA", "C-THA"]
            mask = strcmp(desc_results.Metric, char(met)) & ...
                   strcmp(desc_results.Condition, char(cond)) & ...
                   strcmp(desc_results.Group, char(grp));
            if ~any(mask), continue; end
            fprintf('  %-22s %-12s %12.4f %12.4f %12.4f %12d\n', ...
                char(cond), char(grp), ...
                desc_results.Mean(mask), desc_results.SD(mask), ...
                desc_results.Median(mask), desc_results.n_surgeons(mask));
        end
    end
end

% ── Shared conditions: effect sizes + test ───────────────────────────────
fprintf('\n');
fprintf('========================================================================\n');
fprintf('  HEDGES'' g + MANN-WHITNEY U — SHARED CONDITIONS\n');
fprintf('========================================================================\n');

for m = 1:length(metrics)
    met = string(metrics{m});
    fprintf('\n  Metric: %s\n', char(met));
    fprintf('  %-22s %6s %6s %12s %12s %12s %7s %7s %8s  %-30s\n', ...
        'Condition', 'n_Rob', 'n_Con', 'Mean_Rob', 'Mean_Con', 'Diff', ...
        'd', 'g', 'p_value', 'Result');
    fprintf('  %s\n', repmat('-', 1, 130));

    for c = 1:length(shared_conditions)
        cond = string(shared_conditions{c});
        mask = strcmp(results_shared.Metric, char(met)) & ...
               strcmp(results_shared.Condition, char(cond));
        if ~any(mask), continue; end
        row = results_shared(mask, :);

        if isnan(row.p_value)
            p_val_str = '     N/A ';
        else
            p_val_str = sprintf('%8.4f', row.p_value);
        end

        fprintf('  %-22s %6d %6d %12.4f %12.4f %12.4f %7.3f %7.3f %s  %s\n', ...
            char(cond), row.n_rob, row.n_conv, ...
            row.Mean_Robotic, row.Mean_Conventional, row.Difference, ...
            row.Cohens_d, row.Hedges_g, p_val_str, ...
            char(row.Interpretation));
    end
end

% ── Unique conditions ────────────────────────────────────────────────────
fprintf('\n');
fprintf('========================================================================\n');
fprintf('  UNIQUE CONDITIONS (descriptive only)\n');
fprintf('========================================================================\n');
fprintf('  %-10s %-22s %12s %12s %12s %6s\n', ...
    'Metric', 'Condition', 'Mean', 'SD', 'Median', 'n');
fprintf('  %s\n', repmat('-', 1, 80));
for i = 1:height(results_unique)
    fprintf('  %-10s %-22s %12.4f %12.4f %12.4f %6d\n', ...
        results_unique.Metric{i}, results_unique.Condition{i}, ...
        results_unique.Mean(i), results_unique.SD(i), ...
        results_unique.Median(i), results_unique.n_surgeons(i));
end

% ── Significant / notable findings ───────────────────────────────────────
fprintf('\n');
fprintf('========================================================================\n');
fprintf('  NOTABLE FINDINGS (|g| >= 0.5)\n');
fprintf('========================================================================\n');
notable_mask = results_shared.Abs_g >= 0.5;
if any(notable_mask)
    T_note = results_shared(notable_mask, :);
    [~, sidx] = sort(T_note.Abs_g, 'descend');
    T_note = T_note(sidx, :);
    for i = 1:height(T_note)
        fprintf('  [%s | %s]  g=%.3f  p=%.4f  — %s\n', ...
            T_note.Metric{i}, T_note.Condition{i}, ...
            T_note.Hedges_g(i), T_note.p_value(i), ...
            char(T_note.Interpretation(i)));
    end
else
    fprintf('  No effects >= 0.5 in shared conditions.\n');
end

fprintf('\n========================================================================\n');
fprintf('  INTERPRETATION NOTES\n');
fprintf('========================================================================\n');
fprintf('  Unit of analysis: surgeon mean (averaged across acquisitions per surgeon x tech).\n');
fprintf('  SD within physiology_by_group.m reflects within-epoch sample variance —\n');
fprintf('  not used here. Between-surgeon SD used for Hedges'' g.\n');
fprintf('  RMSSD = sqrt(mean(diff(RR)^2)) — approximated from HR-derived RR intervals.\n');
fprintf('  Higher RMSSD = greater HRV = lower physiological stress.\n');
fprintf('\nDone.\n');
