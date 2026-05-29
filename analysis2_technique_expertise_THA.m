%% compute_cohens_d_grouped.m
% Companion to synthesise_betas_grouped_SD.m
%
% Computes Hedges' g (bias-corrected Cohen's d) for robotic vs conventional
% comparison within each surgeon x joint group, per brain region and stage.
%
% Hedges' g formula:
%   g = d * J   where J = 1 - 3/(4*df - 1),   df = n_ch_rob + n_ch_conv - 2
%   d = (Mean_Robotic - Mean_Conventional) / SD_pooled
%   SD_pooled = sqrt(((n1-1)*SD1^2 + (n2-1)*SD2^2) / df)
%
% Note on n and SD:
%   SD here is SPATIAL variability across channels within each region.
%   n used for Hedges' correction = N_Channels per group per region.
%   This is not independent replication — effect sizes are descriptive only.
%   For single-surgeon groups (N_Surgeries=1), Hedges' g is still reported
%   but flagged — the correction uses channel n as a proxy.
%
% Benchmarks: Trivial<0.2, Small 0.2-0.5, Medium 0.5-0.8, Large>=0.8
%
% Dylan Patel | SENSE Project

clear all; close all; clc;

%% -------------------------------------------------------------------------
%  CONFIGURATION
%% -------------------------------------------------------------------------

signal      = 'HbDiff';
output_dir  = '/Volumes/Dylan SSD/DYLAN/Results/Synthesised_Grouped';
master_file = fullfile(output_dir, sprintf('master_grouped_region_%s.csv', signal));
out_file    = fullfile(output_dir, sprintf('hedges_g_grouped_%s.csv', signal));

comparisons = {
    'shared_same', 'shared_same';
    'shared_diff', 'shared_diff';
};

%% -------------------------------------------------------------------------
%  LOAD MASTER REGION TABLE
%% -------------------------------------------------------------------------

if ~exist(master_file, 'file')
    error('Master region file not found:\n  %s\nRun synthesise_betas_grouped_SD.m first.', master_file);
end

T = readtable(master_file, 'TextType', 'string');
fprintf('Loaded master region table: %d rows\n', height(T));

str_cols = {'Surgeon','Joint','Tech','Group_Key','Condition','Region','Region_Label'};
for sc = 1:length(str_cols)
    if ismember(str_cols{sc}, T.Properties.VariableNames)
        T.(str_cols{sc}) = string(T.(str_cols{sc}));
    end
end

sj_keys   = strcat(T.Surgeon, '_', T.Joint);
unique_sj = unique(sj_keys);
regions   = unique(T.Region);
n_regions = length(regions);

fprintf('\nFound %d surgeon x joint groups\n', length(unique_sj));
fprintf('Regions: %d\n', n_regions);

%% -------------------------------------------------------------------------
%  COMPUTE HEDGES' G FOR EACH SURGEON x JOINT x STAGE x REGION
%% -------------------------------------------------------------------------

results = table();

for sj = 1:length(unique_sj)
    sj_key  = unique_sj{sj};
    sj_mask = strcmp(sj_keys, sj_key);
    T_sj    = T(sj_mask, :);

    sj_parts = strsplit(char(sj_key), '_');
    surgeon  = sj_parts{1};
    joint    = sj_parts{2};

    fprintf('\n===== %s =====\n', sj_key);

    for comp = 1:size(comparisons, 1)
        stage_label = comparisons{comp, 1};

        rob_mask  = strcmp(T_sj.Tech, 'R') & strcmp(T_sj.Condition, stage_label);
        conv_mask = strcmp(T_sj.Tech, 'C') & strcmp(T_sj.Condition, stage_label);
        T_rob     = T_sj(rob_mask,  :);
        T_conv    = T_sj(conv_mask, :);

        if isempty(T_rob)
            fprintf('  [%s] No robotic data — skipping\n', stage_label);
            continue;
        end
        if isempty(T_conv)
            fprintf('  [%s] No conventional data — skipping\n', stage_label);
            continue;
        end

        fprintf('  [%s] Robotic: %d region rows | Conventional: %d region rows\n', ...
            stage_label, height(T_rob), height(T_conv));

        for r = 1:n_regions
            reg      = regions{r};
            rob_row  = T_rob(strcmp(T_rob.Region,  reg), :);
            conv_row = T_conv(strcmp(T_conv.Region, reg), :);

            if isempty(rob_row) || isempty(conv_row)
                fprintf('    Region %s: missing data for one condition\n', reg);
                continue;
            end

            mean_rob  = mean(rob_row.Beta);
            mean_conv = mean(conv_row.Beta);
            sd_rob    = mean(rob_row.Beta_SD);
            sd_conv   = mean(conv_row.Beta_SD);
            n_ch_rob  = mean(rob_row.N_Channels);
            n_ch_conv = mean(conv_row.N_Channels);

            difference = mean_rob - mean_conv;

            % Pooled SD weighted by channel counts
            df = n_ch_rob + n_ch_conv - 2;
            if df > 0
                sd_pooled = sqrt(((n_ch_rob-1)*sd_rob^2 + (n_ch_conv-1)*sd_conv^2) / df);
            else
                sd_pooled = sqrt((sd_rob^2 + sd_conv^2) / 2);
                df = 1;  % fallback for correction factor
            end

            % Cohen's d
            if sd_pooled > 0
                cohens_d = difference / sd_pooled;
            else
                cohens_d = NaN;
                fprintf('    WARNING: SD pooled = 0 for region %s [%s]\n', reg, stage_label);
            end

            % Hedges' g bias correction
            % J = 1 - 3/(4*df - 1)
            J        = 1 - 3 / (4*df - 1);
            hedges_g = cohens_d * J;
            abs_g    = abs(hedges_g);

            % Flag if single-surgeon (correction uses channel n as proxy)
            n_surg_rob  = mean(rob_row.N_Surgeries);
            n_surg_conv = mean(conv_row.N_Surgeries);
            if n_surg_rob == 1 || n_surg_conv == 1
                g_note = 'Channel-n proxy (single surgeon)';
            else
                g_note = 'Multi-surgeon';
            end

            % Interpretation using standard Cohen benchmarks
            if isnan(abs_g)
                interpretation = 'N/A';
            elseif abs_g < 0.2
                interpretation = 'Trivial';
            elseif abs_g < 0.5
                interpretation = 'Small';
            elseif abs_g < 0.8
                interpretation = 'Medium';
            else
                interpretation = 'Large';
            end

            if isnan(hedges_g)
                direction = 'N/A';
            elseif hedges_g > 0
                direction = 'Robotic > Conventional';
            elseif hedges_g < 0
                direction = 'Conventional > Robotic';
            else
                direction = 'No difference';
            end

            reg_label = char(rob_row.Region_Label(1));

            new_row = {surgeon, joint, stage_label, reg, reg_label, ...
                mean_rob, mean_conv, sd_rob, sd_conv, sd_pooled, ...
                difference, cohens_d, J, hedges_g, abs_g, ...
                interpretation, direction, g_note, ...
                n_ch_rob, n_ch_conv, df, n_surg_rob, n_surg_conv};

            results = [results; cell2table(new_row, 'VariableNames', { ...
                'Surgeon','Joint','Stage','Region','Region_Label', ...
                'Mean_Robotic','Mean_Conventional','SD_Robotic','SD_Conventional','SD_Pooled', ...
                'Difference','Cohens_d','J_correction','Hedges_g','Abs_g', ...
                'Interpretation','Direction','Note', ...
                'N_Channels_Rob','N_Channels_Conv','df', ...
                'N_Surgeries_Rob','N_Surgeries_Conv'})];
        end
    end
end

%% -------------------------------------------------------------------------
%  SAVE RESULTS
%% -------------------------------------------------------------------------

writetable(results, out_file);
fprintf('\nSaved Hedges'' g results: %s\n', out_file);
fprintf('Total rows: %d\n', height(results));

%% -------------------------------------------------------------------------
%  DISPLAY SUMMARY
%% -------------------------------------------------------------------------

fprintf('\n========================================\n');
fprintf('  HEDGES'' G SUMMARY\n');
fprintf('========================================\n');
fprintf('  Benchmarks: Trivial<0.2, Small 0.2-0.5, Medium 0.5-0.8, Large>=0.8\n');
fprintf('  d = Cohen''s d (uncorrected)\n');
fprintf('  J = bias correction factor = 1 - 3/(4*df-1)\n');
fprintf('  g = d * J  (Hedges'' g, bias-corrected)\n');

sj_list = unique(strcat(results.Surgeon, '_', results.Joint));

for sj = 1:length(sj_list)
    sj_key   = sj_list{sj};
    sj_parts = strsplit(char(sj_key), '_');
    fprintf('\n--- %s ---\n', sj_key);

    for comp = 1:size(comparisons, 1)
        stage = comparisons{comp, 1};
        mask  = strcmp(results.Surgeon, sj_parts{1}) & ...
                strcmp(results.Joint,   sj_parts{2}) & ...
                strcmp(results.Stage,   stage);
        T_sub = results(mask, :);

        if isempty(T_sub), continue; end

        fprintf('\n  Stage: %s\n', stage);
        fprintf('  %-40s %12s %8s %6s %8s %8s  %-22s  %s\n', ...
            'Region', 'Difference', 'd', 'J', 'g', '|g|', 'Interpretation', 'Note');
        fprintf('  %s\n', repmat('-', 1, 120));

        for r = 1:height(T_sub)
            fprintf('  %-40s %12.4e %8.3f %6.3f %8.3f %8.3f  %-22s  %s (%s)\n', ...
                char(T_sub.Region_Label(r)), T_sub.Difference(r), ...
                T_sub.Cohens_d(r), T_sub.J_correction(r), ...
                T_sub.Hedges_g(r), T_sub.Abs_g(r), ...
                char(T_sub.Interpretation(r)), ...
                char(T_sub.Direction(r)), char(T_sub.Note(r)));
        end
    end
end

fprintf('\n========================================\n');
fprintf('  NOTE ON INTERPRETATION\n');
fprintf('========================================\n');
fprintf('  SD = spatial variability across channels within each region.\n');
fprintf('  n for Hedges'' correction = N_Channels per group per region.\n');
fprintf('  Channels are not independent — effect sizes are descriptive only.\n');
fprintf('  Single-surgeon groups flagged as "Channel-n proxy".\n');
fprintf('\nDone.\n');