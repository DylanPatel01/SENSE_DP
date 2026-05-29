%% synthesise_physiology_grouped.m
% Synthesises grouped physiology data across acquisitions:
%   Step 1: Load physiology_by_group.csv for each acquisition
%   Step 2: Average across multiple surgeries for the same
%           surgeon/joint/tech combination
%   Step 3: Output separate CSV files per surgeon x joint x tech
%
% Reads: physiology_by_group.csv from each acquisition's results folder
% Measures: HR (mean, SD, min, max), BR (mean, SD, min, max), RMSSD

clear all; close all; clc;

%% -------------------------------------------------------------------------
%  CONFIGURATION
% -------------------------------------------------------------------------

acq_path = '/Volumes/Dylan SSD/acquisitions.mat';
load(acq_path, 'acquisitions');

results_base = 'physiology grouped';
output_dir   = '/Volumes/Dylan SSD/DYLAN/Results/Synthesised_Physiology_Grouped';

if ~isfolder(output_dir), mkdir(output_dir); end

% The 4 possible grouped GLM conditions (canonical order)
all_group_conditions = {'shared_same', 'shared_diff', 'conventional_only', 'robotic_only'};

% Columns to synthesise
measures = {'HR_mean','HR_sd','HR_min','HR_max', ...
            'BR_mean','BR_sd','BR_min','BR_max', ...
            'RMSSD'};

%% -------------------------------------------------------------------------
%  PARSE ALL ACQUISITIONS
% -------------------------------------------------------------------------

n_acq = length(acquisitions);

acq_info = struct();
for i = 1:n_acq
    parts = strsplit(acquisitions(i).name, '_');
    acq_info(i).name    = acquisitions(i).name;
    acq_info(i).surgeon = parts{1};
    acq_info(i).joint   = upper(parts{2});
    acq_info(i).tech    = upper(parts{3});
    acq_info(i).key     = sprintf('%s_%s_%s', parts{1}, upper(parts{2}), upper(parts{3}));

    if strcmp(acq_info(i).joint, 'H'), jf = 'hip'; else, jf = 'knee'; end
    if strcmp(acq_info(i).tech, 'C'),  tf = 'c';   else, tf = 'r';    end

    acq_info(i).datafold = fullfile('/Volumes/Dylan SSD/DYLAN/Data', jf, tf, acq_info(i).surgeon);
    acq_info(i).phys_csv = fullfile(acq_info(i).datafold, results_base, ...
        acquisitions(i).name, 'physiology_by_group.csv');
end

all_keys    = {acq_info.key};
unique_keys = unique(all_keys);
fprintf('Found %d unique surgeon x joint x tech combinations:\n', length(unique_keys));
for k = 1:length(unique_keys)
    n = sum(strcmp(all_keys, unique_keys{k}));
    fprintf('  %s : %d surgeries\n', unique_keys{k}, n);
end

%% -------------------------------------------------------------------------
%  MASTER TABLE
% -------------------------------------------------------------------------

master_table = table();

%% -------------------------------------------------------------------------
%  PROCESS EACH SURGEON x JOINT x TECH GROUP
% -------------------------------------------------------------------------

for k = 1:length(unique_keys)
    key = unique_keys{k};
    idx = find(strcmp(all_keys, key));
    n_surgeries = length(idx);

    key_parts = strsplit(key, '_');
    surgeon   = key_parts{1};
    joint     = key_parts{2};
    tech      = key_parts{3};

    fprintf('\n===== Processing %s (%d surgeries) =====\n', key, n_surgeries);

    % ------------------------------------------------------------------
    %  Load physiology CSVs
    % ------------------------------------------------------------------
    phys_stack = {};

    for s = 1:n_surgeries
        csv_path = acq_info(idx(s)).phys_csv;
        if ~exist(csv_path, 'file')
            fprintf('  WARNING: Physiology CSV not found, skipping:\n    %s\n', csv_path);
            continue;
        end

        T = readtable(csv_path);

        % Handle different possible column names for the group/condition
        if ismember('Group', T.Properties.VariableNames)
            % already named Group
        elseif ismember('Condition', T.Properties.VariableNames)
            T.Group = T.Condition;
        elseif ismember('group', T.Properties.VariableNames)
            T.Group = T.group;
        else
            fprintf('  WARNING: No Group/Condition column found in:\n    %s\n', csv_path);
            fprintf('  Columns: %s\n', strjoin(T.Properties.VariableNames, ', '));
            continue;
        end

        % Ensure Group column is a cell array of strings
        if ~iscell(T.Group)
            T.Group = cellstr(string(T.Group));
        end
        fprintf('  Loaded: %s (%d conditions)\n', acq_info(idx(s)).name, height(T));
        phys_stack{end+1} = T;
    end

    if isempty(phys_stack)
        fprintf('  No physiology files found for %s, skipping.\n', key);
        continue;
    end

    % ------------------------------------------------------------------
    %  Find all conditions across surgeries (union)
    % ------------------------------------------------------------------
    all_conditions = {};
    for s = 1:length(phys_stack)
        all_conditions = union(all_conditions, phys_stack{s}.Group, 'stable');
    end

    % Sort in canonical order
    [~, order] = ismember(all_conditions, all_group_conditions);
    order(order == 0) = length(all_group_conditions) + 1;
    [~, sort_idx] = sort(order);
    conditions = all_conditions(sort_idx);

    fprintf('  Conditions: %s\n', strjoin(conditions, ', '));

    % ------------------------------------------------------------------
    %  Average each measure across surgeries per condition
    % ------------------------------------------------------------------
    n_conds    = length(conditions);
    n_measures = length(measures);

    % Collect values: conditions x measures x surgeries
    vals_all = NaN(n_conds, n_measures, length(phys_stack));

    for s = 1:length(phys_stack)
        T = phys_stack{s};
        t_groups = T.Group;

        for g = 1:n_conds
            row_idx = find(strcmp(t_groups, conditions{g}));
            if isempty(row_idx), continue; end

            for m = 1:n_measures
                if ismember(measures{m}, T.Properties.VariableNames)
                    val = T.(measures{m})(row_idx);
                    if ~isnan(val) && val ~= 0
                        vals_all(g, m, s) = val;
                    end
                end
            end
        end
    end

    % Average across surgeries
    avg_vals   = mean(vals_all, 3, 'omitnan');
    sd_vals    = std(vals_all, 0, 3, 'omitnan');
    n_valid    = sum(~isnan(vals_all), 3);

    % ------------------------------------------------------------------
    %  Build output table for this surgeon x joint x tech
    % ------------------------------------------------------------------
    out_table = table();

    for g = 1:n_conds
        row = table();
        row.Surgeon    = {surgeon};
        row.Joint      = {joint};
        row.Tech       = {tech};
        row.Group_Key  = {key};
        row.Condition  = conditions(g);
        row.N_Surgeries = length(phys_stack);

        for m = 1:n_measures
            row.(measures{m}) = avg_vals(g, m);
        end

        % Add cross-surgery SD for the key measures (HR_mean, BR_mean, RMSSD)
        % to indicate variability across surgeries
        row.HR_mean_across_surgery_SD  = sd_vals(g, strcmp(measures, 'HR_mean'));
        row.BR_mean_across_surgery_SD  = sd_vals(g, strcmp(measures, 'BR_mean'));
        row.RMSSD_across_surgery_SD    = sd_vals(g, strcmp(measures, 'RMSSD'));

        out_table = [out_table; row];
    end

    % ------------------------------------------------------------------
    %  SAVE: Per surgeon x joint x tech
    % ------------------------------------------------------------------
    out_file = fullfile(output_dir, sprintf('phys_grouped_%s.csv', key));
    writetable(out_table, out_file);
    fprintf('  Saved: %s\n', out_file);

    % Append to master
    master_table = [master_table; out_table];

    % ------------------------------------------------------------------
    %  DISPLAY PREVIEW
    % ------------------------------------------------------------------
    fprintf('\n  --- Physiology preview for %s ---\n', key);
    disp(out_table(:, {'Condition','HR_mean','HR_sd','BR_mean','BR_sd','RMSSD','N_Surgeries'}));
end

%% -------------------------------------------------------------------------
%  SAVE MASTER TABLE
% -------------------------------------------------------------------------

out_master = fullfile(output_dir, 'master_physiology_grouped.csv');
writetable(master_table, out_master);
fprintf('\nSaved master physiology table: %s\n', out_master);

%% -------------------------------------------------------------------------
%  SUMMARY
% -------------------------------------------------------------------------

fprintf('\n========================================\n');
fprintf('  GROUPED PHYSIOLOGY SYNTHESIS COMPLETE\n');
fprintf('========================================\n');
fprintf('Acquisitions : %d total\n', n_acq);
fprintf('Groups       : %d unique surgeon x joint x tech\n', length(unique_keys));
fprintf('Conditions   : %s\n', strjoin(all_group_conditions, ', '));
fprintf('Measures     : %s\n', strjoin(measures, ', '));
fprintf('\nPer-group files saved to:\n  %s\n', output_dir);
fprintf('\nFiles per group:\n');
fprintf('  phys_grouped_<key>.csv — averaged physiology per condition\n');
fprintf('\nMaster file:\n');
fprintf('  master_physiology_grouped.csv — all data combined (long format)\n');
fprintf('\nDone.\n');