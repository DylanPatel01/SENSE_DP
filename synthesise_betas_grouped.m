%% synthesise_betas_grouped.m
% Synthesises fNIRS GROUPED GLM beta values across acquisitions:
%   Step 1: Average betas across multiple surgeries for the same
%           surgeon/joint/tech combination
%   Step 2: Group channels into brain regions (Brodmann areas)
%   Step 3: Output separate CSV files per surgeon x joint x tech
%
% The grouped GLM has 4 possible conditions as rows:
%   shared_same, shared_diff, conventional_only, robotic_only
% (not all conditions appear in every acquisition)
%
% Reads beta CSVs produced by extract_betas_grouped.m

clear all; close all; clc;

%% -------------------------------------------------------------------------
%  CONFIGURATION
% -------------------------------------------------------------------------

acq_path = '/Volumes/Dylan SSD/acquisitions.mat';
load(acq_path, 'acquisitions');

signal       = 'HbDiff';
results_base = 'GLM TDD Phys Grouped';
beta_prefix  = 'betas_grouped_';
output_dir   = '/Volumes/Dylan SSD/DYLAN/Results/Synthesised_Grouped';

if ~isfolder(output_dir), mkdir(output_dir); end

%% -------------------------------------------------------------------------
%  BRAIN REGION DEFINITIONS (channel -> region mapping)
% -------------------------------------------------------------------------

regions = struct();
regions(1).name     = 'DLPFC_BA9_46';
regions(1).label    = 'Dorsolateral PFC (BA9, BA46)';
regions(1).channels = [1, 2, 3, 4, 5, 7, 20, 22, 23, 24, 26];

regions(2).name     = 'OFC_BA11_47';
regions(2).label    = 'Orbitofrontal Cortex (BA11, BA47)';
regions(2).channels = [2, 11, 17, 25];

regions(3).name     = 'MPFC_BA9_10';
regions(3).label    = 'Medial Prefrontal Cortex (BA9, BA10)';
regions(3).channels = [8, 9, 13, 14, 15, 18, 19];

regions(4).name     = 'SFG_BA8_9';
regions(4).label    = 'Superior Frontal Gyrus (BA8, BA9)';
regions(4).channels = [6, 10, 12, 16, 21];

regions(5).name     = 'PSC_BA1_2_3';
regions(5).label    = 'Primary Somatosensory Cortex (BA1, 2, 3)';
regions(5).channels = [33, 34, 35, 40, 41, 42];

regions(6).name     = 'SPL_BA5_7';
regions(6).label    = 'Superior Parietal Lobule (BA5, 7)';
regions(6).channels = [27, 31, 34, 37, 39, 41, 43, 45];

regions(7).name     = 'Angular_Gyrus_BA39';
regions(7).label    = 'Angular Gyrus (BA39)';
regions(7).channels = [28, 29, 30, 32, 36, 38, 44, 46, 47, 48];

regions(8).name     = 'SMG_BA40';
regions(8).label    = 'Supramarginal Gyrus (BA40)';
regions(8).channels = [29, 30, 33, 40, 47, 48];

n_regions = length(regions);

% The 4 possible grouped GLM conditions
all_group_conditions = {'shared_same', 'shared_diff', 'conventional_only', 'robotic_only'};

%% -------------------------------------------------------------------------
%  PARSE ALL ACQUISITIONS INTO SURGEON x JOINT x TECH GROUPS
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
    acq_info(i).results  = fullfile(acq_info(i).datafold, results_base, acquisitions(i).name);
    acq_info(i).beta_csv = fullfile(acq_info(i).results, ...
        sprintf('%s%s_bf1.csv', beta_prefix, signal));
end

all_keys    = {acq_info.key};
unique_keys = unique(all_keys);
fprintf('Found %d unique surgeon x joint x tech combinations:\n', length(unique_keys));
for k = 1:length(unique_keys)
    n = sum(strcmp(all_keys, unique_keys{k}));
    fprintf('  %s : %d surgeries\n', unique_keys{k}, n);
end

%% -------------------------------------------------------------------------
%  MASTER TABLES
% -------------------------------------------------------------------------

master_channel = table();
master_region  = table();

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
    %  STEP 1: Load betas from each surgery
    % ------------------------------------------------------------------
    beta_stack = {};
    condition_names_all = {};

    for s = 1:n_surgeries
        csv_path = acq_info(idx(s)).beta_csv;
        if ~exist(csv_path, 'file')
            fprintf('  WARNING: Beta CSV not found, skipping:\n    %s\n', csv_path);
            continue;
        end

        T = readtable(csv_path, 'ReadRowNames', true);
        fprintf('  Loaded: %s (%d conditions x %d channels)\n', ...
            acq_info(idx(s)).name, size(T,1), size(T,2));

        beta_stack{end+1} = T;
        condition_names_all{end+1} = T.Properties.RowNames;
    end

    if isempty(beta_stack)
        fprintf('  No beta files found for %s, skipping.\n', key);
        continue;
    end

    % ------------------------------------------------------------------
    %  Find common conditions and channels across surgeries
    % ------------------------------------------------------------------
    % For grouped GLM, conditions are shared_same, shared_diff, etc.
    % Use union of all conditions (not intersection) since different
    % surgeries might have different active groups
    all_conditions = {};
    for s = 1:length(condition_names_all)
        all_conditions = union(all_conditions, condition_names_all{s}, 'stable');
    end

    % Sort conditions in canonical order
    [~, order] = ismember(all_conditions, all_group_conditions);
    order(order == 0) = length(all_group_conditions) + 1;  % unknown conditions go last
    [~, sort_idx] = sort(order);
    common_conditions = all_conditions(sort_idx);

    fprintf('  Conditions found: %s\n', strjoin(common_conditions, ', '));

    % Common channels
    common_channels = beta_stack{1}.Properties.VariableNames;
    for s = 2:length(beta_stack)
        common_channels = intersect(common_channels, ...
            beta_stack{s}.Properties.VariableNames, 'stable');
    end
    fprintf('  Common channels: %d\n', length(common_channels));

    % ------------------------------------------------------------------
    %  Average betas across surgeries (NaN for missing)
    % ------------------------------------------------------------------
    n_conds = length(common_conditions);
    n_ch    = length(common_channels);
    beta_all = NaN(n_conds, n_ch, length(beta_stack));

    for s = 1:length(beta_stack)
        T = beta_stack{s};
        for g = 1:n_conds
            row_idx = find(strcmp(T.Properties.RowNames, common_conditions{g}));
            if isempty(row_idx), continue; end  % this surgery didn't have this condition
            for c = 1:n_ch
                col_idx = find(strcmp(T.Properties.VariableNames, common_channels{c}));
                if ~isempty(col_idx)
                    val = T{row_idx, col_idx};
                    if val ~= 0  % treat zero as missing/excluded
                        beta_all(g, c, s) = val;
                    end
                end
            end
        end
    end
    beta_avg = mean(beta_all, 3, 'omitnan');

    % Create channel-level table
    avg_table = array2table(beta_avg, ...
        'RowNames', common_conditions, ...
        'VariableNames', common_channels);

    % ------------------------------------------------------------------
    %  SAVE: Channel-level averaged betas
    % ------------------------------------------------------------------
    out_ch = fullfile(output_dir, sprintf('grouped_betas_avg_%s_%s.csv', signal, key));
    writetable(avg_table, out_ch, 'WriteRowNames', true);
    fprintf('  Saved channel-level: %s\n', out_ch);

    % ------------------------------------------------------------------
    %  STEP 2: Group channels into brain regions
    % ------------------------------------------------------------------
    ch_nums = cellfun(@(x) str2double(x(3:end)), common_channels);

    region_avg = NaN(n_conds, n_regions);
    region_labels = cell(1, n_regions);
    region_ch_count = zeros(1, n_regions);

    for r = 1:n_regions
        region_labels{r} = regions(r).name;

        [~, col_idx] = ismember(regions(r).channels, ch_nums);
        col_idx = col_idx(col_idx > 0);

        if ~isempty(col_idx)
            region_avg(:, r) = mean(beta_avg(:, col_idx), 2, 'omitnan');
            region_ch_count(r) = sum(~isnan(beta_avg(1, col_idx)));
        else
            fprintf('  WARNING: No channels found for region %s\n', regions(r).label);
        end
    end

    region_table = array2table(region_avg, ...
        'RowNames', common_conditions, ...
        'VariableNames', region_labels);

    % ------------------------------------------------------------------
    %  SAVE: Region-level averaged betas
    % ------------------------------------------------------------------
    out_reg = fullfile(output_dir, sprintf('grouped_betas_region_%s_%s.csv', signal, key));
    writetable(region_table, out_reg, 'WriteRowNames', true);
    fprintf('  Saved region-level: %s\n', out_reg);

    % ------------------------------------------------------------------
    %  APPEND TO MASTER TABLES (long format)
    % ------------------------------------------------------------------
    % Channel-level
    for g = 1:n_conds
        for c = 1:n_ch
            if isnan(beta_avg(g, c)), continue; end  % skip missing
            row_data = {surgeon, joint, tech, key, common_conditions{g}, ...
                common_channels{c}, ch_nums(c), beta_avg(g, c), n_surgeries};
            master_channel = [master_channel; cell2table(row_data, ...
                'VariableNames', {'Surgeon','Joint','Tech','Group_Key', ...
                'Condition','Channel','Channel_Num','Beta','N_Surgeries'})];
        end
    end

    % Region-level
    for g = 1:n_conds
        for r = 1:n_regions
            if isnan(region_avg(g, r)), continue; end  % skip missing
            row_data = {surgeon, joint, tech, key, common_conditions{g}, ...
                regions(r).name, regions(r).label, region_avg(g, r), ...
                region_ch_count(r), n_surgeries};
            master_region = [master_region; cell2table(row_data, ...
                'VariableNames', {'Surgeon','Joint','Tech','Group_Key', ...
                'Condition','Region','Region_Label','Beta','N_Channels','N_Surgeries'})];
        end
    end

    % ------------------------------------------------------------------
    %  DISPLAY PREVIEW
    % ------------------------------------------------------------------
    fprintf('\n  --- Region beta preview for %s ---\n', key);
    disp(region_table);
end

%% -------------------------------------------------------------------------
%  SAVE MASTER TABLES
% -------------------------------------------------------------------------

out_master_ch = fullfile(output_dir, sprintf('master_grouped_channel_%s.csv', signal));
writetable(master_channel, out_master_ch);
fprintf('\nSaved master channel table: %s\n', out_master_ch);

out_master_reg = fullfile(output_dir, sprintf('master_grouped_region_%s.csv', signal));
writetable(master_region, out_master_reg);
fprintf('Saved master region table: %s\n', out_master_reg);

%% -------------------------------------------------------------------------
%  SUMMARY
% -------------------------------------------------------------------------

fprintf('\n========================================\n');
fprintf('  GROUPED SYNTHESIS COMPLETE\n');
fprintf('========================================\n');
fprintf('Signal       : %s\n', signal);
fprintf('Acquisitions : %d total\n', n_acq);
fprintf('Groups       : %d unique surgeon x joint x tech\n', length(unique_keys));
fprintf('Conditions   : %s\n', strjoin(all_group_conditions, ', '));
fprintf('Regions      : %d brain regions\n', n_regions);
fprintf('\nPer-group files saved to:\n  %s\n', output_dir);
fprintf('\nFiles per group:\n');
fprintf('  grouped_betas_avg_<signal>_<key>.csv    — channel-level\n');
fprintf('  grouped_betas_region_<signal>_<key>.csv — region-level\n');
fprintf('\nMaster files:\n');
fprintf('  master_grouped_channel_<signal>.csv — all channel data (long)\n');
fprintf('  master_grouped_region_<signal>.csv  — all region data (long)\n');
fprintf('\nRegion definitions:\n');
for r = 1:n_regions
    ch_str = sprintf('%d,', regions(r).channels);
    fprintf('  %d. %-45s channels: [%s]\n', r, regions(r).label, ch_str(1:end-1));
end
fprintf('\nDone.\n');