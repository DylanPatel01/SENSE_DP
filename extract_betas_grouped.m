%% extract_betas_grouped.m
% Extracts GLM beta values (bf1 = canonical HRF) from the GROUPED
% spm_fnirs results and saves them as CSV files.
%
% USAGE: Run this script after standard_glm_grouped.m completes.
%
% OUTPUT:
%   betas_grouped_HbDiff_bf1.csv           - all channels, rows = groups
%   betas_grouped_HbDiff_bf1_validch.csv   - valid channels only
%   beta_grouped_summary.txt               - summary
%
% Dylan Patel | SENSE Project

clear all; close all; clc;

%% -------------------------------------------------------------------------
%  INPUT PARAMETERS — reads from acquisitions struct
% -------------------------------------------------------------------------

acq_path = '/Volumes/Dylan SSD/acquisitions.mat';
load(acq_path, 'acquisitions');

acq_idx = 28;  % <-- CHANGE THIS to match the acquisition you ran the grouped GLM on

% --- Alternative: select by name ---
% acq_name = '1_H_C_101025_1';
% acq_idx  = find(strcmp({acquisitions.name}, acq_name));
% if isempty(acq_idx), error('Acquisition "%s" not found.', acq_name); end

A = acquisitions(acq_idx);
fprintf('Extracting grouped betas for: %s\n', A.name);

% Parse joint and tech
name_parts = strsplit(A.name, '_');
joint_str  = upper(name_parts{2});
tech_str   = upper(name_parts{3});
ptp        = name_parts{1};

if strcmpi(joint_str, 'H'),     joint_folder = 'hip';
elseif strcmpi(joint_str, 'K'), joint_folder = 'knee';
else, error('Unknown joint code: %s', joint_str); end

if strcmpi(tech_str, 'C'),     tech_folder = 'c';
elseif strcmpi(tech_str, 'R'), tech_folder = 'r';
else, error('Unknown tech code: %s', tech_str); end

signal         = 'HbDiff';
results_folder = ['GLM TDD Phys Grouped' filesep A.name];
n_nuisance     = 4;  % hr, br, rvt, constant (set to 1 if add_regs=0)

datafold   = fullfile('/Volumes/Dylan SSD/DYLAN/Data', joint_folder, tech_folder, ptp);
spm_mat    = fullfile(datafold, results_folder, results_folder, signal, 'SPM.mat');
nirs_mat   = fullfile(datafold, results_folder, 'NIRS.mat');
pos_mat    = fullfile(datafold, 'Digitisation', 'POS.mat');
output_dir = fullfile(datafold, results_folder);

%% -------------------------------------------------------------------------
%  LOAD DATA
% -------------------------------------------------------------------------

fprintf('Loading SPM.mat from:\n  %s\n', spm_mat);
if ~exist(spm_mat, 'file'), error('SPM.mat not found:\n  %s', spm_mat); end
load(spm_mat);

fprintf('Loading P from NIRS.mat...\n');
if ~exist(nirs_mat, 'file'), error('NIRS.mat not found:\n  %s', nirs_mat); end
load(nirs_mat, 'P');

fprintf('Loading R from POS.mat...\n');
if ~exist(pos_mat, 'file'), error('POS.mat not found:\n  %s', pos_mat); end
load(pos_mat, 'R');

%% -------------------------------------------------------------------------
%  EXTRACT BETA VALUES
% -------------------------------------------------------------------------

fprintf('\nFound %d regressors, %d channels.\n', size(SPM.beta,1), size(SPM.beta,2));

% bf(1) indices = canonical HRF only (every 3rd regressor, excluding nuisance)
n_event_regs = size(SPM.beta,1) - n_nuisance;
bf1_idx      = 1:3:n_event_regs;

event_names = SPM.xX.name(bf1_idx)';
fprintf('Found %d group conditions (bf1 regressors).\n', length(bf1_idx));

% Beta matrix: groups x channels
betas_bf1 = SPM.beta(bf1_idx, :);

% Clean event names (remove 'Sn(1) ' prefix and '*bf(1)' suffix)
clean_events = regexprep(event_names, 'Sn\(1\) ', '');
clean_events = regexprep(clean_events, '\*bf\(1\)', '');

% Channel labels
ch_labels_all = arrayfun(@(x) sprintf('ch%02d', x), R.ch.label, 'UniformOutput', false);

%% -------------------------------------------------------------------------
%  SAVE: ALL CHANNELS
% -------------------------------------------------------------------------

beta_table_all = array2table(betas_bf1, ...
    'RowNames',      clean_events, ...
    'VariableNames', ch_labels_all);

out_all = fullfile(output_dir, sprintf('betas_grouped_%s_bf1.csv', signal));
writetable(beta_table_all, out_all, 'WriteRowNames', true);
fprintf('\nSaved (all channels):\n  %s\n', out_all);

%% -------------------------------------------------------------------------
%  SAVE: VALID CHANNELS ONLY
% -------------------------------------------------------------------------

valid_ch     = find(R.ch.mask == 1);
valid_labels = R.ch.label(valid_ch);
betas_valid  = betas_bf1(:, valid_ch);

ch_labels_valid = arrayfun(@(x) sprintf('ch%02d', x), valid_labels, 'UniformOutput', false);

beta_table_valid = array2table(betas_valid, ...
    'RowNames',      clean_events, ...
    'VariableNames', ch_labels_valid);

out_valid = fullfile(output_dir, sprintf('betas_grouped_%s_bf1_validch.csv', signal));
writetable(beta_table_valid, out_valid, 'WriteRowNames', true);
fprintf('Saved (valid channels only):\n  %s\n', out_valid);

%% -------------------------------------------------------------------------
%  SAVE: LONG-FORMAT TABLE (for easy stats in R / Python / Excel)
% -------------------------------------------------------------------------
% Creates a tall table with columns:
%   Acquisition | Group | Channel | Beta
% This makes it straightforward to pivot, filter, or run stats.

n_groups_found = length(clean_events);
n_valid_ch     = length(valid_ch);

acq_col   = repmat({A.name}, n_groups_found * n_valid_ch, 1);
surgeon_col = repmat(name_parts(1), n_groups_found * n_valid_ch, 1);
joint_col = repmat({joint_str}, n_groups_found * n_valid_ch, 1);
tech_col  = repmat({tech_str}, n_groups_found * n_valid_ch, 1);
group_col = {};
ch_col    = {};
beta_col  = [];

row = 0;
for g = 1:n_groups_found
    for c = 1:n_valid_ch
        row = row + 1;
        group_col{row,1} = clean_events{g};
        ch_col{row,1}    = ch_labels_valid{c};
        beta_col(row,1)  = betas_valid(g, c);
    end
end

long_table = table(acq_col, surgeon_col, joint_col, tech_col, group_col, ch_col, beta_col, ...
    'VariableNames', {'Acquisition', 'Surgeon', 'Joint', 'Tech', 'Group', 'Channel', 'Beta'});

out_long = fullfile(output_dir, sprintf('betas_grouped_%s_bf1_long.csv', signal));
writetable(long_table, out_long);
fprintf('Saved (long format):\n  %s\n', out_long);

%% -------------------------------------------------------------------------
%  SAVE: SUMMARY TEXT FILE
% -------------------------------------------------------------------------

out_summary = fullfile(output_dir, 'beta_grouped_summary.txt');
fid = fopen(out_summary, 'w');
fprintf(fid, 'Grouped Beta Extraction Summary\n');
fprintf(fid, '===============================\n');
fprintf(fid, 'Acquisition : %s\n', A.name);
fprintf(fid, 'Surgeon     : %s\n', ptp);
fprintf(fid, 'Joint       : %s\n', joint_str);
fprintf(fid, 'Tech        : %s\n', tech_str);
fprintf(fid, 'Signal      : %s\n', signal);
fprintf(fid, 'Date        : %s\n\n', datestr(now));
fprintf(fid, 'Total channels    : %d\n', size(SPM.beta,2));
fprintf(fid, 'Valid channels    : %d\n', n_valid_ch);
fprintf(fid, 'Valid channel IDs : %s\n\n', num2str(valid_labels));
fprintf(fid, 'Groups extracted (%d):\n', n_groups_found);
for i = 1:n_groups_found
    fprintf(fid, '  %d. %s\n', i, clean_events{i});
end
fprintf(fid, '\nNuisance regressors (not extracted): hr, br, rvt, constant\n');
fprintf(fid, '\nOutput files:\n');
fprintf(fid, '  Wide (all ch):   %s\n', out_all);
fprintf(fid, '  Wide (valid ch): %s\n', out_valid);
fprintf(fid, '  Long format:     %s\n', out_long);
fclose(fid);
fprintf('Saved summary:\n  %s\n', out_summary);

%% -------------------------------------------------------------------------
%  DISPLAY PREVIEW
% -------------------------------------------------------------------------

fprintf('\n--- Grouped beta preview (valid channels) ---\n');
disp(beta_table_valid(:, 1:min(8,end)));

fprintf('\n--- Long format preview (first 10 rows) ---\n');
disp(long_table(1:min(10,end), :));

fprintf('\nDone. Files saved to:\n  %s\n', output_dir);