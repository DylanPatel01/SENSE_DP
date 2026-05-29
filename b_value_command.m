%% extract_betas.m
% Extracts GLM beta values (bf1 = canonical HRF) from spm_fnirs results
% and saves them as CSV files after running standard_glm_analysis_phys_Dylan.m
%
% USAGE: Run this script after standard_glm_analysis_phys completes.
% Edit the INPUT PARAMETERS section below before running.
%
% OUTPUT:
%   betas_HbDiff_bf1.csv          - all channels
%   betas_HbDiff_bf1_validch.csv  - valid (spatially registered) channels only
%   beta_summary.txt              - summary of channels and events
%
% Dylan Patel | SENSE Project

clear all; close all; clc;

%% -------------------------------------------------------------------------
%  INPUT PARAMETERS — now reads from acquisitions struct
% -------------------------------------------------------------------------

acq_path = '/Volumes/Dylan SSD/acquisitions.mat';
load(acq_path, 'acquisitions');

acq_idx = 28;  % <-- CHANGE THIS to pick participant

% --- Alternative: select by name ---
% acq_name = '1_H_C_101025_1';
% acq_idx  = find(strcmp({acquisitions.name}, acq_name));
% if isempty(acq_idx), error('Acquisition "%s" not found.', acq_name); end

A = acquisitions(acq_idx);
fprintf('Extracting betas for: %s\n', A.name);

% Parse joint and tech from name (format: surgeon_joint_tech_date_timebin)
name_parts = strsplit(A.name, '_');
joint_str  = name_parts{2};  % 'H' or 'K'
tech_str   = name_parts{3};  % 'C' or 'R'

if strcmpi(joint_str, 'H'),     joint_folder = 'hip';
elseif strcmpi(joint_str, 'K'), joint_folder = 'knee';
else, error('Unknown joint code: %s', joint_str); end

if strcmpi(tech_str, 'C'),     tech_folder = 'c';
elseif strcmpi(tech_str, 'R'), tech_folder = 'r';
else, error('Unknown tech code: %s', tech_str); end

ptp            = name_parts{1};  % surgeon ID
signal         = 'HbDiff';       % HbDiff, HbO, or HbR
results_folder = ['GLM TDD Phys' filesep A.name];
n_nuisance     = 4;              % nuisance regressors at end: hr, br, rvt, constant
                                 % set to 1 if add_regs=0 (only constant)

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

event_names  = SPM.xX.name(bf1_idx)';
fprintf('Found %d events (bf1 regressors).\n', length(bf1_idx));

% Beta matrix: events x channels
betas_bf1 = SPM.beta(bf1_idx, :);

% Clean event names
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

out_all = fullfile(output_dir, sprintf('betas_%s_bf1.csv', signal));
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

out_valid = fullfile(output_dir, sprintf('betas_%s_bf1_validch.csv', signal));
writetable(beta_table_valid, out_valid, 'WriteRowNames', true);
fprintf('Saved (valid channels only):\n  %s\n', out_valid);

%% -------------------------------------------------------------------------
%  SAVE: SUMMARY TEXT FILE
% -------------------------------------------------------------------------

out_summary = fullfile(output_dir, 'beta_summary.txt');
fid = fopen(out_summary, 'w');
fprintf(fid, 'Beta Extraction Summary\n');
fprintf(fid, '=======================\n');
fprintf(fid, 'Participant : %s\n', ptp);
fprintf(fid, 'Joint       : %s\n', joint_str);
fprintf(fid, 'Tech        : %s\n', tech_str);
fprintf(fid, 'Signal      : %s\n', signal);
fprintf(fid, 'Date        : %s\n\n', datestr(now));
fprintf(fid, 'Total channels    : %d\n', size(SPM.beta,2));
fprintf(fid, 'Valid channels    : %d\n', length(valid_ch));
fprintf(fid, 'Valid channel IDs : %s\n\n', num2str(valid_labels));
fprintf(fid, 'Events extracted (%d):\n', length(clean_events));
for i = 1:length(clean_events)
    fprintf(fid, '  %2d. %s\n', i, clean_events{i});
end
fprintf(fid, '\nNuisance regressors (not extracted): hr, br, rvt, constant\n');
fclose(fid);
fprintf('Saved summary:\n  %s\n', out_summary);

%% -------------------------------------------------------------------------
%  DISPLAY PREVIEW
% -------------------------------------------------------------------------

fprintf('\n--- Beta preview (valid channels, first 5 events) ---\n');
disp(beta_table_valid(1:min(5,end), 1:min(8,end)));

fprintf('\nDone. Files saved to:\n  %s\n', output_dir);