%% physiology_by_event.m
% Extracts mean, SD, min, max of HR, BR, RVT, and approximate RMSSD
% for each individual surgical event (onset to next onset).
% Saves results as a CSV spreadsheet.
%
% USAGE: Set acq_idx and run. No GLM needed — works directly from
%        the acquisitions struct.
%
% OUTPUT:
%   physiology_by_event.csv   — one row per event, columns for each metric

clear all; close all; clc;

%% -------------------------------------------------------------------------
%  INPUT PARAMETERS
% -------------------------------------------------------------------------
acq_path = '/Volumes/Dylan SSD/acquisitions.mat';
load(acq_path, 'acquisitions');

acq_idx = 28;  % <-- CHANGE THIS

% --- Alternative: select by name ---
% acq_name = '1_H_C_101025_1';
% acq_idx  = find(strcmp({acquisitions.name}, acq_name));
% if isempty(acq_idx), error('Acquisition "%s" not found.', acq_name); end

A = acquisitions(acq_idx);
fprintf('Physiology extraction for: %s\n', A.name);

% Parse acquisition name
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

results_base = ['Physiology Events' filesep A.name];
datafold     = fullfile('/Volumes/Dylan SSD/DYLAN/Data', joint_folder, tech_folder, ptp);
output_dir   = fullfile(datafold, results_base);
if ~isfolder(output_dir), mkdir(output_dir); end

%% -------------------------------------------------------------------------
%  EXTRACT PHYSIOLOGY AND TIME VECTORS
% -------------------------------------------------------------------------
hr  = A.heartrate(:);
br  = A.breathingrate(:);
rvt = A.respvoltime(:);
nirs_time = A.nirstime(:);

% Truncate all vectors to the shortest length (physiology may be slightly
% shorter than NIRS time if trimmed to overlap window during preprocessing)
n_samples = min([length(nirs_time), length(hr), length(br), length(rvt)]);
nirs_time = nirs_time(1:n_samples);
hr        = hr(1:n_samples);
br        = br(1:n_samples);
rvt       = rvt(1:n_samples);

% Sampling rate from time vector
fs = 1 / (nirs_time(2) - nirs_time(1));

% Approximate RR intervals from HR (in ms)
rr = 60000 ./ hr;

%% -------------------------------------------------------------------------
%  SORT EVENTS BY ONSET TIME
% -------------------------------------------------------------------------
[sorted_times, sort_idx] = sort(A.stagetime);
sorted_names = A.stage(sort_idx);
n_events = length(sorted_names);

%% -------------------------------------------------------------------------
%  COMPUTE PHYSIOLOGY METRICS PER EVENT
% -------------------------------------------------------------------------
% Each event epoch = onset to next onset (last event = onset to end of recording)

results = struct();
for e = 1:n_events
    t_start = sorted_times(e);
    if e < n_events
        t_end = sorted_times(e + 1);
    else
        t_end = nirs_time(end);
    end

    % Find sample indices for this epoch
    idx = find(nirs_time >= t_start & nirs_time < t_end);

    if isempty(idx)
        warning('No samples found for event %d: %s (%.1f - %.1f s)', e, sorted_names{e}, t_start, t_end);
        idx = find(nirs_time >= t_start, 1, 'first');  % at least one sample
    end

    % Duration
    results(e).Event      = sorted_names{e};
    results(e).Onset_s    = t_start;
    results(e).Duration_s = t_end - t_start;
    results(e).N_samples  = length(idx);

    % HR
    results(e).HR_mean = mean(hr(idx), 'omitnan');
    results(e).HR_sd   = std(hr(idx), 'omitnan');
    results(e).HR_min  = min(hr(idx));
    results(e).HR_max  = max(hr(idx));

    % BR
    results(e).BR_mean = mean(br(idx), 'omitnan');
    results(e).BR_sd   = std(br(idx), 'omitnan');
    results(e).BR_min  = min(br(idx));
    results(e).BR_max  = max(br(idx));

    % RVT
    results(e).RVT_mean = mean(rvt(idx), 'omitnan');
    results(e).RVT_sd   = std(rvt(idx), 'omitnan');
    results(e).RVT_min  = min(rvt(idx));
    results(e).RVT_max  = max(rvt(idx));

    % RMSSD (approximate from HR-derived RR intervals)
    rr_epoch = rr(idx);
    rr_diff  = diff(rr_epoch);
    if length(rr_diff) > 1
        results(e).RMSSD = sqrt(mean(rr_diff.^2, 'omitnan'));
    else
        results(e).RMSSD = NaN;
    end
end

%% -------------------------------------------------------------------------
%  BUILD AND SAVE TABLE
% -------------------------------------------------------------------------
T = struct2table(results);

% Add metadata columns
T.Acquisition = repmat({A.name}, height(T), 1);
T.Surgeon     = repmat({ptp}, height(T), 1);
T.Joint       = repmat({joint_str}, height(T), 1);
T.Tech        = repmat({tech_str}, height(T), 1);

% Reorder so metadata comes first
T = T(:, [end-3:end, 1:end-4]);

out_file = fullfile(output_dir, 'physiology_by_event.csv');
writetable(T, out_file);
fprintf('\nSaved: %s\n', out_file);

%% -------------------------------------------------------------------------
%  DISPLAY PREVIEW
% -------------------------------------------------------------------------
fprintf('\n--- Preview (first 8 events) ---\n');
disp(T(1:min(8, end), {'Event', 'Duration_s', 'HR_mean', 'BR_mean', 'RMSSD'}));

fprintf('\nDone. Output saved to:\n  %s\n', output_dir);