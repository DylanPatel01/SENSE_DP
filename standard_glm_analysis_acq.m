%Script to conduct the basic GLM analysis for an event-related
%approximation of the surgeons naturalistic study.
%Now reads all participant data from the acquisitions MAT file.
addpath(genpath('/Volumes/Dylan SSD/DYLAN/Toolboxes/spm12_OK'));
addpath(genpath('/Volumes/Dylan SSD/DYLAN/Toolboxes/spm_fnirs_OK'));
addpath(genpath('/Volumes/Dylan SSD/DYLAN/Toolboxes/Homer3-master'))
addpath('/Volumes/Dylan SSD/SENSE/Code');

%% -------------------------------------------------------------------------
%  INPUT PARAMETERS
% -------------------------------------------------------------------------
% Load the acquisitions struct
acq_path = '/Volumes/Dylan SSD/acquisitions.mat';
load(acq_path, 'acquisitions');

% Select acquisition — by index or by name
acq_idx = 2;  % <-- CHANGE THIS to pick participant

% --- Alternative: select by name ---
% acq_name = '1_H_C_101025_1';
% acq_idx  = find(strcmp({acquisitions.name}, acq_name));
% if isempty(acq_idx), error('Acquisition "%s" not found.', acq_name); end

A = acquisitions(acq_idx);  % working copy of this acquisition
fprintf('Running GLM for: %s\n', A.name);

% Parse joint and tech from the name (format: surgeon_joint_tech_date_timebin)
name_parts = strsplit(A.name, '_');
joint_str  = name_parts{2};  % 'H' or 'K'
tech_str   = name_parts{3};  % 'C' or 'R'

if strcmpi(joint_str, 'H'),     joint = 1;
elseif strcmpi(joint_str, 'K'), joint = 2;
else, error('Unknown joint code: %s', joint_str); end

if strcmpi(tech_str, 'C'),     tech = 1;
elseif strcmpi(tech_str, 'R'), tech = 2;
else, error('Unknown tech code: %s', tech_str); end

ptp = name_parts{1};  % surgeon ID used for folder path

% Extract start times from the struct
nirs_abs_start = A.nirs_abs_start;
phys_abs_start = A.phys_abs_start;

% GLM parameters
signal         = 'HbDiff'; % HbDiff, HbO, or HbR
d_sample       = 1;        % 1 = downsample NIRS to new_fs
add_regs       = 1;        % 1 = include ECG/Resp regressors
results_folder = ['GLM TDD Phys' filesep A.name];
phys_fs        = 2000;     % physiology sampling rate

if d_sample == 1
    new_fs = 1;
else
    new_fs = [];
end

%% -------------------------------------------------------------------------
%  DIRECTORY SETUP
% -------------------------------------------------------------------------
datafold = '/Volumes/Dylan SSD/DYLAN/Data';
if joint == 1,     datafold = [datafold filesep 'hip'];
elseif joint == 2, datafold = [datafold filesep 'knee']; end
if tech == 1,      datafold = [datafold filesep 'c'];
elseif tech == 2,  datafold = [datafold filesep 'r']; end

ptpdatafold = [datafold filesep ptp];
disp(ptpdatafold)

%% -------------------------------------------------------------------------
%  EXTRACT NIRS DATA FROM ACQUISITIONS STRUCT
% -------------------------------------------------------------------------
% A.data is samples x channels x 3  (dim3: HbO, HbR, HbT)
Y.hbo = A.data(:, :, 1);
Y.hbr = A.data(:, :, 2);
Y.hbt = A.data(:, :, 3);
Y.hbd = Y.hbo - Y.hbr;
Y.od  = [];  % placeholder — spmfnirsflow expects this field to exist (removes it internally)

% NIRS time vector
nirs_time = A.nirstime(:);  % ensure column vector

% Bad channels
if ~isempty(A.badchannels)
    exclude_channels = A.badchannels;
else
    disp('No bad channels found, continuing with all channels')
    exclude_channels = [];
end

%% -------------------------------------------------------------------------
%  NIRS TIME VECTOR AND ABSOLUTE TIME
% -------------------------------------------------------------------------
nirs_abs_t   = nirs_abs_start + seconds(nirs_time);
nirs_abs_end = nirs_abs_t(end);

% Compute sampling rate from the time vector
nirs_fs = 1 / (nirs_time(2) - nirs_time(1));

%% -------------------------------------------------------------------------
%  ONSETS FROM ACQUISITIONS STRUCT
% -------------------------------------------------------------------------
% A.stage     = 1 x N cell array of event names
% A.stagetime = 1 x N array of onset times (seconds from NIRS start)
[~, sort_idx] = sort(A.stagetime);
names     = A.stage(sort_idx);
onsets    = num2cell(A.stagetime(sort_idx));
n_times   = length(names);
durations = num2cell(zeros(1, n_times));  % instantaneous events; adjust if needed

Params = cell(1, n_times);
for x = 1:n_times
    Params{x}.Pname = 'none';
    Params{x}.h     = 0;
    Params{x}.P     = [];
end

% Save onsets file to acquisition-specific results folder
results_path = [ptpdatafold filesep results_folder];
if ~isfolder(results_path), mkdir(results_path); end
save([results_path filesep 'onsets.mat'], 'names', 'durations', 'onsets', 'Params', '-mat')
ons_file = [results_path filesep 'onsets'];

%% -------------------------------------------------------------------------
%  EVENT PLOT
% -------------------------------------------------------------------------
[sort_times, idx] = sort(cell2mat(onsets));
names_sort = names(idx)';

figure
plot(nirs_time, Y.hbo(:,1));
hold on
xline(sort_times, '-', string(1:n_times))
title(sprintf('Acquisition: %s', A.name), 'Interpreter', 'none')

%% -------------------------------------------------------------------------
%  PHYSIOLOGY REGRESSORS
% -------------------------------------------------------------------------
% Physiology is already resampled to NIRS rate in the acquisitions struct,
% so no alignment or downsampling is needed — use directly.
hr  = A.heartrate(:);
br  = A.breathingrate(:);
rvt = A.respvoltime(:);

% If downsampling NIRS, also downsample the physiology regressors
if d_sample == 1
    newfs = new_fs;
else
    newfs = nirs_fs;
end

ds_factor   = round(nirs_fs / newfs);
nirs_ds_len = round(size(A.data, 1) / ds_factor);

if ds_factor > 1
    % Decimate physiology to match the downsampled NIRS length
    hr  = decimate(hr,  ds_factor);
    br  = decimate(br,  ds_factor);
    rvt = decimate(rvt, ds_factor);
end

% Trim or pad to exactly match downsampled NIRS length
hr  = hr(1:min(end, nirs_ds_len));
br  = br(1:min(end, nirs_ds_len));
rvt = rvt(1:min(end, nirs_ds_len));
if length(hr) < nirs_ds_len
    hr(end+1:nirs_ds_len)  = hr(end);
    br(end+1:nirs_ds_len)  = br(end);
    rvt(end+1:nirs_ds_len) = rvt(end);
end

regressors{1,1}      = zeros(nirs_ds_len, 3);
regressors{1,1}(:,1) = hr;    regressors{2}{1} = 'hr';
regressors{1,1}(:,2) = br;    regressors{2}{2} = 'br';
regressors{1,1}(:,3) = rvt;   regressors{2}{3} = 'rvt';

if add_regs == 0
    regressors = [];
end

%% -------------------------------------------------------------------------
%  BUILD P STRUCT AND SAVE NIRS.MAT
% -------------------------------------------------------------------------
P.ns   = size(Y.hbo, 1);
P.nch  = size(Y.hbo, 2);
P.fs   = nirs_fs;
fprintf('fs = %f\n', P.fs);
P.mask       = ones(1, P.nch);
P.fname.pos  = [ptpdatafold filesep 'Digitisation' filesep 'POS.mat'];
P.fname.nirs = [ptpdatafold filesep results_folder filesep 'NIRS.mat'];
P.fname.hrf  = 'hrf (with time and dispersion derivatives)';

% Delete stale NIRS.mat from any previous run before saving
if isfile(P.fname.nirs), delete(P.fname.nirs); end
if ~isfolder([ptpdatafold filesep results_folder])
    mkdir([ptpdatafold filesep results_folder]);
end
save(P.fname.nirs, 'Y', 'P');

%% -------------------------------------------------------------------------
%  RUN GLM
% -------------------------------------------------------------------------
figure('Visible', 'on');  % ensure a figure exists for spmfnirsflow's getframe
spmfnirsflow(P.fname.nirs, ons_file, results_folder, signal, regressors, d_sample, new_fs, exclude_channels)