%% standard_glm_grouped.m
% GLM analysis with surgical stages collapsed into 4 groups:
%   1. shared_same           — stages common to C and R, operationally identical
%   2. shared_diff           — stages common to C and R, operationally different
%   3. conventional_only     — stages unique to conventional
%   4. robotic_only          — stages unique to robotic
%
% Separate group definitions for hip and knee.
% Stage names are normalised from the raw acquisitions struct.

addpath(genpath('/Volumes/Dylan SSD/DYLAN/Toolboxes/spm12_OK'));
addpath(genpath('/Volumes/Dylan SSD/DYLAN/Toolboxes/spm_fnirs_OK'));
addpath(genpath('/Volumes/Dylan SSD/DYLAN/Toolboxes/Homer3-master'))
addpath('/Volumes/Dylan SSD/SENSE/Code');

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
fprintf('Running grouped GLM for: %s\n', A.name);

% Parse joint and tech
name_parts = strsplit(A.name, '_');
joint_str  = upper(name_parts{2});  % 'H' or 'K'
tech_str   = upper(name_parts{3});  % 'C' or 'R'
ptp        = name_parts{1};

if strcmpi(joint_str, 'H'),     joint_folder = 'hip';
elseif strcmpi(joint_str, 'K'), joint_folder = 'knee';
else, error('Unknown joint code: %s', joint_str); end

if strcmpi(tech_str, 'C'),     tech_folder = 'c';
elseif strcmpi(tech_str, 'R'), tech_folder = 'r';
else, error('Unknown tech code: %s', tech_str); end

nirs_abs_start = A.nirs_abs_start;
phys_abs_start = A.phys_abs_start;

% GLM parameters
signal         = 'HbDiff';
d_sample       = 1;
add_regs       = 1;
results_folder = ['GLM TDD Phys Grouped' filesep A.name];
phys_fs        = 2000;

if d_sample == 1
    new_fs = 1;
else
    new_fs = [];
end

%% -------------------------------------------------------------------------
%  DIRECTORY SETUP
% -------------------------------------------------------------------------
datafold    = fullfile('/Volumes/Dylan SSD/DYLAN/Data', joint_folder, tech_folder, ptp);
ptpdatafold = datafold;
disp(ptpdatafold)

%% -------------------------------------------------------------------------
%  STAGE NAME NORMALISATION
% -------------------------------------------------------------------------
% Maps every raw stage name variant to a canonical name.
% Returns empty '' for stages to EXCLUDE (notes, corrections, control stages
% that don't map to a surgical step).

function canon = normalise_stage(raw)
    s = lower(strtrim(raw));

    % ---- CONTROL / EXCLUDE ----
    % Glove changes, notes, corrections, robot issues, misc
    if contains(s, 'glove') || contains(s, 'gc ') || ...
       contains(s, 'GLOVES', 'IgnoreCase', true)
        canon = 'glove_change'; return
    end
    if contains(s, 'ignore') || contains(s, 'last one was') || ...
       contains(s, 'last is') || contains(s, 'problem') || ...
       contains(s, 'problen') || contains(s, 'robot stopped') || ...
       contains(s, 'robot not working') || contains(s, 'issue with') || ...
       contains(s, 'surgeron stopped') || contains(s, 'ricci take over') || ...
       contains(s, 'ricci box cuts') || contains(s, 'the gap before') || ...
       contains(s, 'struggled to verify') || contains(s, 'switch to manual') || ...
       contains(s, 'robot problems') || contains(s, 'robot fixed') || ...
       contains(s, 'waiting for cement') || contains(s, 'redo fem') || ...
       contains(s, 'checkpoint reg again') || ...
       contains(s, 'adjustments') || contains(s, 'corrections on mako') || ...
       contains(s, 'corrections and deformities') || ...
       strcmp(s, 'a') || strcmp(s, 'corrections') || ...
       contains(s, 'back to') || contains(s, 'continuation of') || ...
       contains(s, 'continue with') || ...
       strcmp(s, 'brite') || ...
       contains(s, 'finish setting up') || ...
       contains(s, 'checking on mako') || ...
       contains(s, 'injection of local')
        canon = ''; return
    end

    % ---- END / START markers ----
    if strcmp(s, 'end') || strcmp(s, 'end of surgery') || ...
       strcmp(s, 'end surgery') || strcmp(s, 'start surgery')
        canon = ''; return
    end

    % ---- APPROACH ----
    if contains(s, 'approach') || contains(s, 'incision') || ...
       contains(s, 'skin incision')
        canon = 'approach'; return
    end

    % ---- MAKO PINS / ARRAYS ----
    if contains(s, 'insert mako') || contains(s, 'insertion of mako') || ...
       contains(s, 'insert pins') || contains(s, 'insertion of pins') || ...
       contains(s, 'pin insertion') || contains(s, 'applying arrays') || ...
       contains(s, 'checkpoint insertion') || ...
       contains(s, 'screw/checkpoint insertion')
        canon = 'insert_mako_pins'; return
    end

    % ---- MAKO SETUP ----
    if contains(s, 'mako set up') || contains(s, 'mako setup')
        canon = 'mako_setup'; return
    end

    % ---- REMOVE MAKO PINS / CHECKPOINTS ----
    if contains(s, 'removal of mako') || contains(s, 'removal mako') || ...
       contains(s, 'removal of checkpoint') || contains(s, 'checkpoint removal') || ...
       contains(s, 'checkpoint out') || contains(s, 'pins out') || ...
       contains(s, 'removal of pins') || contains(s, 'removal of rest of mako') || ...
       contains(s, 'remove mako') || ...
       contains(s, 'acetabular checkpoint removal')
        canon = 'remove_mako_pins'; return
    end

    % ---- FEMORAL REGISTRATION / VERIFICATION ----
    if contains(s, 'femoral reg') || contains(s, 'femoral verification') || ...
       contains(s, 'fem reg') || ...
       (contains(s, 'registration') && contains(s, 'verification') && ~contains(s, 'acetab') && ~contains(s, 'tibial')) || ...
       (strcmp(s, 'registration and verification')) || ...
       (strcmp(s, 'registration')) || ...
       (strcmp(s, 'reg and veri')) || ...
       (strcmp(s, 'verification')) || ...
       contains(s, 'capture centre of rotation') || ...
       contains(s, 'insertion of femoral checkpoints')
        canon = 'femoral_registration'; return
    end

    % ---- ACETABULAR REGISTRATION / VERIFICATION ----
    if contains(s, 'acetab reg') || contains(s, 'acetabulum reg') || ...
       contains(s, 'acetabular reg') || ...
       contains(s, 'insertion of acetabular checkpoints') || ...
       contains(s, 'placing the acetabular checkpoint')
        canon = 'acetabular_registration'; return
    end

    % ---- TIBIAL REGISTRATION / VERIFICATION ----
    if contains(s, 'tibial reg') || contains(s, 'tibial verification')
        canon = 'tibial_registration'; return
    end

    % ---- DISLOCATION ----
    if contains(s, 'dislocation')
        canon = 'dislocation'; return
    end

    % ---- NECK OSTEOTOMY / MARKING ----
    if contains(s, 'neck osteotomy') || contains(s, 'neck cut') || ...
       contains(s, 'marking neck') || contains(s, 'performing neck')
        canon = 'neck_osteotomy'; return
    end

    % ---- ACETABULUM / ACETABULAR EXPOSURE ----
    if contains(s, 'acetab') && contains(s, 'expos')
        canon = 'acetabular_exposure'; return
    end
    if contains(s, 'femoral exposure')
        canon = 'acetabular_exposure'; return  % likely same operative step context
    end

    % ---- ACETABULAR PREP ----
    if contains(s, 'acetab') && (contains(s, 'prep') || contains(s, 'ream'))
        canon = 'acetabular_prep'; return
    end

    % ---- CUP IMPACTION / INSERTION / TRIAL ----
    if contains(s, 'cup trial') || contains(s, 'acetabular cup trial')
        canon = 'cup_trial'; return
    end
    if contains(s, 'cup imp') || contains(s, 'cup insertion') || ...
       contains(s, 'cup insert')
        canon = 'cup_impaction'; return
    end

    % ---- DRILL / INSERT SCREWS ----
    if contains(s, 'drill') && (contains(s, 'screw') || contains(s, 'insert'))
        canon = 'drill_insert_screws'; return
    end
    if contains(s, 'screw insertion') || contains(s, 'screw insert') || ...
       contains(s, 'insert screws')
        canon = 'drill_insert_screws'; return
    end

    % ---- LINER INSERTION ----
    if contains(s, 'liner')
        canon = 'liner_insertion'; return
    end

    % ---- FEMORAL PREP ----
    if contains(s, 'femoral prep') || contains(s, 'fem prep') || ...
       contains(s, 'drying the femoral') || contains(s, 'wash') || ...
       contains(s, 'femoral preparation') || ...
       (contains(s, 'back to femur'))
        canon = 'femoral_prep'; return
    end

    % ---- STEM TRIAL / REDUCTION ----
    if contains(s, 'stem trial') || contains(s, 'stem trialling') || ...
       contains(s, 'trial reduction')
        canon = 'stem_trial_reduction'; return
    end

    % ---- STEM INSERTION ----
    if contains(s, 'stem insert') || contains(s, 'insert stem') || ...
       contains(s, 'real stem')
        canon = 'stem_insertion'; return
    end

    % ---- HEAD IMPACTION ----
    if contains(s, 'head impaction')
        canon = 'head_impaction'; return
    end

    % ---- REDUCTION (final) ----
    if contains(s, 'reduction') && ~contains(s, 'trial') && ~contains(s, 'stem')
        canon = 'reduction'; return
    end

    % ---- IMPLANT CHECK ----
    if contains(s, 'implant check')
        canon = 'implant_check'; return
    end

    % ---- CAPSULAR REPAIR ----
    if contains(s, 'capsul') || contains(s, 'repair of the')
        canon = 'capsular_repair'; return
    end

    % ---- CLOSURE ----
    if contains(s, 'closur') || contains(s, 'closing')
        canon = 'closure'; return
    end

    % ---- CEMENT (hip context: stem cementation) ----
    if contains(s, 'cement') && ~contains(s, 'tibial') && ~contains(s, 'impaction of implant')
        canon = 'cementation'; return
    end

    % ---- ASSESSMENT OF STABILITY ----
    if contains(s, 'stability') || contains(s, 'atbility')
        canon = 'stability_assessment'; return
    end

    % ---- ASSESSMENT OF LEG LENGTH ----
    if contains(s, 'assessment of leg length') || contains(s, 'final check') || ...
       contains(s, 'assessing stem')
        canon = 'final_check'; return
    end

    % ===================== KNEE-SPECIFIC =====================

    % ---- OSTEOPHYTE REMOVAL ----
    if contains(s, 'osteoph') || contains(s, 'osteo removal') || ...
       contains(s, 'osteph') || contains(s, 'osetoph') || ...
       contains(s, 'removal of soft tissue')
        canon = 'osteophyte_removal'; return
    end

    % ---- ASSESSMENT OF DEFORMITY ----
    if contains(s, 'deform') || contains(s, 'assessment of deformity') || ...
       contains(s, 'assessment of deformities')
        canon = 'assessment_of_deformity'; return
    end

    % ---- DISTAL FEMORAL CUT ----
    if contains(s, 'distal fem') || ...
       (contains(s, 'bone cuts') && contains(s, 'distal'))
        canon = 'distal_femoral_cut'; return
    end

    % ---- PROXIMAL TIBIAL CUT ----
    if contains(s, 'proximal tib') || contains(s, 'tibia cut') || ...
       contains(s, 'tibial cut') || contains(s, 'tibia prep') || ...
       contains(s, 'sawing of cuts')
        canon = 'proximal_tibial_cut'; return
    end

    % ---- BONE CUTS (generic — often distal femoral + proximal tibial) ----
    if contains(s, 'bone cut') || contains(s, 'fem bone cut') || ...
       contains(s, 'end of cut') || contains(s, 'cuts end') || ...
       contains(s, 'completion of cuts') || contains(s, 'cuts complete') || ...
       contains(s, 'completion of rest') || contains(s, 'continuation of cuts') || ...
       contains(s, 'completion of femoral') || contains(s, 'completion/end')
        canon = 'bone_cuts'; return
    end

    % ---- MENISCAL RESECTION ----
    if contains(s, 'menisc')
        canon = 'meniscal_resection'; return
    end

    % ---- FLEXION / EXTENSION GAPS ----
    if contains(s, 'flexion') || contains(s, 'extension') || ...
       contains(s, 'gap balancing') || contains(s, 'assessment of gaps') || ...
       contains(s, 'assessment of flexion')
        canon = 'gap_assessment'; return
    end

    % ---- FEMORAL SIZING ----
    if contains(s, 'femoral sizing')
        canon = 'femoral_sizing'; return
    end

    % ---- BOX CUTS / CHAMFER ----
    if contains(s, 'box cut') || contains(s, 'chamfer') || ...
       contains(s, 'posterior femur cut') || contains(s, 'posterior femur cut')
        canon = 'box_cuts'; return
    end

    % ---- POSTERIOR CLEARANCE ----
    if contains(s, 'posterior clearance') || contains(s, 'posterior clearence')
        canon = 'posterior_clearance'; return
    end

    % ---- TRIALLING (knee) ----
    if contains(s, 'trialling') || contains(s, 'trailling') || ...
       contains(s, 'trial') || contains(s, 'tibial trial')
        canon = 'trialling'; return
    end

    % ---- PATELLA PREPARATION ----
    if contains(s, 'patella')
        canon = 'patella_prep'; return
    end

    % ---- TIBIAL PREPARATION ----
    if contains(s, 'tibial prep') || contains(s, 'tibia prep') || ...
       contains(s, 'preparation/insertion of tibial')
        canon = 'tibial_prep'; return
    end

    % ---- CEMENT AND IMPACTION OF IMPLANTS (knee) ----
    if contains(s, 'cement') && (contains(s, 'impaction of implant') || ...
       contains(s, 'impact implant') || contains(s, 'impaction'))
        canon = 'cement_and_impaction'; return
    end
    if contains(s, 'cement') || contains(s, 'cementing')
        canon = 'cement_and_impaction'; return
    end

    % ---- INSERTION OF TIBIAL LINER ----
    if contains(s, 'tibial liner') || contains(s, 'insert liner')
        canon = 'insert_tibial_liner'; return
    end

    % ---- FEMORAL ANTIVERSION ----
    if contains(s, 'antiversion')
        canon = 'assessment_femoral_antiversion'; return
    end

    % ---- CATCH-ALL: unmatched ----
    canon = '';  % exclude unmatched stages
end

%% -------------------------------------------------------------------------
%  GROUP DEFINITIONS
% -------------------------------------------------------------------------

% HIP groups
hip_groups = struct();
hip_groups.shared_same = {
    'approach'
    'drill_insert_screws'
    'liner_insertion'
    'stem_insertion'
    'reduction'
    'capsular_repair'
    'closure'
};
hip_groups.shared_diff = {
    'femoral_prep'
    'acetabular_prep'
    'acetabular_exposure'
    'cup_impaction'
    'stem_trial_reduction'
};
hip_groups.conventional_only = {
    'dislocation'
    'neck_osteotomy'
    'cementation'
    'head_impaction'
};
hip_groups.robotic_only = {
    'insert_mako_pins'
    'mako_setup'
    'femoral_registration'
    'acetabular_registration'
    'remove_mako_pins'
};

% KNEE groups
% NOTE: All knee acquisitions in this dataset are robotic (K_R).
% Stages are still grouped by their theoretical classification for when
% conventional knee data is collected. For now, conventional_only will
% produce no onsets.
knee_groups = struct();
knee_groups.shared_same = {
    'approach'
    'meniscal_resection'
    'gap_assessment'
    'trialling'
    'patella_prep'
    'tibial_prep'
    'posterior_clearance'
    'cement_and_impaction'
    'insert_tibial_liner'
    'closure'
};
knee_groups.shared_diff = {
    'distal_femoral_cut'
    'proximal_tibial_cut'
    'box_cuts'
    'stability_assessment'
    'bone_cuts'
};
knee_groups.conventional_only = {
    % Currently no conventional knee acquisitions in dataset
};
knee_groups.robotic_only = {
    'insert_mako_pins'
    'mako_setup'
    'femoral_registration'
    'tibial_registration'
    'femoral_sizing'
    'assessment_of_deformity'
    'osteophyte_removal'
    'remove_mako_pins'
};

% Select group definitions based on joint
if strcmpi(joint_str, 'H')
    groups = hip_groups;
else
    groups = knee_groups;
end

group_names  = fieldnames(groups);
n_groups     = length(group_names);

%% -------------------------------------------------------------------------
%  EXTRACT NIRS DATA FROM ACQUISITIONS STRUCT
% -------------------------------------------------------------------------
Y.hbo = A.data(:, :, 1);
Y.hbr = A.data(:, :, 2);
Y.hbt = A.data(:, :, 3);
Y.hbd = Y.hbo - Y.hbr;
Y.od  = [];  % placeholder for spmfnirsflow

nirs_time = A.nirstime(:);

if ~isempty(A.badchannels)
    exclude_channels = A.badchannels;
else
    disp('No bad channels found, continuing with all channels')
    exclude_channels = [];
end

nirs_abs_t   = nirs_abs_start + seconds(nirs_time);
nirs_abs_end = nirs_abs_t(end);
nirs_fs      = 1 / (nirs_time(2) - nirs_time(1));

%% -------------------------------------------------------------------------
%  NORMALISE STAGES AND ASSIGN TO GROUPS
% -------------------------------------------------------------------------
raw_stages  = A.stage;
raw_times   = A.stagetime;
n_raw       = length(raw_stages);

% Normalise all stage names
norm_stages = cell(1, n_raw);
for i = 1:n_raw
    norm_stages{i} = normalise_stage(raw_stages{i});
end

% Cross-check: exclude stages that belong to the wrong technique
% e.g., if a conventional case has a mislabelled robotic-only stage
robotic_only_list = groups.robotic_only;
conv_only_list    = groups.conventional_only;

for i = 1:n_raw
    if isempty(norm_stages{i}), continue; end
    if strcmpi(tech_str, 'C') && ismember(norm_stages{i}, robotic_only_list)
        fprintf('  Excluding "%s" (-> %s): robotic-only stage in conventional case\n', ...
            raw_stages{i}, norm_stages{i});
        norm_stages{i} = '';
    elseif strcmpi(tech_str, 'R') && ismember(norm_stages{i}, conv_only_list)
        fprintf('  Excluding "%s" (-> %s): conventional-only stage in robotic case\n', ...
            raw_stages{i}, norm_stages{i});
        norm_stages{i} = '';
    end
end

% Report unmatched stages
unmatched = raw_stages(cellfun(@isempty, norm_stages));
if ~isempty(unmatched)
    fprintf('\n--- Excluded/unmatched stages (%d) ---\n', length(unmatched));
    for i = 1:length(unmatched)
        fprintf('  - %s\n', unmatched{i});
    end
end

% Also report stages that normalised but don't appear in any group
all_group_stages = {};
for g = 1:n_groups
    all_group_stages = [all_group_stages; groups.(group_names{g})];
end

matched_but_ungrouped = {};
for i = 1:n_raw
    if ~isempty(norm_stages{i}) && ~ismember(norm_stages{i}, all_group_stages)
        matched_but_ungrouped{end+1} = sprintf('%s -> %s', raw_stages{i}, norm_stages{i});
    end
end
if ~isempty(matched_but_ungrouped)
    fprintf('\n--- Normalised but not in any group (excluded from GLM) ---\n');
    for i = 1:length(matched_but_ungrouped)
        fprintf('  - %s\n', matched_but_ungrouped{i});
    end
end

% Build grouped onsets: one onset per group, merging all constituent stage times
names     = {};
onsets    = {};
durations = {};
for g = 1:n_groups
    gname       = group_names{g};
    member_list = groups.(gname);

    % Find all normalised stages that belong to this group
    group_times = [];
    for i = 1:n_raw
        if ismember(norm_stages{i}, member_list)
            group_times(end+1) = raw_times(i);
        end
    end

    if ~isempty(group_times)
        names{end+1}     = gname;
        onsets{end+1}    = sort(group_times);  % multiple onsets for this condition
        durations{end+1} = zeros(size(group_times));  % instantaneous
        fprintf('Group %-25s : %d onsets\n', gname, length(group_times));
    else
        fprintf('Group %-25s : NO onsets found (skipping)\n', gname);
    end
end

n_times = length(names);

% Ensure row vectors — spmfnirsflow uses size(names, 2) to count conditions
names     = names(:)';
onsets    = onsets(:)';
durations = durations(:)';

Params = cell(1, n_times);
for x = 1:n_times
    Params{x}.Pname = 'none';
    Params{x}.h     = 0;
    Params{x}.P     = [];
end

% Save onsets file to acquisition-specific results folder
results_path = [ptpdatafold filesep results_folder];
if ~isfolder(results_path), mkdir(results_path); end
save([results_path filesep 'onsets_grouped.mat'], 'names', 'durations', 'onsets', 'Params', '-mat')
ons_file = [results_path filesep 'onsets_grouped'];

%% -------------------------------------------------------------------------
%  EVENT PLOT
% -------------------------------------------------------------------------
colours = lines(n_groups);
figure; hold on
plot(nirs_time, Y.hbo(:,1), 'Color', [0.7 0.7 0.7]);
for g = 1:n_times
    for t = 1:length(onsets{g})
        xline(onsets{g}(t), '-', names{g}, 'Color', colours(g,:), ...
            'LabelOrientation', 'horizontal', 'FontSize', 7);
    end
end
title(sprintf('Grouped onsets: %s', A.name), 'Interpreter', 'none')
hold off

%% -------------------------------------------------------------------------
%  PHYSIOLOGY REGRESSORS
% -------------------------------------------------------------------------
hr  = A.heartrate(:);
br  = A.breathingrate(:);
rvt = A.respvoltime(:);

if d_sample == 1
    newfs = new_fs;
else
    newfs = nirs_fs;
end

ds_factor   = round(nirs_fs / newfs);
nirs_ds_len = round(size(A.data, 1) / ds_factor);

if ds_factor > 1
    hr  = decimate(hr,  ds_factor);
    br  = decimate(br,  ds_factor);
    rvt = decimate(rvt, ds_factor);
end

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

fprintf('\n=== Grouped GLM complete for %s ===\n', A.name);