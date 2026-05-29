%% physiology_by_group.m
% Extracts mean, SD, min, max of HR, BR, RVT, and approximate RMSSD
% for surgical events grouped into 4 categories:
%   shared_same, shared_diff, conventional_only, robotic_only
%
% Uses the same normalisation and group definitions as standard_glm_grouped.m
%
% Each group's epoch = all samples belonging to any event in that group
% (onset to next onset for each constituent event, then pooled).
%
% OUTPUT:
%   physiology_by_group.csv        — one row per group, summary stats
%   physiology_by_group_long.csv   — long format for multi-acquisition stats


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
fprintf('Grouped physiology extraction for: %s\n', A.name);

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

results_base = ['Physiology Grouped' filesep A.name];
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

% Truncate all vectors to the shortest length
n_samples = min([length(nirs_time), length(hr), length(br), length(rvt)]);
nirs_time = nirs_time(1:n_samples);
hr        = hr(1:n_samples);
br        = br(1:n_samples);
rvt       = rvt(1:n_samples);

% Sampling rate
fs = 1 / (nirs_time(2) - nirs_time(1));

% Approximate RR intervals from HR (in ms)
rr = 60000 ./ hr;

%% -------------------------------------------------------------------------
%  STAGE NAME NORMALISATION
% -------------------------------------------------------------------------

function canon = normalise_stage(raw)
    s = lower(strtrim(raw));

    % ---- CONTROL / EXCLUDE ----
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
        canon = 'acetabular_exposure'; return
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

    % ---- CEMENT (hip context) ----
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

    % ---- BONE CUTS (generic) ----
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
       contains(s, 'posterior femur cut')
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

    % ---- CEMENT AND IMPACTION (knee) ----
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

    % ---- CATCH-ALL ----
    canon = '';
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
knee_groups.conventional_only = {};
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

if strcmpi(joint_str, 'H')
    groups = hip_groups;
else
    groups = knee_groups;
end

group_names = fieldnames(groups);
n_groups    = length(group_names);

%% -------------------------------------------------------------------------
%  NORMALISE STAGES AND BUILD EVENT EPOCHS
% -------------------------------------------------------------------------
raw_stages = A.stage;
raw_times  = A.stagetime;
n_raw      = length(raw_stages);

% Sort all events chronologically (needed to define epoch boundaries)
[all_sorted_times, all_sort_idx] = sort(raw_times);
all_sorted_names = raw_stages(all_sort_idx);

% Normalise
norm_sorted = cell(1, n_raw);
for i = 1:n_raw
    norm_sorted{i} = normalise_stage(all_sorted_names{i});
end

% Cross-check: exclude wrong-technique stages
robotic_only_list = groups.robotic_only;
conv_only_list    = groups.conventional_only;
for i = 1:n_raw
    if isempty(norm_sorted{i}), continue; end
    if strcmpi(tech_str, 'C') && ismember(norm_sorted{i}, robotic_only_list)
        norm_sorted{i} = '';
    elseif strcmpi(tech_str, 'R') && ismember(norm_sorted{i}, conv_only_list)
        norm_sorted{i} = '';
    end
end

% Build all group stage lists
all_group_stages = {};
for g = 1:n_groups
    all_group_stages = [all_group_stages; groups.(group_names{g})];
end

%% -------------------------------------------------------------------------
%  COMPUTE PHYSIOLOGY PER GROUP
% -------------------------------------------------------------------------
% For each group, collect all sample indices from constituent event epochs,
% then compute summary stats on the pooled samples.

results = struct();
row = 0;

for g = 1:n_groups
    gname       = group_names{g};
    member_list = groups.(gname);

    % Collect all sample indices for this group
    group_idx = [];
    n_events_in_group = 0;
    total_duration = 0;

    for e = 1:n_raw
        if ~ismember(norm_sorted{e}, member_list), continue; end

        n_events_in_group = n_events_in_group + 1;
        t_start = all_sorted_times(e);
        if e < n_raw
            t_end = all_sorted_times(e + 1);
        else
            t_end = nirs_time(end);
        end

        idx = find(nirs_time >= t_start & nirs_time < t_end);
        group_idx = [group_idx; idx(:)];
        total_duration = total_duration + (t_end - t_start);
    end

    if isempty(group_idx)
        fprintf('Group %-25s : NO samples (skipping)\n', gname);
        continue;
    end

    fprintf('Group %-25s : %d events, %d samples, %.0f s total\n', ...
        gname, n_events_in_group, length(group_idx), total_duration);

    row = row + 1;
    results(row).Group          = gname;
    results(row).N_events       = n_events_in_group;
    results(row).Total_duration = total_duration;
    results(row).N_samples      = length(group_idx);

    % HR
    results(row).HR_mean = mean(hr(group_idx), 'omitnan');
    results(row).HR_sd   = std(hr(group_idx), 'omitnan');
    results(row).HR_min  = min(hr(group_idx));
    results(row).HR_max  = max(hr(group_idx));

    % BR
    results(row).BR_mean = mean(br(group_idx), 'omitnan');
    results(row).BR_sd   = std(br(group_idx), 'omitnan');
    results(row).BR_min  = min(br(group_idx));
    results(row).BR_max  = max(br(group_idx));

    % RVT
    results(row).RVT_mean = mean(rvt(group_idx), 'omitnan');
    results(row).RVT_sd   = std(rvt(group_idx), 'omitnan');
    results(row).RVT_min  = min(rvt(group_idx));
    results(row).RVT_max  = max(rvt(group_idx));

    % RMSSD (approximate from HR-derived RR intervals)
    rr_epoch = rr(group_idx);
    rr_diff  = diff(rr_epoch);
    if length(rr_diff) > 1
        results(row).RMSSD = sqrt(mean(rr_diff.^2, 'omitnan'));
    else
        results(row).RMSSD = NaN;
    end
end

%% -------------------------------------------------------------------------
%  BUILD AND SAVE WIDE TABLE
% -------------------------------------------------------------------------
T = struct2table(results);

T.Acquisition = repmat({A.name}, height(T), 1);
T.Surgeon     = repmat({ptp}, height(T), 1);
T.Joint       = repmat({joint_str}, height(T), 1);
T.Tech        = repmat({tech_str}, height(T), 1);

% Reorder metadata first
T = T(:, [end-3:end, 1:end-4]);

out_wide = fullfile(output_dir, 'physiology_by_group.csv');
writetable(T, out_wide);
fprintf('\nSaved (wide): %s\n', out_wide);

%% -------------------------------------------------------------------------
%  BUILD AND SAVE LONG TABLE (for multi-acquisition concatenation)
% -------------------------------------------------------------------------
% One row per group per metric — easier for stats

metrics = {'HR', 'BR', 'RVT', 'RMSSD'};
long_rows = {};
r = 0;

for i = 1:height(T)
    for m = 1:length(metrics)
        r = r + 1;
        long_rows{r, 1} = A.name;
        long_rows{r, 2} = ptp;
        long_rows{r, 3} = joint_str;
        long_rows{r, 4} = tech_str;
        long_rows{r, 5} = T.Group{i};
        long_rows{r, 6} = metrics{m};

        switch metrics{m}
            case 'HR'
                long_rows{r, 7} = T.HR_mean(i);
                long_rows{r, 8} = T.HR_sd(i);
                long_rows{r, 9} = T.HR_min(i);
                long_rows{r,10} = T.HR_max(i);
            case 'BR'
                long_rows{r, 7} = T.BR_mean(i);
                long_rows{r, 8} = T.BR_sd(i);
                long_rows{r, 9} = T.BR_min(i);
                long_rows{r,10} = T.BR_max(i);
            case 'RVT'
                long_rows{r, 7} = T.RVT_mean(i);
                long_rows{r, 8} = T.RVT_sd(i);
                long_rows{r, 9} = T.RVT_min(i);
                long_rows{r,10} = T.RVT_max(i);
            case 'RMSSD'
                long_rows{r, 7} = T.RMSSD(i);
                long_rows{r, 8} = NaN;
                long_rows{r, 9} = NaN;
                long_rows{r,10} = NaN;
        end
    end
end

long_table = cell2table(long_rows, 'VariableNames', ...
    {'Acquisition', 'Surgeon', 'Joint', 'Tech', 'Group', 'Metric', ...
     'Mean', 'SD', 'Min', 'Max'});

out_long = fullfile(output_dir, 'physiology_by_group_long.csv');
writetable(long_table, out_long);
fprintf('Saved (long): %s\n', out_long);

%% -------------------------------------------------------------------------
%  DISPLAY PREVIEW
% -------------------------------------------------------------------------
fprintf('\n--- Grouped physiology summary ---\n');
disp(T(:, {'Group', 'N_events', 'Total_duration', 'HR_mean', 'BR_mean', 'RMSSD'}));

fprintf('\nDone. Output saved to:\n  %s\n', output_dir);