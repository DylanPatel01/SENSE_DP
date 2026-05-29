%% analysis2_physiology_THA.m
% Analysis 2 (Physiology): Within-surgeon comparison of RA-THA vs C-THA
%
% Only includes surgeons who have performed BOTH techniques:
%   Surgeons with both RA-THA and C-THA: 1, 5, 6 (numbers from master file)
%
% Design:
%   - Within-surgeon comparison (individual differences controlled)
%   - Uses physiology_by_event.csv (per-event output) to compute
%     group-level means per surgeon x tech, matching physiology_by_group.m
%     grouping definitions
%   - For surgeons with multiple acquisitions per tech: average across
%     acquisitions first so each surgeon contributes one value per tech
%   - Effect size: Hedges' g ONLY — no inferential stats (underpowered)
%   - Positive difference = RA-THA higher than C-THA
%
% Stage group definitions mirror physiology_by_group.m (Hip):
%   shared_same, shared_diff, robotic_only, conventional_only
%
% Input:  physiology_by_event.csv files (from physiology_by_event.m)
% Output: analysis2_physiology_within_surgeon.csv
%
% Dylan Patel | SENSE Project

clear all; close all; clc;

%% -------------------------------------------------------------------------
%  CONFIGURATION
%% -------------------------------------------------------------------------

data_base  = '/Volumes/Dylan SSD/DYLAN/Data';
output_dir = '/Volumes/Dylan SSD/DYLAN/Results/Analysis';
if ~isfolder(output_dir), mkdir(output_dir); end

% Surgeons with BOTH RA-THA and C-THA
both_tech_surgeons = ["1", "5", "6"];

metrics = {'HR', 'BR', 'RVT', 'RMSSD'};

%% -------------------------------------------------------------------------
%  STAGE GROUP DEFINITIONS (Hip) — mirrors physiology_by_group.m
%% -------------------------------------------------------------------------

hip_groups = struct();
hip_groups.shared_same = {
    'approach'
    'dislocation'
    'neck_osteotomy'
    'acetabular_exposure'
    'acetabular_prep'
    'cup_impaction'
    'liner_insertion'
    'femoral_prep'
    'stem_trial'
    'stem_impaction'
    'reduction'
    'closure'
};
hip_groups.shared_diff = {
    'cup_trial'
    'drill_insert_screws'
    'femoral_sizing'
    'trialling'
};
hip_groups.conventional_only = {
    'freehand_cup'
};
hip_groups.robotic_only = {
    'insert_mako_pins'
    'mako_setup'
    'femoral_registration'
    'acetabular_registration'
    'remove_mako_pins'
};

group_names = fieldnames(hip_groups);

%% -------------------------------------------------------------------------
%  STAGE NORMALISATION FUNCTION
%  (mirrors normalise_stage from physiology_by_group.m)
%% -------------------------------------------------------------------------

function canon = normalise_stage(raw)
    s = lower(strtrim(raw));
    if contains(s,'glove') || contains(s,'gc ') || contains(s,'ignore') || ...
       contains(s,'problem') || contains(s,'robot stop') || contains(s,'issue with') || ...
       contains(s,'switch to manual') || contains(s,'waiting for cement') || ...
       contains(s,'adjustments') || contains(s,'corrections') || ...
       strcmp(s,'end') || strcmp(s,'end of surgery') || strcmp(s,'start surgery') || ...
       strcmp(s,'a') || contains(s,'back to') || contains(s,'continuation of') || ...
       contains(s,'checking on mako') || contains(s,'injection of local')
        canon = ''; return
    end
    if contains(s,'approach') || contains(s,'incision'),           canon = 'approach'; return; end
    if contains(s,'insert mako') || contains(s,'insertion of mako') || ...
       contains(s,'pin insertion') || contains(s,'applying arrays') || ...
       contains(s,'checkpoint insertion'),                          canon = 'insert_mako_pins'; return; end
    if contains(s,'mako set up') || contains(s,'mako setup'),      canon = 'mako_setup'; return; end
    if contains(s,'removal of mako') || contains(s,'checkpoint removal') || ...
       contains(s,'pins out') || contains(s,'remove mako'),        canon = 'remove_mako_pins'; return; end
    if contains(s,'femoral reg') || contains(s,'femoral verif') || ...
       contains(s,'fem reg') || strcmp(s,'registration') || ...
       strcmp(s,'verification') || strcmp(s,'reg and veri') || ...
       contains(s,'capture centre of rotation'),                    canon = 'femoral_registration'; return; end
    if contains(s,'acetab reg') || contains(s,'acetabular reg') || ...
       contains(s,'placing the acetabular checkpoint'),             canon = 'acetabular_registration'; return; end
    if contains(s,'dislocation'),                                   canon = 'dislocation'; return; end
    if contains(s,'neck osteotomy') || contains(s,'neck cut'),     canon = 'neck_osteotomy'; return; end
    if contains(s,'acetab') && contains(s,'expos'),                canon = 'acetabular_exposure'; return; end
    if contains(s,'femoral exposure'),                             canon = 'acetabular_exposure'; return; end
    if contains(s,'acetab') && (contains(s,'prep') || contains(s,'ream')), canon = 'acetabular_prep'; return; end
    if contains(s,'cup trial'),                                    canon = 'cup_trial'; return; end
    if contains(s,'cup imp') || contains(s,'cup insert'),          canon = 'cup_impaction'; return; end
    if contains(s,'freehand'),                                     canon = 'freehand_cup'; return; end
    if contains(s,'drill') && contains(s,'screw'),                 canon = 'drill_insert_screws'; return; end
    if contains(s,'screw insertion') || contains(s,'insert screws'), canon = 'drill_insert_screws'; return; end
    if contains(s,'liner'),                                        canon = 'liner_insertion'; return; end
    if contains(s,'femoral prep') || contains(s,'fem prep') || ...
       contains(s,'femoral preparation'),                          canon = 'femoral_prep'; return; end
    if contains(s,'femoral siz') || contains(s,'fem siz'),         canon = 'femoral_sizing'; return; end
    if contains(s,'stem trial') || contains(s,'stem trialling'),   canon = 'stem_trial'; return; end
    if contains(s,'trialling') || contains(s,'trial reduction'),   canon = 'trialling'; return; end
    if contains(s,'stem imp') || contains(s,'stem insert'),        canon = 'stem_impaction'; return; end
    if contains(s,'reduction'),                                    canon = 'reduction'; return; end
    if contains(s,'closure') || contains(s,'wound'),               canon = 'closure'; return; end
    canon = s;
end

%% -------------------------------------------------------------------------
%  LOAD ALL physiology_by_event.csv FILES (Hip only)
%% -------------------------------------------------------------------------

fprintf('Searching for physiology_by_event.csv files (Hip)...\n');

csv_files = dir(fullfile(data_base, 'hip', '**', 'physiology_by_event.csv'));

if isempty(csv_files)
    error('No physiology_by_event.csv files found under:\n  %s/hip\nRun physiology_by_event.m first.', data_base);
end

fprintf('Found %d files.\n', length(csv_files));

master = table();
for f = 1:length(csv_files)
    fpath = fullfile(csv_files(f).folder, csv_files(f).name);
    T = readtable(fpath, 'TextType', 'string');
    master = [master; T];
end

% Force string columns
str_cols = {'Acquisition','Surgeon','Joint','Tech','Event'};
for s = 1:length(str_cols)
    if ismember(str_cols{s}, master.Properties.VariableNames)
        master.(str_cols{s}) = string(master.(str_cols{s}));
    end
end

% Filter to surgeons with both techniques
master = master(ismember(master.Surgeon, both_tech_surgeons), :);
fprintf('Rows after filtering to both-tech surgeons: %d\n', height(master));

%% -------------------------------------------------------------------------
%  MAP EVENTS TO STAGE GROUPS
%% -------------------------------------------------------------------------

master.Group = repmat("", height(master), 1);

for i = 1:height(master)
    canon = normalise_stage(char(master.Event(i)));
    if isempty(canon), continue; end

    for g = 1:length(group_names)
        gname = group_names{g};
        if ismember(canon, hip_groups.(gname))
            master.Group(i) = string(gname);
            break;
        end
    end
end

% Remove unmapped events
master = master(master.Group ~= "", :);
fprintf('Rows after group mapping: %d\n', height(master));

%% -------------------------------------------------------------------------
%  AGGREGATE: EVENT -> ACQUISITION -> SURGEON LEVEL
%  Step 1: Mean within each acquisition x group (pool events in same group)
%  Step 2: Mean across acquisitions per surgeon x tech x group
%% -------------------------------------------------------------------------

fprintf('\n--- Aggregating to surgeon level ---\n');

surg_level = table();
all_groups = [group_names; {'shared_same';'shared_diff';'robotic_only';'conventional_only'}];
all_groups = unique(group_names);

for s = 1:length(both_tech_surgeons)
    surg = both_tech_surgeons(s);
    for tk = ["R", "C"]
        for g = 1:length(group_names)
            grp = string(group_names{g});

            % Cross-check: skip robotic_only for conventional and vice versa
            if grp == "robotic_only"      && tk == "C", continue; end
            if grp == "conventional_only" && tk == "R", continue; end

            mask = master.Surgeon == surg & master.Tech == tk & master.Group == grp;
            rows = master(mask, :);
            if isempty(rows), continue; end

            % Unique acquisitions for this surgeon x tech
            acqs = unique(rows.Acquisition);

            % Step 1: mean per acquisition
            acq_means = NaN(length(acqs), length(metrics));
            for a = 1:length(acqs)
                acq_mask = rows.Acquisition == acqs(a);
                acq_rows = rows(acq_mask, :);
                for m = 1:length(metrics)
                    met = string(metrics{m});
                    if met == "RMSSD"
                        acq_means(a, m) = mean(acq_rows.RMSSD, 'omitnan');
                    else
                        col = strcat(char(met), '_mean');
                        acq_means(a, m) = mean(acq_rows.(col), 'omitnan');
                    end
                end
            end

            % Step 2: mean across acquisitions
            for m = 1:length(metrics)
                met = string(metrics{m});
                surg_mean  = mean(acq_means(:, m), 'omitnan');
                n_acq_used = sum(~isnan(acq_means(:, m)));

                new_row = {char(surg), char(tk), char(grp), char(met), ...
                           surg_mean, n_acq_used};
                surg_level = [surg_level; cell2table(new_row, 'VariableNames', ...
                    {'Surgeon','Tech','Group','Metric','Mean','N_Acquisitions'})];
            end
        end
    end
end

% Force string columns
for s = 1:4
    col = surg_level.Properties.VariableNames{s};
    surg_level.(col) = string(surg_level.(col));
end

fprintf('Surgeon-level rows: %d\n', height(surg_level));

%% -------------------------------------------------------------------------
%  COMPUTE HEDGES' G (WITHIN-SURGEON)
%  For each surgeon: d_i = (RA - C) / SD_pooled_across_surgeons
%  Then Hedges' g = d * J
%
%  With only 3 surgeons, we:
%    1. Compute individual difference scores (RA - C) per surgeon
%    2. Compute Hedges' g using the 3 difference scores
%       (pooled SD = SD of differences across surgeons)
%% -------------------------------------------------------------------------

fprintf('\n--- Computing Hedges'' g (within-surgeon) ---\n');

results = table();

shared_conditions = {'shared_same', 'shared_diff'};
unique_conditions = {'robotic_only', 'conventional_only'};
all_conditions    = [shared_conditions, unique_conditions];

% Storage for pooled SD and J per metric x condition (used later for
% individual surgeon Hedges' g)
pooled_sd_store = containers.Map('KeyType','char','ValueType','double');
J_store         = containers.Map('KeyType','char','ValueType','double');

for m = 1:length(metrics)
    met = string(metrics{m});
    for c = 1:length(all_conditions)
        cond     = string(all_conditions{c});
        is_shared = ismember(char(cond), shared_conditions);

        % Collect RA and C values per surgeon
        vals_rob  = NaN(length(both_tech_surgeons), 1);
        vals_conv = NaN(length(both_tech_surgeons), 1);

        for s = 1:length(both_tech_surgeons)
            surg = both_tech_surgeons(s);

            rob_mask  = surg_level.Surgeon == surg & surg_level.Tech == "R" & ...
                        surg_level.Group == cond & surg_level.Metric == met;
            conv_mask = surg_level.Surgeon == surg & surg_level.Tech == "C" & ...
                        surg_level.Group == cond & surg_level.Metric == met;

            if any(rob_mask),  vals_rob(s)  = surg_level.Mean(rob_mask); end
            if any(conv_mask), vals_conv(s) = surg_level.Mean(conv_mask); end
        end

        % Remove surgeons with missing data for this condition
        valid = ~isnan(vals_rob) & ~isnan(vals_conv);
        vr    = vals_rob(valid);
        vc    = vals_conv(valid);
        n_s   = sum(valid);

        % For unique conditions only one tech has data — report descriptively
        if ~is_shared
            if cond == "robotic_only"
                vals_only = vals_rob(~isnan(vals_rob));
            else
                vals_only = vals_conv(~isnan(vals_conv));
            end
            n_only = length(vals_only);

            new_row = {char(met), char(cond), ...
                mean(vals_only,'omitnan'), std(vals_only,'omitnan'), ...
                NaN, NaN, NaN, NaN, NaN, NaN, ...
                "Descriptive only", n_only, 0, "Unique condition"};
            results = [results; cell2table(new_row, 'VariableNames', ...
                {'Metric','Condition', ...
                 'Mean_Robotic','Mean_Conventional', ...
                 'Mean_Difference','SD_Difference', ...
                 'Cohens_d','J_correction','Hedges_g','Abs_g', ...
                 'Interpretation','n_rob','n_conv','g_note'})];
            continue;
        end

        % Shared conditions: compute Hedges' g from difference scores
        diff_val = NaN; sd_diff = NaN; d = NaN; J = NaN; g = NaN;

        if n_s >= 2
            diffs    = vr - vc;
            diff_val = mean(diffs);
            sd_diff  = std(diffs);
            % Use pooled SD of raw values (not differences) for Cohen's d
            df   = n_s*2 - 2;
            sd_p = sqrt(((n_s-1)*std(vr)^2 + (n_s-1)*std(vc)^2) / df);
            if sd_p > 0
                d = (mean(vr) - mean(vc)) / sd_p;
                J = 1 - 3 / (4*df - 1);
                g = d * J;
                % Store for individual surgeon g calculation
                key = char(met + "_" + cond);
                pooled_sd_store(key) = sd_p;
                J_store(key)         = J;
            end
            g_note = sprintf('n=%d surgeons, pooled SD', n_s);
        elseif n_s == 1
            diff_val = vr(1) - vc(1);
            g_note   = 'n=1 surgeon — g not calculable';
        else
            g_note = 'No matched data';
        end

        % Direction
        if ~isnan(g)
            if g > 0,      dir = "RA-THA > C-THA";
            elseif g < 0,  dir = "C-THA > RA-THA";
            else,          dir = "No difference";
            end
            if abs(g) < 0.2,      eff = "Trivial";
            elseif abs(g) < 0.5,  eff = "Small";
            elseif abs(g) < 0.8,  eff = "Medium";
            else,                  eff = "Large";
            end
            interp = string(eff) + " | " + string(dir);
        else
            interp = "Descriptive only";
        end

        new_row = {char(met), char(cond), ...
            mean(vr,'omitnan'), mean(vc,'omitnan'), ...
            diff_val, sd_diff, d, J, g, abs(g), ...
            interp, n_s, n_s, g_note};
        results = [results; cell2table(new_row, 'VariableNames', ...
            {'Metric','Condition', ...
             'Mean_Robotic','Mean_Conventional', ...
             'Mean_Difference','SD_Difference', ...
             'Cohens_d','J_correction','Hedges_g','Abs_g', ...
             'Interpretation','n_rob','n_conv','g_note'})];
    end
end

%% -------------------------------------------------------------------------
%  INDIVIDUAL SURGEON PROFILES (for reporting)
%  Individual Hedges' g = (RA_i - C_i) / SD_pooled_group * J
%  Uses the group-level pooled SD and J correction already computed above.
%  This standardises each surgeon's difference on the same scale as the
%  group effect size, enabling per-surgeon bars on a Hedges' g axis.
%% -------------------------------------------------------------------------

surg_profiles = table();

for s = 1:length(both_tech_surgeons)
    surg = both_tech_surgeons(s);
    for m = 1:length(metrics)
        met = string(metrics{m});
        for c = 1:length(shared_conditions)
            cond = string(shared_conditions{c});

            rob_mask  = surg_level.Surgeon == surg & surg_level.Tech == "R" & ...
                        surg_level.Group == cond & surg_level.Metric == met;
            conv_mask = surg_level.Surgeon == surg & surg_level.Tech == "C" & ...
                        surg_level.Group == cond & surg_level.Metric == met;

            val_r = NaN; val_c = NaN;
            if any(rob_mask),  val_r = surg_level.Mean(rob_mask); end
            if any(conv_mask), val_c = surg_level.Mean(conv_mask); end

            raw_diff = val_r - val_c;

            % Compute individual Hedges' g using group pooled SD and J
            ind_g = NaN;
            key   = char(met + "_" + cond);
            if isKey(pooled_sd_store, key) && ~isnan(raw_diff)
                sd_p_grp = pooled_sd_store(key);
                J_grp    = J_store(key);
                if sd_p_grp > 0
                    ind_g = (raw_diff / sd_p_grp) * J_grp;
                end
            end

            new_row = {char(surg), char(met), char(cond), val_r, val_c, raw_diff, ind_g};
            surg_profiles = [surg_profiles; cell2table(new_row, 'VariableNames', ...
                {'Surgeon','Metric','Condition', ...
                 'Mean_Robotic','Mean_Conventional','Difference','Individual_Hedges_g'})];
        end
    end
end

%% -------------------------------------------------------------------------
%  SAVE RESULTS
%% -------------------------------------------------------------------------

out_main    = fullfile(output_dir, 'analysis2_physiology_within_surgeon.csv');
out_profile = fullfile(output_dir, 'analysis2_physiology_surgeon_profiles.csv');

writetable(results,       out_main);
writetable(surg_profiles, out_profile);

fprintf('\nSaved: %s\n', out_main);
fprintf('Saved: %s\n', out_profile);

%% -------------------------------------------------------------------------
%  PRINT SUMMARY
%% -------------------------------------------------------------------------

fprintf('\n');
fprintf('========================================================================\n');
fprintf('  ANALYSIS 2 (PHYSIOLOGY): Within-Surgeon RA-THA vs C-THA\n');
fprintf('========================================================================\n');
fprintf('  Surgeons with both techniques: %s\n', strjoin(both_tech_surgeons, ', '));
fprintf('  Effect size: Hedges'' g ONLY — no inferential stats (n=%d surgeons)\n', ...
    length(both_tech_surgeons));
fprintf('  Positive g/difference = RA-THA > C-THA\n');
fprintf('  Benchmarks: Trivial<0.2, Small 0.2-0.5, Medium 0.5-0.8, Large>=0.8\n');
fprintf('------------------------------------------------------------------------\n');

% ── Individual profiles ───────────────────────────────────────────────────
fprintf('\n  INDIVIDUAL SURGEON PROFILES (Hedges'' g, standardised by group pooled SD)\n');
for m = 1:length(metrics)
    met = string(metrics{m});
    fprintf('\n  Metric: %s\n', char(met));
    fprintf('  %-22s', 'Condition');
    for s = 1:length(both_tech_surgeons)
        fprintf('  Surg_%s (g)', char(both_tech_surgeons(s)));
    end
    fprintf('\n  %s\n', repmat('-', 1, 80));

    for c = 1:length(shared_conditions)
        cond = string(shared_conditions{c});
        fprintf('  %-22s', char(cond));
        for s = 1:length(both_tech_surgeons)
            surg = both_tech_surgeons(s);
            mask = strcmp(surg_profiles.Surgeon, char(surg)) & ...
                   strcmp(surg_profiles.Metric,  char(met))  & ...
                   strcmp(surg_profiles.Condition, char(cond));
            if any(mask)
                fprintf('  %+14.4f    ', surg_profiles.Individual_Hedges_g(mask));
            else
                fprintf('  %14s    ', 'N/A');
            end
        end
        fprintf('\n');
    end
end

% ── Hedges' g summary ────────────────────────────────────────────────────
fprintf('\n');
fprintf('========================================================================\n');
fprintf('  HEDGES'' g SUMMARY — SHARED CONDITIONS\n');
fprintf('========================================================================\n');
fprintf('  %-10s %-22s %12s %12s %12s %7s %7s %7s  %-25s\n', ...
    'Metric', 'Condition', 'Mean_RA', 'Mean_C', 'Diff', 'd', 'J', 'g', 'Interpretation');
fprintf('  %s\n', repmat('-', 1, 115));

for m = 1:length(metrics)
    met = string(metrics{m});
    for c = 1:length(shared_conditions)
        cond = string(shared_conditions{c});
        mask = strcmp(results.Metric, char(met)) & strcmp(results.Condition, char(cond));
        if ~any(mask), continue; end
        row = results(mask, :);
        fprintf('  %-10s %-22s %12.4f %12.4f %12.4f %7.3f %7.3f %7.3f  %s\n', ...
            char(met), char(cond), ...
            row.Mean_Robotic, row.Mean_Conventional, row.Mean_Difference, ...
            row.Cohens_d, row.J_correction, row.Hedges_g, ...
            char(row.Interpretation));
    end
end

fprintf('\n');
fprintf('========================================================================\n');
fprintf('  UNIQUE CONDITIONS (descriptive only)\n');
fprintf('========================================================================\n');
fprintf('  %-10s %-22s %12s %12s\n', 'Metric', 'Condition', 'Mean', 'SD');
fprintf('  %s\n', repmat('-', 1, 60));
for i = 1:height(results)
    if ~ismember(char(results.Condition(i)), shared_conditions)
        fprintf('  %-10s %-22s %12.4f %12.4f\n', ...
            results.Metric{i}, results.Condition{i}, ...
            results.Mean_Robotic(i), results.SD_Difference(i));
    end
end

fprintf('\n========================================================================\n');
fprintf('  NOTE\n');
fprintf('========================================================================\n');
fprintf('  Hedges'' g computed from pooled SD across %d matched surgeons.\n', ...
    length(both_tech_surgeons));
fprintf('  No inferential statistics — study not powered for within-surgeon testing.\n');
fprintf('  Individual difference scores shown to illustrate consistency of direction.\n');
fprintf('\nDone.\n');