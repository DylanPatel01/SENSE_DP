%% analysis3_physiology_RATKA.m
% Analysis 3 (Physiology): Expert vs Novice in RA-TKA
%
% Design:
%   - Independent groups: Expert vs Novice surgeons
%   - RA-TKA only (knee, robotic)
%   - Unit of analysis: surgeon (averaged across acquisitions)
%   - Metrics: HR, BR, RVT, RMSSD
%   - Conditions: shared_same, shared_diff, robotic_only
%   - Effect size: Hedges' g (primary outcome)
%   - Test: Mann-Whitney U (reported with power caveat — n=4 vs n=2)
%
% Surgeon classification (RA-TKA):
%   Expert : 1, 2, 3, 4  (n=4 surgeons)
%   Novice : 5, 7         (n=2 surgeons)
%
% Stage group definitions mirror physiology_by_group.m (Knee):
%   shared_same, shared_diff, robotic_only
%
% Input:  physiology_by_event.csv files (from physiology_by_event.m)
% Output: analysis3_physiology_expertise_RATKA.csv
%         analysis3_physiology_descriptive_RATKA.csv
%
% Dylan Patel | SENSE Project

clear all; close all; clc;

%% -------------------------------------------------------------------------
%  CONFIGURATION
%% -------------------------------------------------------------------------

data_base  = '/Volumes/Dylan SSD/DYLAN/Data';
output_dir = '/Volumes/Dylan SSD/DYLAN/Results/Analysis';
if ~isfolder(output_dir), mkdir(output_dir); end

alpha = 0.05;

expert_ids = ["1", "2", "3", "4"];
novice_ids = ["5", "7"];

metrics = {'HR', 'BR', 'RVT', 'RMSSD'};

%% -------------------------------------------------------------------------
%  STAGE GROUP DEFINITIONS (Knee) — mirrors physiology_by_group.m
%% -------------------------------------------------------------------------

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

group_names       = fieldnames(knee_groups);
shared_conditions = {'shared_same', 'shared_diff'};
unique_conditions = {'robotic_only'};
all_conditions    = [shared_conditions, unique_conditions];

%% -------------------------------------------------------------------------
%  STAGE NORMALISATION (Knee-relevant subset)
%% -------------------------------------------------------------------------

function canon = normalise_stage(raw)
    s = lower(strtrim(raw));
    if contains(s,'glove') || contains(s,'ignore') || contains(s,'problem') || ...
       contains(s,'robot stop') || contains(s,'adjustments') || ...
       contains(s,'corrections') || strcmp(s,'end') || strcmp(s,'start surgery') || ...
       strcmp(s,'a') || contains(s,'back to') || contains(s,'continuation of')
        canon = ''; return
    end
    if contains(s,'approach') || contains(s,'incision'),             canon = 'approach'; return; end
    if contains(s,'insert mako') || contains(s,'pin insertion') || ...
       contains(s,'checkpoint insertion'),                           canon = 'insert_mako_pins'; return; end
    if contains(s,'mako set up') || contains(s,'mako setup'),        canon = 'mako_setup'; return; end
    if contains(s,'removal of mako') || contains(s,'checkpoint removal') || ...
       contains(s,'pins out') || contains(s,'remove mako'),          canon = 'remove_mako_pins'; return; end
    if contains(s,'femoral reg') || contains(s,'fem reg') || ...
       strcmp(s,'registration') || strcmp(s,'verification'),         canon = 'femoral_registration'; return; end
    if contains(s,'tibial reg') || contains(s,'tibial verif'),       canon = 'tibial_registration'; return; end
    if contains(s,'femoral siz') || contains(s,'fem siz'),           canon = 'femoral_sizing'; return; end
    if contains(s,'meniscal'),                                       canon = 'meniscal_resection'; return; end
    if contains(s,'gap assess'),                                     canon = 'gap_assessment'; return; end
    if contains(s,'distal fem') && contains(s,'cut'),                canon = 'distal_femoral_cut'; return; end
    if contains(s,'proximal tib') && contains(s,'cut'),              canon = 'proximal_tibial_cut'; return; end
    if contains(s,'box cut'),                                        canon = 'box_cuts'; return; end
    if contains(s,'bone cut'),                                       canon = 'bone_cuts'; return; end
    if contains(s,'stability'),                                      canon = 'stability_assessment'; return; end
    if contains(s,'osteophyte'),                                     canon = 'osteophyte_removal'; return; end
    if contains(s,'assessment of deformity') || contains(s,'deformity assess'), canon = 'assessment_of_deformity'; return; end
    if contains(s,'patella'),                                        canon = 'patella_prep'; return; end
    if contains(s,'tibial prep') || contains(s,'tib prep'),          canon = 'tibial_prep'; return; end
    if contains(s,'posterior clear'),                                canon = 'posterior_clearance'; return; end
    if contains(s,'trialling') || contains(s,'trial'),               canon = 'trialling'; return; end
    if contains(s,'cement') || contains(s,'impaction'),              canon = 'cement_and_impaction'; return; end
    if contains(s,'tibial liner') || contains(s,'insert liner'),     canon = 'insert_tibial_liner'; return; end
    if contains(s,'closure') || contains(s,'wound'),                 canon = 'closure'; return; end
    canon = s;
end

%% -------------------------------------------------------------------------
%  LOAD ALL physiology_by_event.csv FILES (Knee, Robotic only)
%% -------------------------------------------------------------------------

fprintf('Searching for physiology_by_event.csv files (Knee, Robotic)...\n');

csv_files = dir(fullfile(data_base, 'knee', 'r', '**', 'physiology_by_event.csv'));

if isempty(csv_files)
    error('No physiology_by_event.csv files found under:\n  %s/knee/r\n', data_base);
end

fprintf('Found %d files.\n', length(csv_files));

master = table();
for f = 1:length(csv_files)
    fpath = fullfile(csv_files(f).folder, csv_files(f).name);
    T = readtable(fpath, 'TextType', 'string');
    master = [master; T];
end

str_cols = {'Acquisition','Surgeon','Joint','Tech','Event'};
for s = 1:length(str_cols)
    if ismember(str_cols{s}, master.Properties.VariableNames)
        master.(str_cols{s}) = string(master.(str_cols{s}));
    end
end

% Filter to classified surgeons only
all_ids = [expert_ids, novice_ids];
master  = master(ismember(master.Surgeon, all_ids), :);
fprintf('Rows after surgeon filter: %d\n', height(master));

% Assign expertise
master.Expertise = repmat("Unknown", height(master), 1);
master.Expertise(ismember(master.Surgeon, expert_ids)) = "Expert";
master.Expertise(ismember(master.Surgeon, novice_ids)) = "Novice";

%% -------------------------------------------------------------------------
%  MAP EVENTS TO STAGE GROUPS
%% -------------------------------------------------------------------------

master.Group = repmat("", height(master), 1);

for i = 1:height(master)
    canon = normalise_stage(char(master.Event(i)));
    if isempty(canon), continue; end
    for g = 1:length(group_names)
        gname = group_names{g};
        if ismember(canon, knee_groups.(gname))
            master.Group(i) = string(gname);
            break;
        end
    end
end

master = master(master.Group ~= "", :);
fprintf('Rows after group mapping: %d\n', height(master));

%% -------------------------------------------------------------------------
%  AGGREGATE TO SURGEON LEVEL
%  Step 1: mean within acquisition x group
%  Step 2: mean across acquisitions per surgeon x group
%% -------------------------------------------------------------------------

fprintf('\n--- Aggregating to surgeon level ---\n');

surg_level = table();
all_surgs  = unique(master.Surgeon);

for s = 1:length(all_surgs)
    surg = all_surgs(s);
    exp_grp = char(master.Expertise(find(master.Surgeon == surg, 1)));

    for g = 1:length(group_names)
        grp  = string(group_names{g});
        mask = master.Surgeon == surg & master.Group == grp;
        rows = master(mask, :);
        if isempty(rows), continue; end

        acqs = unique(rows.Acquisition);

        acq_means = NaN(length(acqs), length(metrics));
        for a = 1:length(acqs)
            acq_rows = rows(rows.Acquisition == acqs(a), :);
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

        for m = 1:length(metrics)
            met       = string(metrics{m});
            surg_mean = mean(acq_means(:, m), 'omitnan');
            n_acq_val = sum(~isnan(acq_means(:, m)));

            new_row = {char(surg), exp_grp, char(grp), char(met), surg_mean, n_acq_val};
            surg_level = [surg_level; cell2table(new_row, 'VariableNames', ...
                {'Surgeon','Expertise','Group','Metric','Mean','N_Acquisitions'})];
        end
    end
end

for s = 1:4
    col = surg_level.Properties.VariableNames{s};
    surg_level.(col) = string(surg_level.(col));
end

fprintf('Surgeon-level rows: %d\n', height(surg_level));

%% -------------------------------------------------------------------------
%  DIAGNOSTIC: Show surgeon-level data summary
%% -------------------------------------------------------------------------
fprintf('\n--- DIAGNOSTIC: Surgeon-level data ---\n');
fprintf('  %-10s %-10s %-25s %-8s %s\n', 'Surgeon','Expertise','Group','Metric','Mean');
fprintf('  %s\n', repmat('-', 1, 65));
for i = 1:height(surg_level)
    fprintf('  %-10s %-10s %-25s %-8s %.4f\n', ...
        char(surg_level.Surgeon(i)), char(surg_level.Expertise(i)), ...
        char(surg_level.Group(i)), char(surg_level.Metric(i)), ...
        surg_level.Mean(i));
end

fprintf('\n  Expert surgeons in surg_level: %s\n', ...
    strjoin(unique(surg_level.Surgeon(surg_level.Expertise == "Expert")), ', '));
fprintf('  Novice surgeons in surg_level: %s\n', ...
    strjoin(unique(surg_level.Surgeon(surg_level.Expertise == "Novice")), ', '));
fprintf('  Groups present: %s\n', ...
    strjoin(unique(surg_level.Group), ', '));

%% -------------------------------------------------------------------------
%  DESCRIPTIVE STATISTICS
%% -------------------------------------------------------------------------

desc_results = table();

for m = 1:length(metrics)
    met = string(metrics{m});
    for c = 1:length(all_conditions)
        cond = string(all_conditions{c});
        for grp = ["Expert", "Novice"]
            mask  = surg_level.Expertise == grp & surg_level.Group == cond & surg_level.Metric == met;
            vals  = surg_level.Mean(mask);
            n_s   = length(vals);
            if n_s == 0, continue; end
            sd_val  = NaN; sem_val = NaN;
            if n_s > 1, sd_val = std(vals); sem_val = sd_val/sqrt(n_s); end

            new_row = {char(met), char(cond), char(grp), ...
                mean(vals,'omitnan'), sd_val, sem_val, ...
                median(vals,'omitnan'), min(vals), max(vals), n_s};
            desc_results = [desc_results; cell2table(new_row, 'VariableNames', ...
                {'Metric','Condition','Expertise', ...
                 'Mean','SD','SEM','Median','Min','Max','n_surgeons'})];
        end
    end
end

%% -------------------------------------------------------------------------
%  HEDGES' G + MANN-WHITNEY U — ALL CONDITIONS
%% -------------------------------------------------------------------------

fprintf('\n--- Computing Hedges'' g and Mann-Whitney U ---\n');

results = table();

for m = 1:length(metrics)
    met = string(metrics{m});
    for c = 1:length(all_conditions)
        cond      = string(all_conditions{c});
        is_shared = ismember(char(cond), shared_conditions);

        mask_exp = surg_level.Expertise == "Expert" & surg_level.Group == cond & surg_level.Metric == met;
        mask_nov = surg_level.Expertise == "Novice" & surg_level.Group == cond & surg_level.Metric == met;

        vals_exp = surg_level.Mean(mask_exp);
        vals_nov = surg_level.Mean(mask_nov);
        n_exp    = length(vals_exp);
        n_nov    = length(vals_nov);
        diff_val = mean(vals_exp,'omitnan') - mean(vals_nov,'omitnan');

        % ── Hedges' g (always computed where possible) ───────────────────
        if n_exp >= 2 && n_nov >= 2
            df   = n_exp + n_nov - 2;
            sd_p = sqrt(((n_exp-1)*std(vals_exp)^2 + (n_nov-1)*std(vals_nov)^2) / df);
            g_note = sprintf('n_exp=%d, n_nov=%d', n_exp, n_nov);
        elseif n_exp >= 2 && n_nov == 1
            df = n_exp - 1; sd_p = std(vals_exp);
            g_note = 'Expert SD proxy (n_nov=1)';
        elseif n_exp == 1 && n_nov >= 2
            df = n_nov - 1; sd_p = std(vals_nov);
            g_note = 'Novice SD proxy (n_exp=1)';
        else
            df = NaN; sd_p = NaN; g_note = 'n=1 both — not calculable';
        end

        if ~isnan(sd_p) && sd_p > 0
            d = diff_val / sd_p;
            J = 1 - 3 / (4*df - 1);
            g = d * J;
        else
            d = NaN; J = NaN; g = NaN;
        end

        % ── Mann-Whitney U (shared conditions, n >= 2 per group) ─────────
        if is_shared && n_exp >= 2 && n_nov >= 2
            [p, ~, stats] = ranksum(vals_exp, vals_nov, 'alpha', alpha);
            U = stats.ranksum - n_exp*(n_exp+1)/2;
            interp_str = "Pending";
        else
            p = NaN; U = NaN;
            if ~is_shared
                interp_str = "g only (unique condition)";
            else
                interp_str = "g only (n<2 in one group)";
            end
        end

        new_row = {char(met), char(cond), ...
            mean(vals_exp,'omitnan'), mean(vals_nov,'omitnan'), diff_val, ...
            U, p, d, g, abs(g), J, interp_str, n_exp, n_nov, g_note};
        results = [results; cell2table(new_row, 'VariableNames', ...
            {'Metric','Condition', ...
             'Mean_Expert','Mean_Novice','Difference', ...
             'U_stat','p_value', ...
             'Cohens_d','Hedges_g','Abs_g','J_correction', ...
             'Interpretation','n_exp','n_nov','g_note'})];
    end
end

% Assign interpretations
for i = 1:height(results)
    if isnan(results.p_value(i)) || ~strcmp(char(results.Interpretation(i)), 'Pending')
        % For g-only rows, assign effect label without significance
        if ~isnan(results.Abs_g(i))
            g = results.Abs_g(i);
            dif = results.Difference(i);
            if g < 0.2,      eff = "Trivial";
            elseif g < 0.5,  eff = "Small";
            elseif g < 0.8,  eff = "Medium";
            else,            eff = "Large";
            end
            if dif > 0,      dir = "Expert > Novice";
            elseif dif < 0,  dir = "Novice > Expert";
            else,            dir = "No difference";
            end
            if strcmp(char(results.Interpretation(i)), 'Pending')
                results.Interpretation(i) = "ns | " + eff + " | " + dir;
            else
                results.Interpretation(i) = results.Interpretation(i) + ...
                    " | " + eff + " | " + dir;
            end
        end
        continue;
    end

    p   = results.p_value(i);
    g   = results.Abs_g(i);
    dif = results.Difference(i);

    if p < 0.001,     sig = "***";
    elseif p < 0.01,  sig = "**";
    elseif p < alpha, sig = "*";
    else,             sig = "ns";
    end
    if g < 0.2,      eff = "Trivial";
    elseif g < 0.5,  eff = "Small";
    elseif g < 0.8,  eff = "Medium";
    else,            eff = "Large";
    end
    if dif > 0,      dir = "Expert > Novice";
    elseif dif < 0,  dir = "Novice > Expert";
    else,            dir = "No difference";
    end
    results.Interpretation(i) = sig + " | " + eff + " | " + dir;
end

%% -------------------------------------------------------------------------
%  SAVE RESULTS
%% -------------------------------------------------------------------------

out_main = fullfile(output_dir, 'analysis3_physiology_expertise_RATKA.csv');
out_desc = fullfile(output_dir, 'analysis3_physiology_descriptive_RATKA.csv');

writetable(results,      out_main);
writetable(desc_results, out_desc);

fprintf('\nSaved: %s\n', out_main);
fprintf('Saved: %s\n', out_desc);

%% -------------------------------------------------------------------------
%  PRINT SUMMARY
%% -------------------------------------------------------------------------

fprintf('\n');
fprintf('========================================================================\n');
fprintf('  ANALYSIS 3 (PHYSIOLOGY): Expert vs Novice — RA-TKA\n');
fprintf('========================================================================\n');
fprintf('  Expert: %s (n=%d) | Novice: %s (n=%d)\n', ...
    strjoin(expert_ids,', '), length(expert_ids), ...
    strjoin(novice_ids,', '), length(novice_ids));
fprintf('\n  *** POWER WARNING ***\n');
fprintf('  n=4 vs n=2 — minimum achievable p ~ 0.133. Hedges'' g is primary outcome.\n');
fprintf('  Benchmarks: Trivial<0.2, Small 0.2-0.5, Medium 0.5-0.8, Large>=0.8\n');
fprintf('------------------------------------------------------------------------\n');

% Descriptive stats
fprintf('\n  DESCRIPTIVE STATISTICS (Mean ± SD)\n');
for m = 1:length(metrics)
    met = string(metrics{m});
    fprintf('\n  Metric: %s\n', char(met));
    fprintf('  %-22s %-10s %12s %12s %12s %6s\n', ...
        'Condition','Group','Mean','SD','Median','n');
    fprintf('  %s\n', repmat('-', 1, 75));
    for c = 1:length(all_conditions)
        cond = string(all_conditions{c});
        for grp = ["Expert","Novice"]
            mask = strcmp(desc_results.Metric,char(met)) & ...
                   strcmp(desc_results.Condition,char(cond)) & ...
                   strcmp(desc_results.Expertise,char(grp));
            if ~any(mask), continue; end
            fprintf('  %-22s %-10s %12.4f %12.4f %12.4f %6d\n', ...
                char(cond), char(grp), ...
                desc_results.Mean(mask), desc_results.SD(mask), ...
                desc_results.Median(mask), desc_results.n_surgeons(mask));
        end
    end
end

% Hedges g + Mann-Whitney
fprintf('\n');
fprintf('========================================================================\n');
fprintf('  HEDGES'' g + MANN-WHITNEY U\n');
fprintf('========================================================================\n');
for m = 1:length(metrics)
    met = string(metrics{m});
    fprintf('\n  Metric: %s\n', char(met));
    fprintf('  %-22s %6s %6s %12s %12s %12s %7s %7s %8s  %-30s\n', ...
        'Condition','n_Exp','n_Nov','Mean_Exp','Mean_Nov','Diff', ...
        'd','g','p_value','Result');
    fprintf('  %s\n', repmat('-', 1, 130));
    for c = 1:length(all_conditions)
        cond = string(all_conditions{c});
        mask = strcmp(results.Metric,char(met)) & strcmp(results.Condition,char(cond));
        if ~any(mask), continue; end
        row = results(mask,:);
        if isnan(row.p_value), p_val = '     N/A '; else, p_val = sprintf('%8.4f',row.p_value); end
        fprintf('  %-22s %6d %6d %12.4f %12.4f %12.4f %7.3f %7.3f %s  %s\n', ...
            char(cond), row.n_exp, row.n_nov, ...
            row.Mean_Expert, row.Mean_Novice, row.Difference, ...
            row.Cohens_d, row.Hedges_g, p_val, char(row.Interpretation));
    end
end

fprintf('\n========================================================================\n');
fprintf('  NOTABLE EFFECTS (|g| >= 0.5)\n');
fprintf('========================================================================\n');
notable = results(results.Abs_g >= 0.5, :);
if ~isempty(notable)
    [~,sidx] = sort(notable.Abs_g,'descend');
    notable  = notable(sidx,:);
    for i = 1:height(notable)
        fprintf('  [%s | %s]  g=%.3f  — %s  (%s)\n', ...
            notable.Metric{i}, notable.Condition{i}, ...
            notable.Hedges_g(i), char(notable.Interpretation(i)), char(notable.g_note(i)));
    end
else
    fprintf('  No effects >= 0.5 found.\n');
end

fprintf('\n========================================================================\n');
fprintf('  INTERPRETATION NOTES\n');
fprintf('========================================================================\n');
fprintf('  Positive g/difference = Expert > Novice.\n');
fprintf('  Higher RMSSD = greater HRV = lower stress.\n');
fprintf('  BR and RVT are correlated — treat convergent findings as one evidence stream.\n');
fprintf('  p-values reported for completeness only — not powered for significance.\n');
fprintf('\nDone.\n');
