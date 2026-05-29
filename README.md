# fNIRS Analysis Pipeline

## Overview

This repository contains the MATLAB scripts used in the analysis for the thesis titled *'The Moderating Role of Surgical Expertise on Cortical and Autonomic Demands of Robotic and Conventional Hip and Knee Arthroplasty'*. The pipeline processes fNIRS and physiological data from surgeons performing hip and knee arthroplasty under robotic-assisted and conventional techniques, running GLM-based neural activation analyses and comparing outcomes across groups.

The pipeline proceeds in four broad stages:

1. **GLM fitting** – model neural responses at the individual acquisition level
2. **Beta extraction** – extract canonical HRF beta coefficients from fitted models
3. **Synthesis** – average betas across surgeries and group channels into brain regions
4. **Statistical analysis** – between group comparisons with effect size estimation

---

## Dependencies

- SPM12 with the `spm_fnirs` toolbox
- Homer3
- MATLAB R2021a or later (uses `readtable`, `writetable`, `ranksum`, etc.)
- `acquisitions.mat` – a structured array containing pre-processed fNIRS data, physiological signals, event onsets, and metadata for all acquisitions

---

## Data Structure

All scripts expect an `acquisitions.mat` file at `/Volumes/Dylan SSD/acquisitions.mat`. Each entry in the `acquisitions` struct has the following fields:

| Field | Description |
|---|---|
| `name` | Acquisition ID in the format `surgeonID_joint_tech_date_timebin` |
| `data` | fNIRS signal array: `samples × channels × 3` (HbO, HbR, HbT) |
| `nirstime` | Time vector (seconds) for NIRS samples |
| `nirs_abs_start` | Absolute datetime of NIRS recording start |
| `phys_abs_start` | Absolute datetime of physiology recording start |
| `stage` | Cell array of surgical event names |
| `stagetime` | Onset times for each event (seconds from NIRS start) |
| `heartrate` | Heart rate resampled to NIRS rate |
| `breathingrate` | Breathing rate resampled to NIRS rate |
| `respvoltime` | Respiratory volume-time (RVT) resampled to NIRS rate |
| `badchannels` | Indices of channels to exclude from GLM |

Acquisition names encode key metadata:
- **Joint**: `H` = hip, `K` = knee
- **Tech**: `C` = conventional, `R` = robotic-assisted

---

## Stage 1 – GLM Fitting

### `standard_glm_analysis_acq.m`

Runs the per-event SPM fNIRS GLM for a single acquisition.

- Loads fNIRS data and event onsets from the acquisitions struct
- Optionally downsamples NIRS data to 1 Hz
- Builds physiological noise regressors (HR, BR, RVT)
- Fits the GLM using `spmfnirsflow`

**Key parameters to set:**

| Parameter | Default | Description |
|---|---|---|
| `acq_idx` | `2` | Index into acquisitions structure |
| `signal` | `'HbDiff'` | `'HbDiff'`, `'HbO'`, or `'HbR'` |
| `d_sample` | `1` | `1` = downsample to 1 Hz |
| `add_regs` | `1` | `1` = include ECG/Resp regressors |

**Output:** `GLM TDD Phys/<acq_name>/` containing `NIRS.mat`, `onsets.mat`, and the fitted SPM results.

### `standard_glm_grouped.m`

Variant of the GLM that groups surgical events into higher-level conditions (`shared_same`, `shared_diff`, `conventional_only`, `robotic_only`) before fitting.

**Output:** `GLM TDD Phys Grouped/<acq_name>/`

---

## Stage 2 – Beta Extraction

### `b_value_command.m`

Extracts GLM beta values (canonical HRF regressor, `bf(1)`) from a per-event SPM result.

- Selects the `bf(1)` regressor for each event (every 3rd regressor, excluding nuisance terms)
- Filters to spatially registered (valid) channels using the `POS.mat` digitisation file
- Saves wide-format and a summary text file

**Key parameters to set:**

| Parameter | Default | Description |
|---|---|---|
| `acq_idx` | `2` | Acquisition index |
| `n_nuisance` | `4` | Nuisance regressors: HR, BR, RVT, constant (set to `1` if `add_regs = 0`) |

**Outputs** (saved to the acquisition's results folder):
- `betas_HbDiff_bf1.csv` – all channels
- `betas_HbDiff_bf1_validch.csv` – valid channels only
- `beta_summary.txt`

### `extract_betas_grouped.m`

Equivalent to `b_value_command.m` but for the grouped GLM. Extracts one beta per grouped condition rather than one per surgical event.

**Outputs:**
- `betas_grouped_HbDiff_bf1.csv`
- `betas_grouped_HbDiff_bf1_validch.csv`
- `betas_grouped_HbDiff_bf1_long.csv` – long-format with `Acquisition`, `Surgeon`, `Joint`, `Tech`, `Group`, `Channel`, `Beta` columns, for easy import into R or Python

---

## Stage 3 – Synthesis Across Surgeries

Aggregates per-event beta CSVs across multiple surgeries performed by the same surgeon under the same joint/technique combination.

**Steps:**
1. Identifies all surgeon × joint × tech groups
2. Finds common events and channels across surgeries within each group
3. Averages betas across surgeries (zeros treated as missing/excluded)
4. Maps channels to 8 predefined brain regions
5. Saves per-group and master output files

**Outputs** (saved to `/Results/Synthesised/`):
- `betas_avg_HbDiff_<key>.csv` – channel-level averaged betas per group
- `betas_region_HbDiff_<key>.csv` – region-level averaged betas per group
- `master_channel_HbDiff.csv` – long-format master table, all groups, channel level
- `master_region_HbDiff.csv` – long-format master table, all groups, region level

### `synthesise_betas_grouped.m`

Equivalent to `synthesis_betas.m` for the grouped GLM. Handles the union of conditions across surgeries (since not all grouped conditions appear in every surgery).

**Outputs** (saved to `/Results/Synthesised_Grouped/`):
- `grouped_betas_avg_HbDiff_<key>.csv`
- `grouped_betas_region_HbDiff_<key>.csv`
- `master_grouped_channel_HbDiff.csv`
- `master_grouped_region_HbDiff.csv`

### `Synthesis_betas_grouped_SD.m`

Variant of `synthesise_betas_grouped.m` that additionally computes the spatial standard deviation of betas within each brain region (used as the SD input for Hedges' g calculation in Analysis 2).

### `Synthesise_physiology_grouped.m`

Aggregates per-group physiology summaries (from `physiology_by_group.csv`) across surgeries.

**Measures synthesised:** HR mean/SD/min/max, BR mean/SD/min/max, RMSSD, cross-surgery SD for HR, BR, and RMSSD.

**Outputs** (saved to `/Results/Synthesised_Physiology_Grouped/`):
- `phys_grouped_<key>.csv` – per surgeon × joint × tech
- `master_physiology_grouped.csv`

---

## Physiology Extraction

### `physiology_by_event.m`

Extracts HR, BR, RVT, and approximate RMSSD for each individual surgical event epoch (onset to next onset). Works directly from the acquisitions struct — no GLM required.

**Output:** `Physiology Events/<acq_name>/physiology_by_event.csv`

**Columns:** `Event`, `Onset_s`, `Duration_s`, `N_samples`, `HR_mean`, `HR_sd`, `HR_min`, `HR_max`, `BR_mean`, `BR_sd`, `BR_min`, `BR_max`, `RVT_mean`, `RVT_sd`, `RVT_min`, `RVT_max`, `RMSSD`

### `physiology_by_group.m`

Equivalent to `physiology_by_event.m` but aggregated by grouped condition (`shared_same`, `shared_diff`, etc.).

---

## Stage 4 – Statistical Analyses

### `analysis1_technique_THA.m`

**Analysis 1:** Compares RA-THA vs C-THA neural activation across brain regions for shared surgical conditions.

- **Test:** Mann-Whitney U (non-parametric, appropriate for small n)
- **Effect size:** Hedges' g (bias-corrected Cohen's d)
- **Unit of analysis:** Surgeon (betas averaged across procedures per surgeon)
- **Scope:** THA only (`Joint == "H"`)

**Conditions tested:** `shared_same`, `shared_diff` (inferential); `robotic_only`, `conventional_only` (descriptive only)

**Outputs** (saved to `/Results/Analysis/`):
- `analysis1_shared_conditions_HbDiff.csv`
- `analysis1_unique_conditions_HbDiff.csv`

### `analysis1_physiology_THA.m`

Equivalent Analysis 1 for physiological measures (HR, BR, RMSSD) rather than fNIRS betas.

### `analysis2_technique_expertise_THA.m`

**Analysis 2:** Computes Hedges' g for robotic vs conventional comparison within each surgeon × joint group, using spatial SD across channels as the variability estimate.

**Output:** `hedges_g_grouped_HbDiff.csv`

### `analysis2_physiology_THA.m`

Physiology equivalent of `analysis2_technique_expertise_THA.m`.

### `analysis3_expertise_RATKA.m`

**Analysis 3:** Expertise effects within the RA-TKA group.

### `analysis3_physiology_RATKA.m`

Physiology equivalent of Analysis 3.

---

## Brain Region Definitions

Channels are mapped to 8 frontal and parietal regions based on fNIRS optode digitisation and Brodmann area registration.

| # | Name | Brodmann Areas | Channels |
|---|---|---|---|
| 1 | Dorsolateral PFC | BA9, BA46 | 1–5, 7, 20, 22–24, 26 |
| 2 | Orbitofrontal Cortex | BA11, BA47 | 2, 11, 17, 25 |
| 3 | Medial PFC | BA9, BA10 | 8, 9, 13–15, 18, 19 |
| 4 | Superior Frontal Gyrus | BA8, BA9 | 6, 10, 12, 16, 21 |
| 5 | Primary Somatosensory Cortex | BA1, BA2, BA3 | 33–35, 40–42 |
| 6 | Superior Parietal Lobule | BA5, BA7 | 27, 31, 34, 37, 39, 41, 43, 45 |
| 7 | Angular Gyrus | BA39 | 28–30, 32, 36, 38, 44, 46–48 |
| 8 | Supramarginal Gyrus | BA40 | 29, 30, 33, 40, 47, 48 |

> Note: Channels can belong to more than one region.

---

## Pipeline Order

```
acquisitions.mat
│
├─ standard_glm_analysis_acq.m       ← per-event GLM
│     └─ b_value_command.m           ← extract per-event betas
│           └─ synthesise_betas.m    ← average + regionalise
│
├─ standard_glm_grouped.m            ← grouped GLM
│     └─ extract_betas_grouped.m     ← extract grouped betas
│           ├─ synthesise_betas_grouped.m      ← average + regionalise
│           └─ synthesise_betas_grouped_SD.m   ← with spatial SD
│
├─ physiology_by_event.m             ← per-event physiology
│     └─ physiology_by_group.m       ← grouped physiology
│           └─ synthesise_physiology_grouped.m
│
└─ Analysis
      ├─ analysis1_technique_THA.m           ← RA-THA vs C-THA (fNIRS)
      ├─ analysis1_physiology_THA.m          ← RA-THA vs C-THA (phys)
      ├─ analysis2_technique_expertise_THA.m ← within-surgeon Hedges' g
      ├─ analysis2_physiology_THA.m
      ├─ analysis3_expertise_RATKA.m
      └─ analysis3_physiology_RATKA.m
```

---

## Output Directory Structure

```
/Volumes/Dylan SSD/DYLAN/
├── Data/
│   ├── hip/
│   │   ├── c/<surgeonID>/
│   │   │   ├── GLM TDD Phys/<acq_name>/         ← per-event GLM
│   │   │   └── GLM TDD Phys Grouped/<acq_name>/ ← grouped GLM
│   │   └── r/<surgeonID>/
│   └── knee/
│       ├── c/<surgeonID>/
│       └── r/<surgeonID>/
└── Results/
    ├── Synthesised/                              ← per-event synthesis outputs
    ├── Synthesised_Grouped/                      ← grouped synthesis outputs
    ├── Synthesised_Physiology_Grouped/
    └── Analysis/                                 ← statistical analysis outputs
```
