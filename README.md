# How much of a sleep-deprivation EEG effect is the pipeline?

Code and results for a specification-curve (multiverse) analysis of resting-state EEG after total sleep deprivation, using the open dataset [OpenNeuro ds004902](https://openneuro.org/datasets/ds004902).

> Luo, V. *How Much of a Sleep-Deprivation EEG Effect Is the Pipeline? A Specification-Curve Analysis of Resting-State EEG in an Open Dataset.* Manuscript in preparation.

The analysis runs every defensible combination of eight analysis decisions (3,888 specifications) on the same recordings. It measures how much the sleep-deprivation effect on theta activity depends on those decisions, and whether the decisions that matter change the estimate or change which participants are analysed.

## Contents

| Path | What it is |
|---|---|
| `R/multiverse_pipeline.R` | Stage 1. Reads the eyes-open recordings, runs all 432 pipelines per recording, and writes the feature table, the 3,888-row specification curve and a variance decomposition. |
| `R/multiverse_inference.R` | Stage 2. Builds the fixed-sample analysis sets, specification curves, variance decompositions, joint permutation tests and the retained-recordings table. |
| `R/sample_composition.R` | Stage 3. Leave-one-participant-out influence analysis of the unrestricted theta-mean multiverse; the six 2-s reference × artifact-handling pipelines run on their own sample and on the 57 common participants; the theta power of the most influential recording. |
| `R/figure1_sample_composition.R` | Figure 1 (three panels) from the Stage 3 tables. |
| `R/figures_manuscript.R` | Figure 2 (specification curves) and Figure 3A (variance shares). |
| `R/headmap.py` | Figure 3B (channel map coloured by the fixed-sample theta effect). Python; needs `mne`, `pandas`, `matplotlib`. |
| `results/` | Every table the paper reports (described below). |
| `figures/` | The manuscript figures as written by the three figure scripts. |

## Data

The EEG recordings are not included. Download ds004902 (Xiang et al., 2024; CC0) from OpenNeuro into a folder named `ds004902` inside this project, or change `data_dir` at the top of `R/multiverse_pipeline.R`.

The dataset has 71 participants, each recorded after normal sleep and after one night of total sleep deprivation, with session order counterbalanced. The analysis uses the eyes-open resting-state recordings (142 recordings). Following the dataset's README, session 1 is coded as normal sleep (NS) and session 2 as sleep deprivation (SD).

## Requirements

R with tidyverse, fs, jsonlite, signal, fastICA, furrr, future and arrow (the analyses in the paper were run in R 4.5.2). patchwork is needed for the combined Figure 1; without it, the panels are saved as separate files. `multiverse_pipeline.R` installs any missing packages except patchwork.

## How to reproduce

1. Open `sleep-eeg-multiverse.Rproj` in RStudio. This sets the working directory to the project folder, which both scripts expect.
2. Run `R/multiverse_pipeline.R`. The full run took 18.6 minutes on four cores (`N_CORES` at the top of the script). Results for each recording are cached in `results/cache/`, so an interrupted run picks up where it stopped. The ICA step starts from random values, so a fresh run can give slightly different numbers for the ICA pipelines than the included `results/features_all_pipelines.parquet`, which is the file the paper used.
3. Run `R/multiverse_inference.R` (under a minute). It writes the analysis tables to `results/`. The permutation tests use a fixed seed (20260927), so a rerun reproduces the reported p-values exactly.
4. Run `R/sample_composition.R` (about a minute) for the influence analysis, then `R/figure1_sample_composition.R`, `R/figures_manuscript.R` and `python R/headmap.py` for the figures. All scripts run from the repository root.

To check the statistics without the raw data, skip step 2: `results/features_all_pipelines.parquet` is included, and `multiverse_inference.R` runs from it alone.

## Specification space

| Decision | Options |
|---|---|
| Reference | average; linked mastoids (TP9/TP10); Cz |
| Artifact handling | amplitude threshold only; threshold plus ocular ICA (fastICA, 20 components; components correlating \|r\| > 0.5 with the Fp1/Fp2 mean removed) |
| Epoch length | 2, 4, 8 or 10 s, non-overlapping |
| Normalisation | absolute; relative to 1–30 Hz; relative to 4–30 Hz |
| Aperiodic component | retained; removed (log–log linear fit over 2–30 Hz, residual floored at zero) |
| Outcome statistic | mean, standard deviation or coefficient of variation of band power across epochs |
| Band | theta 4–8 Hz; alpha 8–13 Hz; beta 13–30 Hz |
| Region | frontal (16 channels); centro-temporal (28); parieto-occipital (17), following Cui et al. (2026) |

The first six decisions give 432 pipelines per recording; with three bands and three regions, 3,888 specifications.

Steps shared by every pipeline:

- 0.5–45 Hz zero-phase Butterworth band-pass
- bad channels removed where the z-score of log variance, computed from the median and the median absolute deviation, exceeds 3
- Hann-windowed FFT power spectra over 1–45 Hz
- epochs with any sample beyond ±100 µV rejected
- a recording dropped from a pipeline if fewer than 10 clean epochs remain

## Analysis sets

| Set | File | Specifications | Participants |
|---|---|---|---|
| Unrestricted | `specification_curve.csv` | 3,888 | 14–71, depending on the pipeline |
| Fixed sample, primary ("core B") | `coreB_epoch2_specification_curve.csv` | 972 (2-s epochs) | 57, complete under every reference × artifact-handling combination |
| Fixed sample, sensitivity ("core A") | `coreA_epoch2and4_specification_curve.csv` | 1,944 (2- and 4-s epochs) | 40 |

## Results files

| File | Rows | Contents |
|---|---|---|
| `features_all_pipelines.parquet` | 424,278 | One value per recording × pipeline × band × region. Columns: the eight decisions (`reference`, `artifact`, `epoch_sec`, `normalise`, `aperiodic`, `outcome`, `band`, `roi`), `subject`, `condition` (NS/SD), `value`, `reject_pct` (% of epochs rejected). |
| `specification_curve.csv` | 3,888 | Unrestricted multiverse, one row per specification. |
| `coreB_epoch2_specification_curve.csv` | 972 | Primary fixed-sample set. |
| `coreA_epoch2and4_specification_curve.csv` | 1,944 | Sensitivity set. |
| `variance_decomposition.csv` | 9 | Main-effects ANOVA of the effect size across all 3,888 specifications. |
| `variance_decomposition_all_sets.csv` | 57 | Share of variance by decision for each set (`set`: `unrestricted_theta_mean`, `coreB_all`, `coreB_all_2way`, `coreB_theta_mean`, `coreA_theta_mean`). |
| `permutation_joint_tests_coreB.csv` | 27 | Joint permutation tests on core B: 9 band × outcome families × 3 statistics. |
| `sample_retained_table.csv` | 6 | Recordings retained (of 142) by reference, artifact handling and epoch length (Table II in the paper). |
| `leave_one_out_influence.csv` | 71 | The unrestricted theta-mean multiverse recomputed with each participant removed in turn: median effect, correlation between sample size and effect, median effect by artifact handling, and each decision's share of variance. |
| `six_pipelines_own_vs_common.csv` | 6 | Each 2-s reference × artifact-handling pipeline on its own sample and on the 57 participants common to all six (Figure 1C). |
| `influential_recording.csv` | 4 | Centro-temporal absolute theta power of the most influential participant (47) in both sessions under both artifact-handling options, with its SD across epochs and, for the other 70 participants under the same pipeline, the median, MAD, 5th and 95th percentiles and the resulting modified z-score (Figure 1B). |

**Specification curves.** The eight decisions, then `n` (participants), `dz` (Cohen's d_z of the SD − NS difference), `p` (paired t-test), `p_fdr` (Benjamini–Hochberg within the set) and `sig` (`p_fdr` < 0.05).

**Permutation tests.** For each family, each participant's SD − NS difference was sign-flipped at random, with the same flip across all specifications for that participant, and the whole curve was recomputed. This was repeated 2,000 times.

- `statistic`: `median_dz`, `share_sig_expected` (share of specifications with p < 0.05 in the expected direction) or `mean_z`
- `null_median`, `null_q95`: median and 95th percentile of the 2,000 permuted values
- `p_perm`: (1 + permutations ≥ observed) / 2,001, one-sided in the direction given by `expected_sign`

## Data source

Xiang, C., Fan, X., Bai, D., Lv, K., & Lei, X. (2024). A resting-state EEG dataset for sleep deprivation. *Scientific Data*, 11, 427. https://doi.org/10.1038/s41597-024-03268-2

## License

Code is released under the MIT License (see `LICENSE`). The EEG data belong to the dataset's authors and are distributed by OpenNeuro under CC0.
