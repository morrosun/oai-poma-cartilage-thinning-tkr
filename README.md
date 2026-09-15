# Cartilage thinning rate and knee replacement in the OAI POMA nested case-control — analysis code

> **EN** — Analysis code for a methodological study of how quantitative cartilage relates to
> subsequent knee replacement (KR), using the Osteoarthritis Initiative (OAI) **POMA
> (Pivotal OAI MRI Analyses)** TKR nested case-control subset. The three questions addressed
> are: (i) how to fit a joint model when both members of a matched pair share one index visit;
> (ii) whether empirical-Bayes shrinkage of the exposure fabricated the observed association;
> (iii) whether a latent-class trajectory solution is a subtype or the tail of a continuum.
>
> **中文** — 本研究「定量软骨与后续膝关节置换（KR）关系」的方法学分析代码，使用骨关节炎
> 计划（OAI）**POMA（Pivotal OAI MRI Analyses）** TKR 巢式病例对照子集。回答三个问题：
> ①配对双方共享同一索引访视时，联合模型应如何拟合；②暴露的 EB 收缩是否人为制造了关联；
> ③潜类轨迹解是亚型，还是连续分布的上尾。

**Code only.** No result files and no patient-level records are distributed here — derived tables
regenerate from the source data under its own data-use agreement (see *Data availability*).

---

## Study design at a glance / 研究设计

| Element | Detail |
|---|---|
| Design | Nested case-control, prospectively matched on the index visit |
| Released subset | 225 knee pairs / 450 knee units |
| Primary analysis set | 191 complete pairs / 382 knees / 1248 observations |
| Follow-up visits | baseline and months 12, 24, 36, 48 (all pre-operative) |
| Exposure | Empirical-Bayes (EBLUP) thinning rate of central medial femorotibial cartilage (cMFTC), mm/yr |
| Shrinkage control | Unshrunk per-knee least-squares slope |
| Outcome | Knee replacement during follow-up |
| Two design facts | The index visit is identical within a pair (225/225), and every analysed measurement precedes it |

The identical index visit is what makes a within-pair conditional likelihood attractive: the
baseline hazard and every within-pair-constant term cancel, so only within-pair contrasts are
identified. It is also what makes a stratified joint model with pair-constant covariates
rank-deficient — see `simulation/12b_sim_jmbayes2_probe.R`.

---

## Repository contents / 仓库内容

```
data_prep/     # dataset acquisition, integrity checks, cross-dataset linkage
poma/          # OAI POMA nested case-control — the main analysis of the study
simulation/    # known-truth simulation studies calibrated to the real POMA design
cth_maps/      # companion trajectory / continuum analysis of the CTh-Maps release
```

### `data_prep/`

| Script | Role |
|---|---|
| `download_oaizibcm_hfmirror.py` | Fetch the OAIZIB-CM imaging release from its public mirror |
| `zenodo_get_resume.py` | Fetch the CTh-Maps release from Zenodo, resumable |
| `unzip_oaidatasets.py` | Unpack both archives and verify the file counts |
| `data_capability.py` | Inventory what each release actually contains: trajectory structure, ID overlap, mask label values |
| `linkage_check.py` | Establish that OAIZIB-CM and CTh-Maps join on the OAI subject ID, and resolve the empirically coded knee side |

`data_capability.py` and `linkage_check.py` write their transcripts to `data_prep/out/`.

### `poma/` — main analysis (run in this order)

| Script | Role |
|---|---|
| `02_long_reshape_and_pilot.py` | Reshape the POMA release to one row per knee and visit (unit = knee); pilot summaries and the first slope estimates. **Supersedes `01_long_reshape_and_pilot.py`, whose conditional-logistic Hessian diverged and which admitted single-member strata.** |
| `03_diag_mixedlm.py` | Convergence and variance-component diagnostics for the per-knee slope model |
| `04_shrinkage_moments.py` | Moment cross-check of the EB shrinkage, decomposing Var(OLS slope) into true-slope variance plus estimation variance |
| `05_merge_baseline_covariates.py` | Merge age, sex, BMI, WOMAC pain, PASE and KL grade onto the knee-level table (≥98% coverage required) |
| `06_clogit_adjusted.R` | Within-pair conditional logistic regression over the M0–M7 adjustment ladder |
| `07_probe_jmbayes2.R` | Records the off-the-shelf failure mode (`JMbayes2` → `chol(): decomposition failed`) |
| `08_joint_model.R` | Stratified `JMbayes2` fit, retained as the comparator |
| `jm_core.R` | Shared sampler: longitudinal submodel, within-pair conditional likelihood |
| `09_joint_conditional.R` | The joint model under the within-pair conditional likelihood |
| `10_jm_multimetric.R` | Compartment-by-compartment joint models for the six cartilage metrics |
| `10b_gradient_test.R` | The three formal medial-to-lateral contrasts, computed on one posterior draw set |
| `10c_clogit_seq_primary.R` | Progressive covariate adjustment for the pre-specified primary metric |
| `10d_verify_lme_primary.R` | Verifies the longitudinal submodel against an independent `nlme::lme` fit |
| `10e_minvisits_primary.R` | Refits both models requiring both knees to have ≥J visits (J = 2, 3, 4) |
| `11_trajectory_lcmm.R`, `31_cfg_poma_lcmm.R` | Latent-class mixed models for K = 1–6, best of eight random starts |
| `28_lcmm_continuum_check.R` | The three tests that distinguish a subtype from a continuum tail |
| `30_class_continuum_check.R` | The same three tests, refactored into a config-driven, dataset-agnostic engine (the transferable form) |
| `33_verify_poma_equivalence.py` | Verifies that the generic engine reproduces the POMA-specific result exactly |
| `22_cth_score_corroboration.R` | Cross-pipeline corroboration against an independent automatic severity score |
| `12_manuscript_tables.py` | Derives the numeric bodies of Tables 1–4 from the fitted objects |

### `simulation/` — known-truth simulation studies

All scenarios are calibrated to the real primary analysis set: the visit schedule of every
simulated pair is **resampled from the 191 observed pairs**, so the visit counts and the
case-versus-control visit-density asymmetry are those of the study rather than an assumption.

| Script | Role |
|---|---|
| `sim_core.R` | Shared engine: data generation, per-knee slope estimation, conditional-likelihood fitting |
| `10_calibrate_poma.py` | Derives the visit schedule and generating parameters from the observed primary set |
| `11_sim_calibrate_aS.R` | Solves numerically for the true exposure coefficient that reproduces the observed EBLUP estimate |
| `12_sim_degeneracy.R` | Identifiability of the stratum baseline hazard, risk-set composition, and the estimator comparison |
| `12b_sim_jmbayes2_probe.R` | Controlled grid isolating the model specification under which the stratified joint model fails |
| `13_sim_shrinkage.R` | Null, primary and truncated-visit scenarios comparing EBLUP with unshrunk slopes |
| `14_sim_continuum.R` | Applies the primary latent-class pipeline to data generated from a single continuous distribution |

### `cth_maps/` — companion trajectory / continuum analysis

These scripts analyse the CTh-Maps release (45,343 knee-timepoints) rather than the POMA pairs,
and support the companion analysis of whether cartilage-score trajectories form subtypes.
The cross-pipeline corroboration reported in the main study uses `poma/22_cth_score_corroboration.R`.

| Script | Role |
|---|---|
| `01_cth_trajectory_lcmm.R` | Latent-class mixed models for the score trajectory, k-means seeded |
| `02_cth_blup_gmm.py` | Two-step alternative: per-knee BLUPs followed by Gaussian mixtures |
| `03_cth_slope_continuum.py` | Modality-bandwidth scan, Silverman critical-bandwidth test, mixture-component geometry |
| `04_cth_lcmm_k2_rescue.R` | K = 2 rescue fits from alternative starts |
| `05_cth_lcmm_variance_partition.R` | How much of the level versus the rate variance each solution explains |
| `06_cth_lcmm_chain.R` | Forward-selection chain across K |
| `07_cth_lcmm_figures.R` | Trajectory, class-rate and model-selection figures |

---

## Environment / 环境

- **R ≥ 4.6.0** with `survival`, `nlme`, `lcmm` (2.2.2), `JMbayes2`, `splines`, `MASS`, `mvtnorm`, `ggplot2`, `dplyr`, `tidyr`
- **Python ≥ 3.10** with `pandas`, `numpy`, `scikit-learn`, `python-docx` (table export only)
- **PostgreSQL / MATLAB are not required** — no database or MATLAB component is used in this analysis

---

## Reproducibility notes / 复现说明

- **Paths.** Every script carries the absolute paths of the original working tree
  (`D:/BaiduSyncdisk/OAI/...`) in its `BASE` / `ROOT` / `OUT` constant near the top. Adjust that
  constant to your clone location before running; nothing else is machine-specific.
- **Intermediate files are not shipped.** The scripts consume and produce knee-level derived tables
  (for example `poma_analysis_long.csv`) that are **deliberately excluded**, because they are
  re-identifiable OAI records. Run the pipeline from the start to regenerate them.
- **Run order matters.** `simulation/10_calibrate_poma.py` depends on the analysis set produced by
  `poma/02_long_reshape_and_pilot.py` plus `poma/05_merge_baseline_covariates.py`;
  `simulation/11_sim_calibrate_aS.R` must run before `12`–`14`.
- **Long-running steps.** `09_joint_conditional.R`, `10_jm_multimetric.R` and `11_trajectory_lcmm.R`
  are the computationally heavy ones; the simulation scripts each take minutes to hours depending
  on the number of replicates set at the top of the file.
- **Random seeds** are set inside each simulation script; the reported figures come from the
  replicate counts recorded in the scripts as shipped.

---

## Data availability / 数据可用性

Source imaging and clinical data are **not** redistributed here.

- The **OAI** (including the POMA TKR nested case-control subset and the chondrometrics
  measurements) is available to researchers through the NDA at <https://nda.nih.gov/oai/> under a
  data-use agreement.
- The **CTh-Maps** cartilage-severity release is distributed under its own Zenodo record
  (<https://doi.org/10.5281/zenodo.18745638>); note its licence does not permit redistribution of
  derivatives, so only the code that consumes it is published here.

This repository contains **code only**. To reproduce the reported numbers, obtain the source data
under the agreements above and run the pipeline in the order given.

---

## Citation / 引用

The analysis code is archived at Zenodo:

| | |
|---|---|
| **Version-specific DOI — cite this one** | <https://doi.org/10.5281/zenodo.22771313> (`v1.0.0`) |
| Concept DOI — always resolves to the latest version | <https://doi.org/10.5281/zenodo.22771312> |

```
Wang K. Cartilage thinning rate and knee replacement in the OAI POMA nested case-control:
analysis code (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.22771313
```

---

## License / 许可

Code released under the **MIT License** (see `LICENSE`). No licence is asserted over OAI or
CTh-Maps data, which remain subject to their own terms.
