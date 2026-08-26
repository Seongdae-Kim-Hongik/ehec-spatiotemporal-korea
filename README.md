# EHEC Spatiotemporal Analysis — South Korea, 2001-2024

Reproducible analysis code for:

> **Changing Epidemiology of Enterohemorrhagic *Escherichia coli* Infections in South Korea, 2001-2024: Age-Standardized Bayesian Spatiotemporal Analysis.** Kim S, Chun BC. Submitted to *JMIR Public Health and Surveillance* (2026).

> **Note on version 1.1 (August 2026).** During pre-submission checks we found that several district covariates had been
> ingested with unmeasured years recorded as zero rather than missing, which produced spurious associations. All covariates
> were recovered from source, every affected model was refitted, and a preregistered multiverse was added. Two associations
> reported in earlier drafts (bars/pubs, fiscal independence) do not survive; groundwater use is the only covariate whose
> credible interval excludes 1 in the main model. The incident, its root causes, and the resulting data-quality gates are
> documented in [`DATA_QUALITY_INCIDENT.md`](DATA_QUALITY_INCIDENT.md); the analysis plan fixed before estimates were
> inspected is in [`PREREGISTRATION.md`](PREREGISTRATION.md).


> **Note on version 2.0 (August 2026).** The covariate stage was rebuilt from a full census of every district-level indicator available in the source files (4,305 columns), with 108 indicators prespecified in 7 transmission-pathway domains and 85 entering a 3,000-specification multiverse (scripts 520-534). The final 20-covariate model, period analyses, and figures were re-estimated (540-544). The beef-cattle and dairy exposure was re-sourced from the quinquennial Census of Agriculture, Forestry and Fisheries (2010/2015/2020 district farm counts; scripts 545-548) so that measured values precede the study window; the quarterly livestock monitoring series is retained as an independent-source sensitivity analysis. The per-indicator adjudication file is included as `COVARIATE_PREREG_POOL.csv`.

## Data

The analyses use aggregated official statistics only (no individual-level or identifiable data):

| Source | Content | Access |
|---|---|---|
| KDCA statutory notifiable-disease system | Annual EHEC notifications by district (si-gun-gu), 10-year age band, and sex, 2001-2024 | https://dportal.kdca.go.kr |
| MFDS foodborne-illness surveillance | Month/week-resolved EHEC series (seasonality analyses) | https://www.foodsafetykorea.go.kr |
| KOSIS (Statistics Korea) | District age-specific population; district covariate panel (2011-2024) | https://kosis.kr |

Input files are expected as CSV extracts of these public statistics; column schemas are documented at the top of each script. District boundary shapefile: harmonized 229-district (si-gun-gu) map, EPSG:5179.

## Running the code

Scripts locate inputs through two environment variables, so no absolute paths are hard-coded:

```bash
export EHEC_ROOT=/path/to/analysis   # analysis root: holds data_cache/ and results_ST/
export EHEC_DATA=/path/to/sources    # source extracts: FBD_DATA_ZIP/ (boundary shapefile) and the CSV extracts
```

Both default to the current directory (`.` and `./data`) if unset. Required R packages: `INLA`, `sf`, `spdep`, `dplyr`, `tidyr`, `stringr`, `e1071`, `car`. Required Python packages: `pandas`, `numpy`, `matplotlib`, `python-docx`.

## Pipeline (run in order)

### Step 0 — Data assembly
- `450_agestd_assemble_generic.py` — builds the analysis panels: district x year x age-band cases, indirectly age-standardized expected counts (E_std), national crude/age-standardized trends.

### Step 1 — Descriptive and temporal analyses (2001-2024) → Figures 1, 3; eFigures 2-15
- `380_sexage_descriptive.py`, `388_kdca_enrich2.py`, `389_kdca_enrich3.py` — sex- and age-specific descriptive epidemiology (eFigures 2-4).
- `370_stl_national.R` — STL decompositions, crude and age-standardized (eFigures 5-6).
- `462_person_rate_season.py` — age-specific incidence rates and seasonal concentration (eFigure 9).
- `483_age_rate_ratio.py` — child-to-older rate ratios + NB GLM age-by-year interaction test (Figure 1C, eFigure 12).
- `490_child_share_test.py` — demography-only null decomposition of the falling child share (eFigure 15).
- `486_summer_share.py` — July-August share; Gini vs uniform and constant-seasonality null floors (Figure 3A, eFigure 13).
- `383_wavelet_split.py` — Morlet CWT, early/late halves (eFigure 8).
- `487_wavelet_changetest.py`, `489_wavelet_multimethod.py`, `489b_amp_null.py` — change-point, permutation, AR(1)-surrogate, count-matched null, harmonic-amplitude tests (eTable 24, eFigure 14).
- `492_null_envelope.py` — yearly annual-band share vs count-matched constant-seasonality envelope (Figure 3B).

### Step 2 — Spatial and spatiotemporal models (2011-2024) → Figure 2; Tables 1-2; eTables 3-23
- `451_agestd_model_generic.R` — generic age-standardized NB BYM2 spatiotemporal engine (R-INLA).
- `418_ehec_table1_house.R` — main 16-covariate model (Table 2), posterior spatial RR, exceedance.
- `412_agestd_maps.R`, `412b_jama_maps.R`, `415_agestd_exceedance.R` — maps (Figure 2A, eFigure 1).
- Emerging hot spot analysis (Getis-Ord Gi* + Mann-Kendall on EB-smoothed age-standardized ratios) — Figure 2B, eTable 19.
- `484_spatial_4yr_split.R`, `485_period_interaction.R` — calendar-window splits + formal covariate-by-period interaction (eTables 22-23, eFigure 16).
- `498_covid_split.R`, `499_covid_map.R` — pre- vs post-COVID-19 split (2011-2019 vs 2020-2024) + interaction test and maps (eFigure 17, eTable 25).
- Robustness: `430_agestd_supp_robust.R`, `322_8graph.R` (eight neighborhood graphs), `323_M1M6_ladder.R` (model ladder), `290_RSR_sensitivity.R` (restricted spatial regression), `292_ZINB_sparse.R` / `296_ZINB_refit_final.R` (zero-inflation), `293_PIT_histograms.R` / `297_randomized_PIT.R` (calibration), `463_fullscreen_transform.py` (2,329-indicator negative screen).

### Step 3 — Publication figures
- `493_jama_fig1_person.py`, `494_jama_fig2_place.py`, `496_jama_fig2_time.py` — assemble Figures 1-3.

### Step 4 — Covariate recovery, preregistered multiverse, and behavioral-domain tests (version 1.1)
- `503_full_covariate_audit.py` — quality audit of all 46 covariates: per-year share of zeros (not just missingness), value ranges, unit table, and confounder correlations.
- `500_recover_covariates.py`, `504_recover_handwashing.py`, `505_recover_remaining.py` — recovery of covariates whose source columns recorded uncollected years as zero, including density-versus-count corrections and coalescing of indicators split across differently spelled columns.
- `515_recover_rds_cache.R` — rebuilds the cached covariate matrix used by the supplementary model scripts, which the CSV-level recovery did not reach.
- `506_multiverse.R`, `507_specification_curve.py` — preregistered multiverse over 400 specifications (domain-representative covariate sets x ruminant-exposure form x transformation x spatial structure x offset x adjacency graph), with decision rules fixed before any estimate was inspected; specification curve (eFigure 16, eTable 26).
- `508_detection_adjusted.R` — reestimation with case-ascertainment proxies entered jointly (eTable 27).
- `510_domain6_health_behaviour.py`, `512_d6_urbanicity_test.R` — health-behavior and preventive-care domain (influenza vaccination, infant check-ups, nutrition-label literacy, walking, smoking, drinking), tested singly and jointly (eTable 28).
- `511_kchs_common_factor.py` — principal-component diagnosis of survey-derived indicators, which share a dominant rural-urban component.
- `513_soap_forensics.py`, `514_soap_signflip.R` — sequential-adjustment and confounding checks for the handwashing-with-soap association (eTable 29).


### Step 5 — Full covariate census, preregistered pool, and final model (version 2.0)
- `520_covariate_triage.py`, `521_covariate_manual_adjudication.py` — census of all 4,305 district-level columns; rule-based triage plus per-indicator manual adjudication against the transmission framework (output: `COVARIATE_PREREG_POOL.csv`).
- `530_build_prereg_panel.py` — assembles the 108-indicator district-year panel with the data-quality gate (per-year zero-share audit, denominator integrity, renormalization prohibition); `531_z6_vif_design.py` — survey-indicator PC diagnosis and within-domain VIF pruning (85 modeled).
- `532_prereg_stage1_domain.R`, `533_prereg_multiverse.R`, `534_prereg_mv_summary.py` — 3,000-specification stratified multiverse (12 parallel workers) and its summary (Tables S13-S15).
- `540_final_main_model.R` — final 20-covariate BYM2+RW1+type I model (Table 2, Table S11, maps); `542_final_period_analyses.R` — 4-year windows, COVID-19 split, and period-interaction tests (Tables S8-S9, S17); `543_final_spec_curve.py`, `544_det_windows_final.py` — Figures S16-S17.
- `545_census_livestock.py` — parses the 2010/2015/2020 agricultural-census district livestock files, maps them to the 229-district frame (with administrative-boundary aliases), and interpolates 2011-2024; a national-total assertion guards against double counting. `548_patch_panel_census.py` — swaps the panel's ruminant farm densities to the census series (monitoring-series values preserved as `*_trend2019`). `546_census_beef_model.R` — exposure-source sensitivity (census time-varying / time-invariant / year-wise standardized vs monitoring series; Table S16).

## Environment

R 4.3 (INLA 24.x, spdep, sf, dplyr, ggplot2, patchwork, e1071, car), Python 3.11 (pandas, numpy, scipy, matplotlib, python-docx, PyWavelets, Pillow). INLA's Laplace approximations are deterministic, but newer INLA/BLAS builds may shift DIC/WAIC by a few points without changing any IRR, credible set, or classification.

## Notes

- Script numbering reflects the lab's internal pipeline ordering; scripts are self-contained and documented in-line (Korean comments; English column dictionaries above).
- Some scripts are generic across diseases (`450`, `451`) and take the disease name as an argument.

## License

MIT — see `LICENSE`.
