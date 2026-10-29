# EHEC Spatiotemporal Analysis - South Korea, 2001-2024

Reproducible analysis code for:

> **Changing Epidemiology of Enterohemorrhagic *Escherichia coli* Infections in South Korea, 2001-2024: Age-Standardized Bayesian Spatiotemporal Analysis.** Kim S, Chun BC. Submitted to *Journal of Korean Medical Science* (2026).

> **Note on version 1.1 (August 2026).** During pre-submission checks we found that several district covariates had been
> ingested with unmeasured years recorded as zero rather than missing, which produced spurious associations. All covariates
> were recovered from source, every affected model was refitted, and a preregistered multiverse was added. Two associations
> reported in earlier drafts (bars/pubs, fiscal independence) do not survive; groundwater use is the only covariate whose
> credible interval excludes 1 in the main model. The incident, its root causes, and the resulting data-quality gates are
> documented in [`DATA_QUALITY_INCIDENT.md`](DATA_QUALITY_INCIDENT.md); the analysis plan fixed before estimates were
> inspected is in [`PREREGISTRATION.md`](PREREGISTRATION.md).


> **Note on version 2.0 (August 2026).** The covariate stage was rebuilt from a full census of every district-level indicator available in the source files (4,305 columns), with 108 indicators prespecified in 7 transmission-pathway domains and 85 entering a 3,000-specification multiverse (scripts 520-534). The final 20-covariate model, period analyses, and figures were re-estimated (540-544). The beef-cattle and dairy exposure was re-sourced from the quinquennial Census of Agriculture, Forestry and Fisheries (2010/2015/2020 district farm counts; scripts 545-548) so that measured values precede the study window; the quarterly livestock monitoring series is retained as an independent-source sensitivity analysis. The per-indicator adjudication file is included as `COVARIATE_PREREG_POOL.csv`.

> **Note on version 2.1 (October 2026).** Every district join now runs on the five-digit administrative code (`sigungu_cd`) through the key and alias tables in `gate/`; joins on district names were retired. The steps that silently filled missing values (minimum-value and interpolation fills of expected counts, mean fills of covariates) were removed and replaced by explicit, logged fills, and a model script now stops when a covariate it uses has missing cells. Cheongju 2011-2013 (before the 2014 city-county merger, absent from the KDCA district cross-tabulation) is treated as missing rather than as zero. Every panel must pass `panel_gate` (229 districts x all study years, finite values, positive population and expected counts) before a model is fitted. The resident-registration population and the census farm-density denominators were rebuilt, the latter on boundary-shapefile areas. All models, the 3,000-specification multiverse, and the supplementary tables were re-estimated. In the final 20-covariate model, beef-cattle farm density has IRR 1.264 (95% CrI 1.082-1.478; v2.0: 1.315), fiscal independence 0.785 (0.657-0.939), and handwashing with soap 1.463 (1.125-1.905); reservoir area and current smoking no longer exclude 1. In the multiverse, beef-cattle farm density is credible in 99.5% of specifications and fiscal independence in 25% (v2.0: 70%). The Table S3 check that adds the 0-4-year population share was refitted on the v2 panel (IRR 0.936, 0.768-1.142). Key result tables are included under `results/`.

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
export EHEC_GATE=/path/to/repo/gate  # (version 2.1) district key, alias table and panel gate; defaults to $EHEC_ROOT/gate
```

The first two default to the current directory (`.` and `./data`) if unset. Scripts whose first path line defines `GD` resolve the analysis folder as `$EHEC_ROOT/G.Downloads/FBD_SpatioTemporal` and the shared map-style helper as `$EHEC_ROOT/G.Downloads/EHP_Spatio/fig_style.R` (not included), mirroring the original folder layout. Required R packages: `INLA`, `sf`, `spdep`, `dplyr`, `tidyr`, `stringr`, `e1071`, `car`. Required Python packages: `pandas`, `numpy`, `matplotlib`, `python-docx`.

## Pipeline (run in order)

### Step 0 - Data assembly
- `450_agestd_assemble_generic.py` - builds the analysis panels: district x year x age-band cases, indirectly age-standardized expected counts (E_std), national crude/age-standardized trends.

### Step 1 - Descriptive and temporal analyses (2001-2024) → Figures 1, 3; eFigures 2-15
- `380_sexage_descriptive.py`, `388_kdca_enrich2.py`, `389_kdca_enrich3.py` - sex- and age-specific descriptive epidemiology (eFigures 2-4).
- `370_stl_national.R` - STL decompositions, crude and age-standardized (eFigures 5-6).
- `462_person_rate_season.py` - age-specific incidence rates and seasonal concentration (eFigure 9).
- `483_age_rate_ratio.py` - child-to-older rate ratios + NB GLM age-by-year interaction test (Figure 1C, eFigure 12).
- `490_child_share_test.py` - demography-only null decomposition of the falling child share (eFigure 15).
- `486_summer_share.py` - July-August share; Gini vs uniform and constant-seasonality null floors (Figure 3A, eFigure 13).
- `383_wavelet_split.py` - Morlet CWT, early/late halves (eFigure 8).
- `487_wavelet_changetest.py`, `489_wavelet_multimethod.py`, `489b_amp_null.py` - change-point, permutation, AR(1)-surrogate, count-matched null, harmonic-amplitude tests (eTable 24, eFigure 14).
- `492_null_envelope.py` - yearly annual-band share vs count-matched constant-seasonality envelope (Figure 3B).

### Step 2 - Spatial and spatiotemporal models (2011-2024) → Figure 2; Tables 1-2; eTables 3-23
- `451_agestd_model_generic.R` - generic age-standardized NB BYM2 spatiotemporal engine (R-INLA).
- `418_ehec_table1_house.R` - main 16-covariate model (Table 2), posterior spatial RR, exceedance.
- `412_agestd_maps.R`, `412b_jama_maps.R`, `415_agestd_exceedance.R` - maps (Figure 2A, eFigure 1).
- Emerging hot spot analysis (Getis-Ord Gi* + Mann-Kendall on EB-smoothed age-standardized ratios) - Figure 2B, eTable 19.
- `484_spatial_4yr_split.R`, `485_period_interaction.R` - calendar-window splits + formal covariate-by-period interaction (eTables 22-23, eFigure 16).
- `498_covid_split.R`, `499_covid_map.R` - pre- vs post-COVID-19 split (2011-2019 vs 2020-2024) + interaction test and maps (eFigure 17, eTable 25).
- Robustness: `430_agestd_supp_robust.R`, `322_8graph.R` (eight neighborhood graphs), `323_M1M6_ladder.R` (model ladder), `290_RSR_sensitivity.R` (restricted spatial regression), `292_ZINB_sparse.R` / `296_ZINB_refit_final.R` (zero-inflation), `293_PIT_histograms.R` / `297_randomized_PIT.R` (calibration), `463_fullscreen_transform.py` (2,329-indicator negative screen).

### Step 3 - Publication figures
- `493_jama_fig1_person.py`, `494_jama_fig2_place.py`, `496_jama_fig2_time.py` - assemble Figures 1-3.

### Step 4 - Covariate recovery, preregistered multiverse, and behavioral-domain tests (version 1.1)
- `503_full_covariate_audit.py` - quality audit of all 46 covariates: per-year share of zeros (not just missingness), value ranges, unit table, and confounder correlations.
- `500_recover_covariates.py`, `504_recover_handwashing.py`, `505_recover_remaining.py` - recovery of covariates whose source columns recorded uncollected years as zero, including density-versus-count corrections and coalescing of indicators split across differently spelled columns.
- `515_recover_rds_cache.R` - rebuilds the cached covariate matrix used by the supplementary model scripts, which the CSV-level recovery did not reach.
- `506_multiverse.R`, `507_specification_curve.py` - preregistered multiverse over 400 specifications (domain-representative covariate sets x ruminant-exposure form x transformation x spatial structure x offset x adjacency graph), with decision rules fixed before any estimate was inspected; specification curve (eFigure 16, eTable 26).
- `508_detection_adjusted.R` - reestimation with case-ascertainment proxies entered jointly (eTable 27).
- `510_domain6_health_behaviour.py`, `512_d6_urbanicity_test.R` - health-behavior and preventive-care domain (influenza vaccination, infant check-ups, nutrition-label literacy, walking, smoking, drinking), tested singly and jointly (eTable 28).
- `511_kchs_common_factor.py` - principal-component diagnosis of survey-derived indicators, which share a dominant rural-urban component.
- `513_soap_forensics.py`, `514_soap_signflip.R` - sequential-adjustment and confounding checks for the handwashing-with-soap association (eTable 29).


### Step 5 - Full covariate census, preregistered pool, and final model (version 2.0)
- `520_covariate_triage.py`, `521_covariate_manual_adjudication.py` - census of all 4,305 district-level columns; rule-based triage plus per-indicator manual adjudication against the transmission framework (output: `COVARIATE_PREREG_POOL.csv`).
- `530_build_prereg_panel.py` - assembles the 108-indicator district-year panel with the data-quality gate (per-year zero-share audit, denominator integrity, renormalization prohibition); `531_z6_vif_design.py` - survey-indicator PC diagnosis and within-domain VIF pruning (85 modeled).
- `532_prereg_stage1_domain.R`, `533_prereg_multiverse.R`, `534_prereg_mv_summary.py` - 3,000-specification stratified multiverse (12 parallel workers) and its summary (Tables S13-S15).
- `540_final_main_model.R` - final 20-covariate BYM2+RW1+type I model (Table 2, Table S11, maps); `542_final_period_analyses.R` - 4-year windows, COVID-19 split, and period-interaction tests (Tables S8-S9, S17); `543_final_spec_curve.py`, `544_det_windows_final.py` - Figures S16-S17.
- `545_census_livestock.py` - parses the 2010/2015/2020 agricultural-census district livestock files, maps them to the 229-district frame (with administrative-boundary aliases), and interpolates 2011-2024; a national-total assertion guards against double counting. `548_patch_panel_census.py` - swaps the panel's ruminant farm densities to the census series (monitoring-series values preserved as `*_trend2019`). `546_census_beef_model.R` - exposure-source sensitivity (census time-varying / time-invariant / year-wise standardized vs monitoring series; Table S16).

### Step 6 - Code-keyed district joins, explicit missingness, and the v2.1 refit (version 2.1)
- `gate/` - district join and panel gate module shared by the lab's spatiotemporal series: `sigungu_key.csv` (229 five-digit codes with canonical labels and shapefile areas), `sigungu_alias.csv` (2,329 observed spellings, merger histories and ambiguous or sub-district names that are refused), `sigungu_key.py` / `sigungu_key.R` (`to_code`, `attach_code`, `explicit_fill`, `read_shp_checked`), and `panel_gate.py` / `panel_gate.R` (the pre-model gate). Rules are summarised in `gate/README_시군구결합_규칙.md` (Korean).
- `_v2_common.R` - common rules for every model script: expected counts, cases and covariates are joined on `sigungu_cd`, no minimum-value or interpolation fill of expected counts, standardisation that stops on missing covariates (`v2_need`, `v2_pool_check`), and the Cheongju 2011-2013 mask (`CJ_MASK=0` turns it off for the sensitivity analysis).
- `v2_fixes.py` - builder-side module: code-keyed joins, shapefile-area denominators, raw-source handwashing indicators, external KOSIS series for the Jeju fiscal indicators and the Sejong gross regional product, and the explicit-fill log.
- `530b_v2_prereg_fix.py` - corrections applied to the preregistered panel after `530` to `531` (survey composite PC1 for all 229 districts, structural zeros of the septic-tank indicator, Sejong and Jeju source fills, all logged).
- `_st_setup.R`, `_hm_setup.R` - setup scripts sourced by the spatial engines, now based on `_v2_common.R`.
- `400_agestd_assemble.py`, `400c_ehec_national_rebuild.py`, `450_agestd_assemble_generic.py`, `450c_national_trend_rebuild.py` - population, expected-count and national-trend builders re-keyed on the administrative code (December resident-registration population per code and year, duplicate notification rows dropped).
- `530_build_prereg_panel.py`, `545_census_livestock.py`, `548_patch_panel_census.py`, `531_z6_vif_design.py` - panel builders with shapefile-area denominators and assertions instead of fills.
- Re-estimated model scripts: `540_final_main_model.R` (Table 2, Table S11, Table S22), `540cj_final_main_model.R` (Cheongju row-removal sensitivity), `542_final_period_analyses.R`, `542w_final_windows_3cut.R`, `542x_ehec_window_interaction.R`, `542c_covid_phi_ci.R` (Tables S8-S9, S17, S20-S21, Figure 4, Figure S15-S16), `553_final_ladder_zinb.R`, `553z_nb_vs_zinb_irr.R` (Tables S18-S19), `546_census_beef_model.R`, `546b_monitor_beef_sens.R` (Table S16), `532_prereg_stage1_domain.R`, `533_prereg_multiverse.R`, `534_prereg_mv_summary.py`, `543_final_spec_curve.py` (Tables S12-S15, Figure S17), `411_agestd_S6.R`, `411b_ehsa_esri.R`, `412c_mapB_esri.R` (emerging hot spot analysis, yearly Moran's I, Figure 3A-3B, Tables S7, S21), `595_ehec_yearly_persistence.R`, `595b_ehec_f46_redraw.R` (persistence, Kendall's W, Figure 3C, Table S23), `415_agestd_exceedance.R` (Figure S14), `451_agestd_model_generic.R`, `499_covid_map.R`, `592_ehec_window_grid.py`, `544_det_windows_final.py`.
- `401b_children_pct_v2.py` - district 0-4 and 0-9 year population shares (2008-2024) from the monthly resident-registration extract, code-keyed with pre-merger units summed, gated before saving (input to the Table S3 refit, `v2d_rerun/s3_child04_refit_v2d.R`).
- `v2d_rerun/` - the scripts that drove the October 2026 re-run: `run_chain_v2d.sh` (the ordered chain with the 12-seed multiverse queue), `mv_pool_v2d.sh` (six-worker multiverse pool), `finalize_after_mv.sh` and `post_finalize_261029.sh` (post-multiverse table, figure and document refresh), `540s_sensitivity.R` (prespecified sensitivity set: PC1 removed, Jeju indicator, Cheongju observed zeros; Table S24), `33_persistence_jkms_100k_rerun.py` (100,000-draw persistence null model and seed stability; Table S23), `s15_axis_robustness.py` (Table S15), `s3_children_pct_v2d.py` and `s3_child04_refit_v2d.R` (Table S3 refit with the 0-4-year population share), `fig3c_f46_v2d.R` and `fig4_grid_v2d.py` (Figure 3C and Figure 4 redraws), and the supplementary-document updaters (`update_table2_v2c.py`, `update_supp_tables_v2c.py`, `update_supp_figures_v2c.py`, `mv_section.py`, `patch_na_counts.py`). Relative paths inside these scripts refer to the lab's working folders (`_v2_rerun_261022/JKMS`, `JKMS_FBD/JKMS_EHEC_v2d_261029`) under `EHEC_ROOT`.
- `results/` - the v2.1 result tables that feed the manuscript: `MASTER_final_table2.csv` and `final_fit.txt` (Table 2, Table S11), `MASTER_final_4yr_split.csv`, `MASTER_final_4yr_fitstat.csv`, `MASTER_final_3win_interaction.csv` (Tables S17, S20, Figure 4), `MASTER_final_covid_*.csv` (Tables S8-S9, S21), `MASTER_census_beef_check.csv`, `MASTER_monitor_beef_sens.csv` (Table S16), `MASTER_agestd_regionSIR.csv` (Table S22, Figure 3A), `MASTER_final_negcontrol.csv` (Table S14), `agestd_moran_yearly.csv`, `EHSA_agestd_cat_esri.csv`, `EHSA_agestd_sum_esri.csv` (Tables S7, S21, Figure 3B), `MASTER_agestd_supp_ladder.csv`, `MASTER_agestd_supp_zinb.csv`, `MASTER_supp_nb_vs_zinb_*.csv` (Tables S18-S19), `PREREG_MV_summary.csv` and `S15_axis_robustness_v2c.csv` (Tables S12-S15), `EHEC_T34c_kendallW.csv`, `JKMS_persistence_S23_B100k_seed20260913.csv`, `JKMS_persistence_seed_stability_B100k.csv` (Table S23), `S24_sensitivity_v2c.csv` (Table S24), and `S3_children_refit_table.csv` with `S3_children_refit_fit.txt` (Table S3). Per-seed multiverse checkpoints, cached panels and figures are not distributed.

## Environment

R 4.3 (INLA 24.x, spdep, sf, dplyr, ggplot2, patchwork, e1071, car), Python 3.11 (pandas, numpy, scipy, matplotlib, python-docx, PyWavelets, Pillow). INLA's Laplace approximations are deterministic, but newer INLA/BLAS builds may shift DIC/WAIC by a few points without changing any IRR, credible set, or classification.

## Notes

- Script numbering reflects the lab's internal pipeline ordering; scripts are self-contained and documented in-line (Korean comments; English column dictionaries above).
- Some scripts are generic across diseases (`450`, `451`) and take the disease name as an argument.

## License

MIT - see `LICENSE`.
