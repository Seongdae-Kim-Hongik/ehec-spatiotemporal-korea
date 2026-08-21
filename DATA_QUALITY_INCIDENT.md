# Data-quality incident: unmeasured covariate years ingested as zeros

*Recorded August 2026, before submission.*

## What happened

The public source file for the district health indicators records years in which an indicator was **not collected as `0` rather than as missing**. The assembly step read those zeros as observed values, and they propagated into the analysis panel and the main model.

Because only *some* years were affected, each contaminated variable retained non-zero variance across districts. The resulting estimates therefore looked ordinary — no infinite values, no degenerate standard errors, no obvious warning signs — while being driven largely by an artifact of which years happened to be measured.

| Variable | Contamination | Years actually measured |
|---|---|---|
| Bars/pubs per 1,000 population | Absolute establishment counts used where a density column existed; 9 of 14 years zero | 2016-2019 |
| Handwashing with soap | 10 of 14 years zero | 2013-2020 (split across differently spelled columns) |
| Handwashing after toilet use | 10 of 14 years zero | 2013-2015, 2017-2020 |
| Basic-livelihood recipients | 10 of 14 years zero | 2015-2018 |
| Foreign residents | All years zero after a failed density conversion | — |
| Groundwater failing quality standards | 8 of 14 years zero | 2011-2016 |
| Septic-tank population | 1 of 14 years zero, plus a units error | 2011-2024 |

Seventeen of 46 covariates were affected in total.

## Why the existing checks did not catch it

Four safeguards were already in place. Each failed for a specific, instructive reason.

1. **Standardization guard.** The guard replaced non-finite standardized values with zero. Partially contaminated columns still had non-zero standard deviation, so they passed as normal.
2. **Estimate fingerprint.** Fully zeroed covariates produce a recognizable signature (IRR exactly 1.000 with a symmetric wide interval). Partial contamination produces ordinary-looking estimates instead.
3. **Coverage audit.** The audit counted missing values only. Because the contamination presented as `0` rather than `NA`, coverage was reported as complete — including for years that were never collected.
4. **High-dimensional negative screen.** Screening 2,329 candidate indicators tests which variables predict incidence; it does not test whether the variables already in the model were ingested correctly.

A fifth failure surfaced a day later: the recovery was applied to the CSV panel but not to a cached matrix (`ST_covariates_panel.rds`) that several supplementary scripts read instead. Those scripts continued to run on contaminated inputs until the cache was rebuilt as well.

## Effect on results

Two associations reported in earlier drafts did not survive recovery.

| | Before recovery | After recovery |
|---|---|---|
| Bars/pubs | 1.228, credible | 0.905 (0.749-1.092) — direction also reverses |
| Fiscal independence | 0.856, credible | 0.842 (0.676-1.050) |
| Groundwater use | — | 0.873 (0.779-0.978) — the only credible covariate |

Descriptive person, time, and place analyses were unaffected: they derive from notification and population data, not from the covariate panel.

## Gates introduced

These are now mandatory before any model is fitted.

| Gate | Check |
|---|---|
| Z1 | Per-year share of **zeros** (not just missingness); a year that is more than 90% zero halts the pipeline |
| Z2 | Inspect all source columns before selection — density versus absolute count, and consistent normalization within a domain |
| Z3 | Produce a units table; a variable whose unit cannot be stated is not understood well enough to model |
| Z4 | Value-range sanity (percentages within [0, 100], rates non-negative) |
| Z5 | Correlation with population size, area, and expected counts; r > 0.5 suggests a scale proxy |
| Z6 | Principal-component check on survey-derived indicators; if the first component explains more than 35% of variance, they are not interpretable as individual behaviors |
| Z7 | Audit for duplicated data caches; recovery must be applied to every copy, and each copy revalidated |
| Z8 | Confirm output ownership before rerunning a script, so one script does not overwrite another's results |

## What this changed in the paper

The determinants section moved from naming two determinants to reporting that no covariate meets a prespecified robustness threshold. To make that claim defensible rather than merely negative, the analysis plan was fixed in advance (see `PREREGISTRATION.md`) and evaluated across 400 model specifications rather than one.
