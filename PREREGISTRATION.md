# Preregistered analysis plan — district-level determinants

*Fixed August 2026, before any incidence rate ratio from the recovered covariate panel was inspected.*

## Rationale

Searching a large covariate space for statistically significant associations is data dredging, and the earlier covariate set had been shaped by that risk: it was assembled while contaminated variables appeared credible (see `DATA_QUALITY_INCIDENT.md`). Rather than repeat a significance-guided search on recovered data, we enumerate the defensible specifications in advance, fit all of them, and report **the share of specifications in which each covariate is credible**. The question shifts from "what is significant" to "how sensitive is the conclusion to specification choice".

## Domains

Covariates are grouped by the transmission pathway they are meant to represent, following the EHEC literature: ruminant reservoir, contamination of agricultural land and produce, water and sanitation, person-to-person and food-handling exposure, and socioeconomic context.

| Domain | Covariates |
|---|---|
| D1 Ruminant reservoir | beef cattle, dairy cattle, combined ruminant density (tertiles), pasture area |
| D2 Agricultural land use | upland field, paddy field, forest area |
| D3 Water and sanitation | piped supply, public sewerage, septic-tank population, groundwater use, household groundwater, groundwater failing standards, wastewater discharge, advanced treatment, organic load |
| D4 Personal hygiene and food handling | handwashing (after going out, before meals, with soap, after toilet use), fast-food and eating-out outlets, bars/pubs |
| D5 Socioeconomic and demographic | sex ratio, fiscal independence, basic-livelihood recipients, elderly living alone, foreign residents, population density, rural and urban population share |
| D6 Health behavior and preventive care *(amendment R7)* | influenza vaccination, infant health check-ups, nutrition-label literacy, walking, smoking, drinking |

Swine and poultry density are retained as **negative controls**: they are livestock exposures but not EHEC reservoirs, so a credible association with either would indicate residual confounding rather than a reservoir effect.

Children aged 0-4 and elderly share are excluded from the main model because age is already absorbed by the age-standardized offset; including them would adjust for age twice.

## Multiverse axes

Every combination of the following is fitted (400 specifications; 247 converged).

| Axis | Levels |
|---|---|
| Covariate set | domain-representative subsets |
| Ruminant exposure form | continuous per-SD, tertiles, ordinal trend |
| Transformation | raw, log1p when skewed |
| Spatial structure | BYM2, ICAR, unstructured |
| Offset | age-standardized expected count, population |
| Adjacency graph | queen, rook, k-nearest neighbors |

## Decision rules (fixed before inspecting estimates)

| Verdict | Criterion |
|---|---|
| Robust | Credible in at least 50% of specifications **and** consistent in direction |
| Suggestive | Credible in 20% to 50% |
| Null | Credible in under 20% |

A covariate that exceeds 50% but reverses direction across specifications does **not** qualify as robust; it is reported as unstable.

## Amendment R7 — health-behavior domain

The original plan excluded influenza vaccination coverage and self-rated health as case-ascertainment proxies. That exclusion rested on analyst judgment rather than evidence, and the same indicators can be read as preventive-behavior exposures — a reading under which the expected direction is opposite. Because the two interpretations are distinguishable by direction, the domain was reinstated and tested, with directional hypotheses recorded before estimation. Decision rules were not changed.

Amendments are recorded by addition only; nothing in this document is deleted after estimates are inspected.

## Handling of survey-derived indicators

Indicators from the community health survey share a dominant common component that tracks the rural-urban gradient (42.0% of variance). They are therefore not interpreted as individual behaviors, and any association among them is assessed with and without adjustment for that gradient before it is reported.
