# Work plan

This is the working plan for OysterFarming. Each task has an ID, a short description and a "Done when" line so it can be handed to a coding agent as-is. Read README.md for the research question and TECH_STACK.md for language choices. Agents must also follow AGENTS.md.

Phases map to thesis chapters:

| Phase | Content | Thesis |
|---|---|---|
| 0 | Repo scaffolding and data contract | — |
| 1 | Data pipeline (Python) | supports all |
| 2 | Growth core (DEB) | Chapter 1 |
| 3 | Farm layer (cohorts, density, nitrogen) | Chapter 1/3 |
| 4 | Risk layer (disease, closures, Vibrio) | Chapter 2 |
| 5 | Decision layer and climate scenarios | Chapter 3 |
| 6 | Farm-data readiness | partnerships |

Phases 1 and 2 can run in parallel once Phase 0 is done.

---

## Data contract

The Python pipeline writes these files to `data/processed/`. The Julia model reads only these files. Changing a schema requires updating this section, the Python writer, the Julia reader and both validators in one change.

All tables are Arrow IPC (`.arrow`). Dates are ISO dates (UTC). Missing values are nulls, never sentinel numbers.

### `sites.arrow`
| column | type | notes |
|---|---|---|
| site_id | string | short unique id, e.g. `dam_lower` |
| name | string | human-readable |
| lat, lon | float64 | WGS84 |
| region | string | `maine` or `chesapeake` |
| growing_area | string, nullable | Maine DMR growing area code, for joining closures |
| notes | string, nullable | |

### `forcing_daily.arrow`
One row per site per day.
| column | type | units / notes |
|---|---|---|
| site_id | string | |
| date | date | |
| temp_c | float64 | water temperature, °C |
| salinity_psu | float64 | practical salinity |
| chl_ug_l | float64, nullable | chlorophyll-a, µg/L (food proxy) |
| pom_mg_l | float64, nullable | particulate organic matter, mg/L |
| source | string | e.g. `nerrs:WELINWQ`, `erddap:<dataset>` |
| qc_flag | int8 | 0 = observed, 1 = gap-filled, 2 = suspect |

### `closures.arrow`
| column | type | notes |
|---|---|---|
| growing_area | string | DMR code |
| area_id | string, nullable | sub-area if available |
| start_date, end_date | date | end_date null if open-ended |
| closure_type | string | `rainfall`, `biotoxin`, `wastewater`, `pollution`, `conservation`, `other` |
| source | string | |

### `growth_obs.arrow`
| column | type | notes |
|---|---|---|
| dataset_id | string | |
| species | string | `C_virginica` or `C_gigas` |
| site_id | string | |
| cohort_id | string | |
| date | date | |
| shell_height_mm | float64, nullable | |
| dry_tissue_g | float64, nullable | |
| n | int | individuals measured |

### `disease_obs.arrow`
| column | type | notes |
|---|---|---|
| site_id | string | |
| date | date | |
| pathogen | string | `P_marinus`, `H_nelsoni` |
| prevalence | float64 | 0–1 |
| intensity_mean | float64, nullable | |
| intensity_scale | string | e.g. `mackin_0_5` |
| n | int, nullable | null when the source gives no per-sample count |
| source | string | |

### `economics.arrow`
| column | type | notes |
|---|---|---|
| year | int | |
| region | string | |
| value_usd | float64 | |
| pounds | float64, nullable | |
| count | float64, nullable | number of oysters |
| source | string | |

---

## Phase 0: Scaffolding

**P0.1 Repo skeleton.** Create the layout in README.md, `.gitignore` (ignore `data/raw/`), and a license.
Done when: the tree exists and `git status` is clean.

**P0.2 Julia package.** The repository root is the Julia package `OysterFarming` (`Project.toml`, `src/OysterFarming.jl`, `test/runtests.jl`) with dependencies pinned via `[compat]`.
Done when: `julia --project -e 'using Pkg; Pkg.test()'` passes with a placeholder test.

**P0.3 Python project.** Create `data-pipeline/` with uv (`pyproject.toml`, `src/oyster_data/`, `tests/`).
Done when: `uv run pytest` passes with a placeholder test.

**P0.4 Data contract validators.** Python: `oyster_data.contract` with one validator per table above (column names, types, units ranges such as 0 ≤ prevalence ≤ 1, −2 < temp_c < 35). Julia: matching `validate_*` functions.
Done when: both validators accept a small hand-made valid fixture and reject fixtures with a missing column, a wrong type, and an out-of-range value. Fixtures live in `data/fixtures/`.

**P0.5 CI.** GitHub Actions running both test suites on push. Tests must not hit the network.
Done when: CI is green.

---

## Phase 1: Data pipeline (Python)

Each fetcher downloads to `data/raw/<source>/`, never overwrites without a flag, and has a separate clean step that writes contract files. Fetchers are tested against saved sample responses, not live servers.

**P1.1 Site registry.** Write `sites.arrow` with 3–5 starting Maine sites that have continuous water data nearby (Damariscotta River is the priority), plus Chesapeake disease-calibration sites chosen in P1.7.
Done when: file validates; each site lists its forcing source in `notes`.

**P1.2 NERRS fetcher (Wells Reserve).** Download water-quality data (temperature, salinity, and chlorophyll or fluorescence if present) from the NERRS Centralized Data Management Office. Check access requirements first (the web service may need registration); fall back to bulk CSV export.
Done when: daily-averaged `forcing_daily` rows for Wells stations validate; QC flags from the source are respected.
*Status (2026-10-11):* blocked on a human step. The CDMO web service answers `Invalid ip` until the caller's IP address is registered with CDMO. Either register the IP or download CSV exports from the CDMO Advanced Query System into `data/raw/nerrs/`; then the parser can be written against a real export.

**P1.3 Buoy fetcher (ERDDAP).** Use `erddapy` to pull temperature, salinity and chlorophyll from regional buoys near the chosen sites (start by searching the NERACOOS ERDDAP server; verify the URL).
Done when: `forcing_daily` rows for at least one Damariscotta-area site validate.
*Status (2026-10-11):* the NERACOOS ERDDAP is now `https://data.neracoos.org/erddap` and no longer serves the E05 Damariscotta buoy; it has only offshore shelf buoys (E01, B01, I01...), and the IOOS sensors ERDDAP has only tide gauges near Damariscotta. Damariscotta forcing instead comes from the Maine EPSCoR LOBO buoys (`oyster_data.sources.lobo`; nodes 65 lower river and 66 upper river at Mook Sea Farm; hourly temperature, salinity, chlorophyll, seasonal deployments 2015–2023). Done-when met with LOBO. An ERDDAP fetcher for offshore buoys is optional.

*Open issue (winter forcing):* no Damariscotta source covers January–February, so P1.4 cannot fill those months from site climatology and multi-year runs break each winter. Candidate fixes need a decision: Maine DMR Boothbay Harbor daily SST (temperature only, 1905–present), satellite products (P1.10), or data from the UMaine group.

**P1.4 Gap filling and merging.** Merge sources per site into one daily series. Fill short gaps (≤ 7 days) by linear interpolation and longer gaps with a day-of-year climatology from the same site; mark both with `qc_flag = 1`.
Done when: a per-site summary (coverage %, filled %) is written to `data/processed/forcing_report.md`; unit tests cover both fill methods.

**P1.5 Maine closures.** Pull closure polygons from the Maine DMR ArcGIS REST services (current closures and the archived biotoxin layer) and the conditional-area rainfall rules. Convert to `closures.arrow` keyed by growing area.
Sources: `https://gis.maine.gov/arcgis/rest/services/dmr/DMR_Public_Health_Shellfish_Closures/MapServer`; archived biotoxin metadata at `https://www1.maine.gov/geolib/catalog/metadata/ARCHIVE_BIOTOXIN.html`.
Done when: file validates; a report lists how many years of history exist per closure type (this answers an open question below).
*Status (2026-10-11):* blocked. Every `gis.maine.gov` ArcGIS endpoint returned HTTP 502, and the ARCHIVE_BIOTOXIN metadata URL now redirects to the GeoLibrary home page. ArcGIS Online mirrors only area boundaries (`BTX_Areas`) and a 2019 NSSP classification layer, with no closure history. Retry the state server; if history is still unavailable, ask DMR for archived closure records. *Update (2026-10-13):* the state server is back but the closure services now return "Token Required" (not public). Public NSSP classification layers exist for 2015-2025 on the DMR ArcGIS Online organization (`services1.arcgis.com/RbMX0mRVOFNTdLzd`, e.g. `MaineDMR_Public_Health_2025_NSSP_Classifications`, fields GA, GA_SECTION, NSSP, CA_TYPE, CA_DETAIL). Both Damariscotta sites are in growing area WQ, classified Approved in 2021-2025 (no conditional rainfall or seasonal rule), so their closures come only from biotoxin or emergency closures; that history is still not available.

**P1.6 Economics.** Extract Maine oyster landings (pounds, count, value by year) from DMR landings reports and historical tables (`https://www.maine.gov/dmr/fisheries/commercial/landings-program/historical-data`).
Done when: `economics.arrow` validates and matches the published PDF totals for every year extracted.
*Status (2026-10-11):* the DMR historical-data page returns 404. `economics.arrow` is built instead from NOAA Fisheries One Stop Shop (`oyster_data.sources.foss`, ACCSP data, 1954–2024, eastern oyster only), and tests check it against the raw API totals. Still to do for the original done-when: compare against DMR's published aquaculture harvest totals, and confirm whether FOSS pounds are meat or live weight.

**P1.7 Chesapeake disease data.** (a) Extract Dermo and MSX prevalence and intensity by bar and year from Maryland DNR Fall Oyster Survey reports (PDF tables, `pdfplumber`). (b) Pull matching temperature and salinity from the Chesapeake Bay Program DataHub (`https://data.chesapeakebay.net`) for the nearest stations. Choose 5–10 disease-monitoring bars with good station coverage and add them to `sites.arrow`.
Done when: `disease_obs.arrow` and Chesapeake `forcing_daily` rows validate; extraction is spot-checked against 3 report tables by hand and recorded in the test suite.
*Status (2026-10-13):* (a) done: `oyster_data.sources.mddnr_fall_survey` extracts Tables 3-4 of the 2024 report (all 43 bars, 1990-2024) into `disease_obs.arrow`, spot-checked against the report in tests. NA/ND are not written; `n` is null (no per-bar counts in the tables; contract column made nullable); rows are dated 1 November of the survey year (midpoint of the 2024 survey window) and `source` marks the date as year-only. Bar coordinates: all 43 from Smith & Jordan 1992 (MD DNR CBRM-OX-93-3, Table 1, fixed sampling positions, verified against the scans) are in `sites.arrow`; (b) water data for the bars not started.

**P1.8 French growth data (method testing only).** Download the IFREMER 26-year Pacific oyster growth and mortality data (Zenodo record 6536065; curated file `DataResco_clean.csv`) and the ECOSCOPA high-frequency environmental data (SEANOE, DOI 10.17882/86131).
Done when: `growth_obs.arrow` (species `C_gigas`) and matching `forcing_daily` rows validate for at least 3 French sites.
*Status (2026-10-13):* superseded. The IFREMER data record only whole-oyster live mass (shell included), not shell height or dry tissue. The model is instead calibrated directly on eastern oyster shell height and dry tissue from Chesapeake cohorts (P1.9, P2.6), so the Pacific oyster method test is dropped rather than changing the contract for total mass.

**P1.9 Eastern oyster growth data.** Collect any public *C. virginica* growth series (shell height over time with matching forcing). Candidates are listed under Open questions.
Done when: whatever is found is in `growth_obs.arrow`, or a short note in `docs/` says what was searched and what is missing.
*Status (2026-10-13):* done for Chesapeake: Kingsley-Smith triploid cohorts (Patuxent, Severn, York; 2005-2007) from the Oyster BMP compilation (Zenodo 10.5281/zenodo.11205266) are in `growth_obs.arrow`. Further validation candidates: Zenodo 11966672 (US Northeast farms incl. Levinton NY 2008-10) and USGS 10.5066/P1C8KQWH (Louisiana 2024-25).

**P1.10 (later) Satellite forcing.** Landsat Collection 2 surface temperature and Sentinel-2-derived chlorophyll for sites without buoys, via STAC (`pystac-client`). This is hard (narrow estuaries, atmospheric correction, chlorophyll algorithms). Before building it, ask the UMaine group whether their processed climatologies can be shared.
Done when: deferred until Phase 2 results show buoy coverage is insufficient.

---

## Phase 2: Growth core (Julia)

**P2.1 Parameter files.** Define a TOML schema in `params/` where every parameter has `value`, `units`, `source` (full citation) and optional `notes`. Write a loader that refuses parameters without a source.
Done when: loader tests pass, including the refusal case.

**P2.2 Standard DEB model.** Implement the standard DEB model for one oyster: state variables reserve E, structural volume V, reproduction buffer E_R (maturity E_H optional). Include Arrhenius temperature correction, a Holling type II functional response on food, a salinity effect on filtration following Lavaud et al. 2017, a spawning rule, and conversion from V to shell height via a shape coefficient. Starting parameters: Add-my-Pet entry for *C. virginica* (verify it exists) and Lavaud et al. 2017.
Done when: the model runs on a `forcing_daily` site and returns daily shell height, dry tissue weight and reproduction buffer.

**P2.3 Model tests.** At minimum: (a) temperature correction equals 1 at the reference temperature; (b) with constant food and temperature, growth approaches the analytical von Bertalanffy curve; (c) no state goes negative; (d) zero food produces shrinking reserve, not growth.
Done when: all pass.

**P2.4 Method test on Pacific oysters.** Implement the *C. gigas* parameter set (Pouvreau et al. 2006) and run it on the French data from P1.8. This tests the machinery, not the eastern oyster model.
Done when: a figure of observed vs simulated growth for each French site is produced by `scripts/p2_4_gigas_check.jl`, with RMSE reported.
*Status (2026-10-13):* superseded by direct eastern oyster calibration (see P1.8 status).

**P2.5 Calibration pipeline.** Bayesian calibration of a small subset of parameters (start with the half-saturation constant for food and one maintenance parameter) using Turing.jl, plus a fast least-squares alternative with Optimization.jl.
Done when: on synthetic data generated from known parameters plus noise, both methods recover the true values within their credible/confidence intervals.

**P2.6 Eastern oyster calibration and validation.** Calibrate on part of the P1.9 data and validate on the rest (hold out whole sites or years, not random points).
Done when: validation metrics and figures are produced by a script; limitations are written up in `docs/`.
*Status (2026-10-13):* done. `scripts/p2_6_virginica_calibration.jl` fits X_K on the Patuxent Kingsley-Smith triploid cohort (Lavaud et al. 2017 parameters otherwise unchanged) and holds out York. X_K = 0.56 ug chl-a/L (least squares); Bayesian median 0.99, 95 % 0.16-3.73. Held-out York shell height RMSE 1.84 mm (2.5 %); dry tissue under-predicted (~30 %). A two-parameter fit (X_K, p_Am) was unidentifiable and worse on York, so rejected. Limitations: Chesapeake food (~16 ug/L) is near saturation, so X_K is bounded mainly from above; triploid calibration applied to diploids; Patuxent shape misfit in 2007 (model keeps growing through an observed plateau); Severn not used (no in-river forcing 2005-07). No plotting package yet, so outputs are tables (`outputs/p2_6/metrics.md`).

**P2.7 Time-to-market.** Function returning days from seeding to a market shell height (default 76 mm, as a parameter), and a script running it for every Maine site.
Done when: a table of time-to-market by site and seeding month is produced. Sanity check: results should broadly match the published UMaine finding that most of the region takes 2–5 years.
*Status (2026-10-13):* done for the two Damariscotta sites (`scripts/p2_7_time_to_market.jl`, 20 mm seed, 76 mm market): with the calibrated X_K, 2.3-2.9 years at dam_lower and 1.8-2.1 years at dam_upper; at the upper end of the X_K interval market size is not reached within the ~4.3-year forcing window. Consistent with the UMaine 2-5 year finding at the calibrated value; the upper bound is open because X_K is weakly identified.

---

## Phase 3: Farm layer (Julia)

**P3.1 Cohorts.** A farm holds cohorts (seeding date, count, individual DEB state). Background mortality as a daily rate parameter.
Done when: a 3-cohort farm runs for 5 years and individual counts are conserved apart from mortality and harvest (tested).

**P3.2 Food depletion.** Box model for the water around the lease: dX/dt = (X_in − X)/τ − (total filtration by all oysters), where X_in is the forcing chlorophyll and τ is a flushing time parameter. Oysters feed on X, not X_in.
Done when: tests show (a) zero oysters gives X = X_in at steady state; (b) higher density lowers per-oyster growth monotonically.

**P3.3 Harvest.** Harvest rules: harvest when individuals reach market size, optionally limited by a weekly labor capacity.
Done when: harvest counts and sizes are logged per day; a test checks no oyster is harvested below market size.

**P3.4 Nitrogen accounting.** Nitrogen removed = harvested tissue (and optionally shell) dry weight × nitrogen fraction. Use the Chesapeake oyster BMP expert panel shell-height-to-weight regressions and nitrogen fractions, with citations in `params/`.
Done when: function tested on a hand-calculated example.

---

## Phase 4: Risk layer (Julia)

**P4.1 Dermo model.** Within-host parasite model driven by temperature and salinity, following the structure of Ragone Calvo et al. 2001 (Journal of Shellfish Research 20:231–241). Map infection intensity to a daily mortality rate.
Done when: the model reproduces qualitative behavior (higher intensity at higher temperature and salinity) in tests.
*Status (2026-10-13):* done. `src/dermo.jl` implements Ragone Calvo et al. 2001 (J. Shellfish Res. 20(1):231-241; no DOI; PDF from the user in `data/raw/literature/`, gitignored) with parameters in `params/dermo_ragone_calvo2001.toml`. Choices: eq. 1 is applied as the sum of flows (as printed it multiplies by P twice); Table 1's S<3 mortality 0.006 d^-1 is used (text says 0.003); salinity values on a bin edge go to the higher bin; eq. 2 and eq. 3 are used as published. Tests show faster build-up at higher temperature and salinity. Open: the paper gives no mapping from parasite density to oyster mortality; that needs another source before P4.6. The Bushek et al. 1994 rank-to-cells regression (`rank_to_log10_cells`) is available for P4.2 against the MD DNR intensity data.

**P4.2 Dermo calibration on Chesapeake.** Calibrate to `disease_obs.arrow` from P1.7, holding out some bars.
Done when: observed vs simulated prevalence/intensity figures and metrics exist for held-out bars.

**P4.3 Dermo in Maine.** Run the calibrated model with Maine forcing under present conditions and warming scenarios (P5.4). Report time to critical infection intensity, if reached.
Done when: a table by site and scenario exists. Caveat to document: Maine strains and hosts may differ from Chesapeake.

**P4.4 Closure model.** From closure history (P1.5), estimate the probability that a growing area is closed on a given day by month and closure type. For rainfall closures, add a model linking daily precipitation (NOAA GHCN-Daily, fetched by a new Python task) to closure events using each area's published rainfall threshold.
Done when: simulated closure-day counts per year match the historical distribution per area (tested on held-out years).

**P4.5 Vibrio rules.** Encode Maine's Vibrio control plan as constraints: the statewide oyster control period (June 1 to October 15 under the 2024 rules; verify current dates) and the time-to-temperature handling requirements, represented as a harvest-day labor cost or capacity limit. Optional extension: a Vibrio abundance index from water temperature using the FDA risk assessment relationship.
Done when: harvests inside the control period carry the configured cost or limit; tests cover the boundary dates.
*Status (2026-10-13):* done. `src/vibrio.jl` + `params/vibrio_maine.toml`: Maine DMR Ch. 115 (as amended 2024-03-19) control period June 1 - October 15 inclusive, statewide for oysters, calendar based (verified against the rule text). Harvests in the period are flagged; `vibrio_capacity` reduces weekly harvest capacity and `vibrio_handling_cost` adds a per-oyster cost. Both magnitudes are scenario inputs: no published compliance-cost estimate was found (only grower comments in the 2024 rulemaking record, e.g. ice over $150/week for ~2,000 lb/week). Boundary dates are tested. The rule is applied to all simulated years as the current regime (it took effect 2016 and became statewide in 2024).

**P4.6 Stochastic driver.** Monte Carlo wrapper that samples closure sequences, disease outcomes and parameter uncertainty, and runs the farm simulator many times. Uses a seeded random number generator for reproducibility.
Done when: running the same seed twice gives identical results; outputs are distributions of yield, mortality and lost harvest days.

---

## Phase 5: Decision layer (Julia)

**P5.1 Economics.** Revenue from harvest count × price (from P1.6, with a scenario range). Costs per oyster seeded and per harvest day (seed, gear, labor); source these from extension enterprise budgets and record citations.
Done when: net present value per lease is computed for a reference farm and checked by hand.

**P5.2 Decision variables and objective.** Decision variables: seeding month, seeding density, harvest size threshold. Objective: expected NPV (with an option for a risk-averse objective such as a lower quantile). Secondary output: nitrogen removed.
Done when: objective function runs end to end on one site.

**P5.3 Optimizer.** Start with a grid search over the decision variables (simple and transparent). Add a derivative-free optimizer from Optimization.jl only if the grid is too coarse.
Done when: the optimal decision per site is reported with its expected NPV and uncertainty.

**P5.4 Climate scenarios.** Delta-change warming: add +0.5, +1.0 and +2.0 °C to each site's temperature series (configurable). Document that this ignores changes in food and salinity.
Done when: scenarios run through Phases 2–5 from one script.

**P5.5 Headline analysis.** For each site and scenario: change in time-to-market, change in risk (mortality, lost harvest days), change in optimal decisions and NPV. This is the core "does warming help or hurt" result.
Done when: `scripts/p5_5_headline.jl` produces the tables and figures for Chapter 3.

---

## Phase 6: Farm-data readiness

**P6.1 Intake templates.** Simple spreadsheet templates for farms: growth checks (date, cohort, shell height sample), mortality counts, seeding and harvest records, and an optional temperature logger file.
Done when: templates exist in `docs/farm_templates/` with a one-page explanation for farmers.

**P6.2 Farm data ingest.** Python converter from those templates into contract files, with anonymized site ids.
Done when: a filled-in example converts and validates.

**P6.3 Site recalibration.** Rerun P2.6 calibration with farm data for that site and compare against the public-data version.
Done when: a script produces the comparison for any farm site.

---

## Open questions

1. **Eastern oyster growth data for validation.** The French data is Pacific oyster. Candidates for *C. virginica*: the UMaine group's validation datasets (ask), Louisiana growth series from Lavaud et al. (check data availability), Texas (Lebreton et al.), and later partner farms.
2. **Overlap with the UMaine time-to-market tool.** Confirm what it covers so this project's growth chapter adds something (risk, decisions, climate) rather than duplicating it.
3. **Closure history depth.** How many years do the DMR GIS layers cover? (Answered by P1.5.)
4. **Vibrio data.** Public Vibrio counts for Maine were not found; the model starts with rules only.
5. **Costs.** Need published enterprise budgets for Maine or Northeast oyster farms.
6. **Language.** Julia core is the default (TECH_STACK.md). Confirm before Phase 2 starts.
7. **Disease transfer.** Is a Chesapeake-calibrated Dermo model defensible for Maine? Plan to treat Maine results as risk scenarios, not predictions.
8. **Winter forcing for Maine sites.** The Damariscotta LOBO buoys never cover January–February (see P1.3 status). Options: Boothbay Harbor DMR daily SST for temperature, satellite products, UMaine climatologies, or a documented winter assumption. Needed before multi-year runs (P2.7).
9. **French growth data measures total live mass.** The IFREMER RESCO data (P1.8) record mean whole-oyster live mass in g (shell included), not shell height or dry tissue. Using it needs (a) a new `growth_obs` column such as `total_mass_g` (contract change) and (b) a sourced conversion from DEB state to total live mass. P1.8/P2.4 are paused until this is decided.
10. **Eastern oyster maturity formulation.** Lavaud et al. (2017) use maturity energy E_H with k_J and a gonad compartment; the current model uses the Pouvreau et al. (2006) volume threshold V_p. Decide whether to implement the Lavaud formulation (unblocks `params/deb_c_virginica.toml`).

## Key references

- Pouvreau et al. 2006, DEB model for *C. gigas*, Journal of Sea Research 56:156–167.
- Lavaud et al. 2017, DEB for *C. virginica* with salinity, Ecological Modelling 363:221–233.
- Kiffney, Brady et al., UMaine satellite-forced eastern oyster DEB for Maine (Aquaculture; WAS 2025 abstract: https://was.org/Meeting/Program/PaperDetail/166316).
- Ragone Calvo et al. 2001, *P. marinus* within-host model, Journal of Shellfish Research 20:231–241.
- Marquis et al. 2020, qPCR survey of oyster parasites in Maine, Pathogens 9(4):256.
- Ferreira et al., FARM model: https://www.marineagronomy.org/node/130
- IFREMER oyster time series data paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC9271045
- Maine DMR Vibrio control plan, Chapter 115 (2024 rules).
- Chesapeake oyster BMP expert panel reports: https://www.chesapeakebay.net
