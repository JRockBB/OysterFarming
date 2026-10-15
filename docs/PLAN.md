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
| region | string | `maine`, `chesapeake` or `great_bay_nh` (Great Bay, New Hampshire: Dermo extrapolation test, added 2026-10-14) |
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

### `mortality_obs.arrow`
Observed oyster mortality for a population (added 2026-10-13 for P4.1/P4.2).
| column | type | notes |
|---|---|---|
| site_id | string | |
| date | date | end of the observation period |
| period_days | int, nullable | length of the period the fraction covers (e.g. 365 for annual) |
| mortality_fraction | float64 | 0–1 |
| method | string | `box_count`, `cohort_count`, ... |
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
*Status (2026-10-11):* blocked on a human step. The CDMO web service answers `Invalid ip` until the caller's IP address is registered with CDMO. Either register the IP or download CSV exports from the CDMO Advanced Query System into `data/raw/nerrs/`; then the parser can be written against a real export. *Update (2026-10-13):* optional, not pursued (user decision). Wells data for 1995-2011 is in the NCEI NERRS archive used by `sources/nerrs_ncei.py` (stations welinwq, welhtwq, wellmwq, welsmwq; nutrient files for chlorophyll) with no registration. Later years need the CDMO zip form, which emails files to the requester. The Damariscotta analysis does not need Wells.

**P1.3 Buoy fetcher (ERDDAP).** Use `erddapy` to pull temperature, salinity and chlorophyll from regional buoys near the chosen sites (start by searching the NERACOOS ERDDAP server; verify the URL).
Done when: `forcing_daily` rows for at least one Damariscotta-area site validate.
*Status (2026-10-11):* the NERACOOS ERDDAP is now `https://data.neracoos.org/erddap` and no longer serves the E05 Damariscotta buoy; it has only offshore shelf buoys (E01, B01, I01...), and the IOOS sensors ERDDAP has only tide gauges near Damariscotta. Damariscotta forcing instead comes from the Maine EPSCoR LOBO buoys (`oyster_data.sources.lobo`; nodes 65 lower river and 66 upper river at Mook Sea Farm; hourly temperature, salinity, chlorophyll, seasonal deployments 2015–2023). Done-when met with LOBO. An ERDDAP fetcher for offshore buoys is optional.

*Open issue (winter forcing):* no Damariscotta source covers January–February, so P1.4 cannot fill those months from site climatology and multi-year runs break each winter. Candidate fixes need a decision: Maine DMR Boothbay Harbor daily SST (temperature only, 1905–present), satellite products (P1.10), or data from the UMaine group.

**P1.4 Gap filling and merging.** Merge sources per site into one daily series. Fill short gaps (≤ 7 days) by linear interpolation and longer gaps with a day-of-year climatology from the same site; mark both with `qc_flag = 1`.
Done when: a per-site summary (coverage %, filled %) is written to `data/processed/forcing_report.md`; unit tests cover both fill methods.

**P1.5 Maine closures.** Pull closure polygons from the Maine DMR ArcGIS REST services (current closures and the archived biotoxin layer) and the conditional-area rainfall rules. Convert to `closures.arrow` keyed by growing area.
Sources: `https://gis.maine.gov/arcgis/rest/services/dmr/DMR_Public_Health_Shellfish_Closures/MapServer`; archived biotoxin metadata at `https://www1.maine.gov/geolib/catalog/metadata/ARCHIVE_BIOTOXIN.html`.
Done when: file validates; a report lists how many years of history exist per closure type (this answers an open question below).
*Status (2026-10-11):* blocked. Every `gis.maine.gov` ArcGIS endpoint returned HTTP 502, and the ARCHIVE_BIOTOXIN metadata URL now redirects to the GeoLibrary home page. ArcGIS Online mirrors only area boundaries (`BTX_Areas`) and a 2019 NSSP classification layer, with no closure history. Retry the state server; if history is still unavailable, ask DMR for archived closure records. *Update (2026-10-13):* the state server is back but the closure services now return "Token Required" (not public). Public NSSP classification layers exist for 2015-2025 on the DMR ArcGIS Online organization (`services1.arcgis.com/RbMX0mRVOFNTdLzd`, e.g. `MaineDMR_Public_Health_2025_NSSP_Classifications`, fields GA, GA_SECTION, NSSP, CA_TYPE, CA_DETAIL). Both Damariscotta sites are in growing area WQ, classified Approved in 2021-2025 (no conditional rainfall or seasonal rule), so their closures come only from biotoxin or emergency closures; that history is still not available. *Update (2026-10-13, later):* bacterial part done. `oyster-data closures-fetch` / `closures` (`sources/dmr_nssp.py`) write `closures.arrow` and `closures_report.md` from the public 2025 NSSP layer for the Maine sites' growing areas (WQ: 8 standing pollution closures since 2008; the WWTP conditional area has no published events). The public current-closures service (gis.maine.gov/mapservices/.../DMR_Public_Health_Current_Shellfish_Closures) holds active closures only, so biotoxin and emergency history is still missing; a search of web archives, DMR notices and toxin-monitoring data is under way. *Update (2026-10-14):* biotoxin and emergency history reconstructed from about 140 Internet Archive captures of DMR closure inventories, hotline pages, 100-A/100-B legal notices and newsletters (`sources/dmr_closure_events.py`). For growing area WQ: one dated event, an oil-spill emergency closure of the upper Damariscotta 2021-06-25 to 2021-06-28 (dam_lower outside; dam_upper on the boundary line). American oysters were never in a biotoxin closure naming the Damariscotta 2009-2021. Whether the 2025 PSP oyster closure (NH border to Pemaquid, late May to 2025-07-11) covered the Damariscotta is unknown; DMR records would settle it.

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
*Status (2026-10-13):* done. `src/dermo.jl` implements Ragone Calvo et al. 2001 (J. Shellfish Res. 20(1):231-241; no DOI; PDF from the user in `data/raw/literature/`, gitignored) with parameters in `params/dermo_ragone_calvo2001.toml`. Choices: eq. 1 is applied as the sum of flows (as printed it multiplies by P twice); Table 1's S<3 mortality 0.006 d^-1 is used (text says 0.003); salinity values on a bin edge go to the higher bin; eq. 2 and eq. 3 are used as published. Tests show faster build-up at higher temperature and salinity. Oyster mortality from Dermo: no published rate function was found (Powell, Klinck & Hofmann 1996 give only a lethal threshold, Mackin 5 ~ 2.3e6 cells/g wet). `scripts/p4_dermo_mortality_fit.jl` instead fits annual instantaneous observed mortality on mean Dermo intensity and MSX prevalence across the MD DNR disease bars (1990-2024, n = 1472): 0.090 yr^-1 per intensity unit (bar-bootstrap 95 % 0.074-0.105), MSX 1.05 yr^-1 per unit prevalence, R^2 0.30 (`params/dermo_mortality_mddnr.toml`). `dermo_excess_mortality_d` converts model density to intensity with the inverse Bushek et al. 1994 relation (assumes the MD 0-7 scale equals Bushek's rank scale) for use as `extra_mortality_d` in `simulate_farm`. Observational fit, so confounded by causes that co-vary with Dermo. The Bushek et al. 1994 rank-to-cells regression (`rank_to_log10_cells`) is available for P4.2 against the MD DNR intensity data.

**P4.2 Dermo calibration on Chesapeake.** Calibrate to `disease_obs.arrow` from P1.7, holding out some bars.
Done when: observed vs simulated prevalence/intensity figures and metrics exist for held-out bars.
*Status (2026-10-13):* run, negative result (`scripts/p4_2_dermo_validation.jl`, outputs/p4_2/validation.md). Forcing for 8 MD DNR bars, 1996-2024: CO-OPS temperature (Solomons, Tolchester, Cambridge, Bishops Head; other gauges as regression fills) and CBP DataHub salinity (CB5.1W, LE1.1, CB3.2, CB3.3W, ET5.2, EE3.1, EE3.0; the API allows at most 10 years per query, so requests are split). Test: one-year-ahead fall mean intensity from 1999 (earlier years used a less sensitive assay), whole bars held out. The published Ragone Calvo model overpredicts by about 3 intensity units (held-out RMSE 3.5 vs 1.0 for persistence) and only ranks years moderately (Spearman 0.55 held out). A division-rate multiplier fit ran to the grid edge (0.1) and is not used. Cause: with the paper's Bushek et al. (1994) conversion, mean intensities 1-3 are 1-4 cells/g, so a single transmission event (100-1000 cells/g) alone gives intensity 4-5. Hackett Point has no usable year pair (CB3.3W salinity gaps). Not done: the model is not validated for Maryland; P4.3 should use the published parameters only as an upper-bound screen (does Maine forcing allow proliferation at all), not as an intensity forecast. A better test needs individual-oyster intensity distributions (not bar means) or a different transmission formulation; see open question 11. *Update (2026-10-14), second model (P4.2b):* the Hofmann et al. (1995) / Powell et al. (1996) infection-class model is implemented in `src/dermo_powell.jl` (parasite and transmission parts as published, with typo corrections documented in `params/dermo_powell1996.toml`; host energetics not coupled; the Ragone Calvo model is kept). Test (`scripts/p4_2b_powell_validation.jl`, outputs/p4_2b/validation.md): continuous runs 1996-2024 on the same bars, scored on fall prevalence 1999-2024 (no 0-7 vs 0-5 intensity mapping is published, so intensity is rank-correlated only). The papers' standard 60 h halving time fitted best of their four values. Result: also negative. Model prevalence is about 0.9 at every bar and year (observed 0.28-0.83 bar means); held-out RMSE 0.46 vs 0.24 for persistence, Pearson -0.08. Cause: the published transmission term has no salinity dependence and saturates, so the model cannot produce the low prevalence at fresher bars (Swan Point, Oyster Shell Pt). Neither published Dermo model is validated for Maryland; Dermo stays a screen (P4.3). *Update (2026-10-14), statistical model (P4.2c):* `scripts/p4_2c_dermo_statistical.jl` (outputs/p4_2c/model.md) fits fall intensity and prevalence on the bars directly: model A = last year's value + degree-days above 20 degC + summer salinity + winter temperature; model B = climate only. Model A is the first Dermo model with held-out skill: intensity RMSE 0.95 vs 1.08 for persistence, prevalence RMSE 0.22 vs 0.24 (Pearson 0.82). Salinity is the dominant predictor. It cannot be applied to Maine: Damariscotta summer salinity (~30) is far above any Maryland bar-year (max 16.5) and its degree-days above 20 degC (5-97) far below (min 288). Data-based conclusion for Maine: even at +2 degC, Damariscotta summers stay at or below the coldest Maryland bar-year (dam_upper 0.89 of the Maryland minimum, dam_lower 0.17). No disease cost is added to the Maine NPV. *Update (2026-10-14), data expansion (phase 1 of the Dermo plan):* bar forcing now covers 42 of 43 MD DNR bars (`chesapeake_sites.BAR_SITES_EXTRA`: nearest long-record CBP station within 15 km, nearest CO-OPS gauges; Lancaster left out at 17 km). Regression fills are now clipped to the observed range at the site (gapfill.py), which removed sub-freezing extrapolations; this can shift other sites' filled days slightly. Rerun of P4.2c on 40 bars (872 bar-years), split 21/19 bars: model A held-out intensity RMSE 0.95 vs 1.12 persistence, prevalence 0.217 vs 0.239; it also beats persistence when each of the upper, middle and lower bay is held out in turn. The climate-only model B transfers poorly (lower-bay prevalence worse than the training mean). Northeast and lab data found (not yet in the pipeline): NH Fish & Game Great Bay RFTM Mackin counts 1996-2014 with NERRS sondes 1995-2011 (Great Bay buoy: 37-169 degree-days above 20 degC, summer salinity 21-30, i.e. Maine-like); Maine PCR prevalence only (Marquis et al. 2020; Damariscotta 34-39 % in 2016-17); lab infection experiments on Damariscotta oysters (Chu & Volety 1997) and others; no lab data at both below 10 degC and above 25 ppt.

**P4.2d Great Bay extrapolation test.** Add Great Bay, New Hampshire (region `great_bay_nh`: contract change in both validators and this file) as the out-of-region test for Dermo models: NH Fish & Game RFTM Mackin counts 1996-2014 for five natural beds (`sources/nh_greatbay.py`, `oyster-data greatbay-fetch`) and NERRS Great Bay sondes 1995-2011 as forcing (measured fills only; no winter data). Fit on Maryland, predict Great Bay prevalence, score against persistence and the Great Bay mean.
Done when: held-out Great Bay skill is reported for every candidate Dermo model.
*Status (2026-10-14):* data in the pipeline (67 bed-surveys, 5 beds; 18,725 forcing rows). Test run (`scripts/p4_2d_greatbay_test.jl`, outputs/p4_2d/greatbay.md), 43 Great Bay bed-years: every Maryland-fitted model fails. Climate-only model: RMSE 0.71, Pearson 0.07 (predicts ~1.0 everywhere: high salinity). Lag + climate: RMSE 0.54 vs 0.18 for persistence. Hofmann/Powell: RMSE 0.52, Pearson 0.18. Great Bay shows an invasion curve, not a climate response: prevalence near 0 in 1996-2003 at all beds, then 0.6-1.0 from 2004-2005, with no matching change in summer degree-days or salinity. Lab constraints point the same way (salinity ~30 not limiting; winter at 2 degC thins but does not clear parasites; northern isolates are not cold-adapted) and Ford & Smolowitz (2007) report the same temperature record giving 100 % prevalence at one Connecticut site and <= 51 % at another. Conclusion: in cool, salty water, whether and when Dermo is established matters more than summer warmth; a Maine Dermo model needs an establishment/transmission state, not only climate. Extracted literature data saved in data/literature/dermo/.

**P4.2e Dermo establishment model.** Infection classes with in-host growth fixed from lab rates, a waterborne parasite pool fed by shedding and deaths, fitted transmission (`src/dermo_establishment.jl`, `params/dermo_establishment_lab.toml`, `scripts/p4_2e_establishment_fit.jl`). Pass rule fixed before running: beat persistence on held-out Great Bay beds (survey-weighted); only then run Maine and couple Dermo mortality to the farm.
Done when: held-out results are recorded and, if it passes, Maine runs and the with/without-Dermo NPV exist.
*Status (2026-10-14):* fails the pass rule; no Maine runs or farm coupling. Fit on 47 sites / 1,141 surveys: beta 5.7e-10 m^2 cell^-1 d^-1, pool decay 0.50 d^-1, host factor at its upper bound (1.0), import 2.1e-4 d^-1 (MD) and 1.1e-3 d^-1 (Great Bay), Great Bay arrival 2003.2 (2002.2-2004.1 in the held-out refits; an earlier report showed '2000' from a rounding bug). Great Bay in-sample RMSE 0.19, Pearson 0.88. Leave one Great Bay bed out: the model beats or ties persistence on four of five beds (Nannie Island 0.13 vs 0.19, Woodman Point 0.24 vs 0.27, Adams Point 0.22 vs 0.22, Oyster River 0.18 vs 0.18; Pearson 0.89-0.95 on three) but loses badly on Squamscott River (0.43 vs 0.22), the freshest, least-infected bed, so the survey-weighted rule fails. It is the first model that tracks the Great Bay invasion curve on held-out beds. Maryland: worse than persistence and near the region mean in every held-out region (0.33-0.38 vs 0.21-0.25): one shared parameter set does not capture the salinity-driven year-to-year swings of an endemic system. Next (needs user decision): the SciML step (learned transmission/salinity term + symbolic regression), or site-level effects.

**P4.2f SciML step (UDE + symbolic regression).** A small network g(S, T) multiplies the waterborne transmission rate of the P4.2e model (`transmission_mod`), trained with LBFGS + ForwardDiff (`scripts/p4_2f_ude_fit.jl`); then symbolic regression distils g into a formula (`scripts/p4_2f_symbolic.jl`, environment `scripts/sciml`). Same data, folds and pass rule as P4.2e.
Done when: held-out results for the network and the formula are recorded.
*Status (2026-10-15), network:* fails. A first run collapsed the network to a constant (weight penalty 1e-3 larger than the data gradient, and a wrong Great Bay arrival year); rerun with penalty 1e-5 and arrival 2003.15. The network learned a strong salinity gate on transmission (about 0 below 12 ppt, rising steeply above 15-20 ppt, consistent with the salinity-dependent transmission of Ragone Calvo et al. 2001) and an implausible temperature shape (dip at 15-22 degC, high at 5 degC). In-sample error fell (0.157 -> 0.139; Maryland Pearson 0.29-0.46). Held out: Great Bay 0.249 vs 0.209 persistence (survey-weighted), worse than P4.2e on four of five beds; Maryland regions 0.32-0.36 vs 0.21-0.25 persistence. The extra flexibility overfits. *Formula (2026-10-15):* symbolic regression of the network over the training conditions (`outputs/p4_2f/symbolic.md`); chosen (simplest within 5 % of best loss, complexity 15): g = max(0, exp(1.587 - exp(T - S + 10.625)) - 0.714 + 0.0790 S). With mechanistic parameters refit per fold, held-out Great Bay 0.207 vs 0.209 persistence: passes the rule, by about 1 %. Wins on Nannie Island (0.149 vs 0.191), Oyster River (0.154 vs 0.178), Woodman Point (0.252 vs 0.269); ties Adams Point; loses on Squamscott (0.389 vs 0.224). Maryland regions still worse than persistence (0.31-0.35 vs 0.21-0.25). Caveats: the formula's form came from the full-data network (held-out data influenced its shape, not its constants); its temperature dependence is implausible (transmission highest in cold water: g = 6.5 at 5 degC and 30 ppt, 1.7 at 22-28 degC), contrary to the lab evidence; the salinity gate (none below ~10 ppt) is consistent with the network and Ragone Calvo et al. (2001). Used for Maine only as a 'with Dermo' scenario next to the existing results. *Maine and farm (2026-10-15):* `scripts/p4_2f_maine.jl` refits the five mechanistic parameters with the formula (`params/dermo_ude_fit.toml`) and runs the Damariscotta sites (parasite present at the start). It predicts ~0.9 fall RFTM prevalence and Mackin ~4 at both sites, far above the observed 2016-17 Maine values (qPCR 34-39 %, light, no die-offs), so it overpredicts; its Dermo mortality (1.7-5.6 % yr^-1 now, 5.2-9.3 % at +2 degC) is used as a pessimistic bound. `DecisionSetup` gained an optional `extra_mortality_d`; P5.4/P5.5 report a 'with Dermo' column. Upper-bound Dermo cost over the 10-year farm: dam_lower 9k (now) to 32k USD (+2 degC), dam_upper 36k to 62k USD; best decision unchanged (April-May, 800,000 seeds); warming still raises NPV at both sites with Dermo included (dam_lower 96k -> 165k, dam_upper 189k -> 246k USD).

**P4.3 Dermo in Maine.** Run the calibrated model with Maine forcing under present conditions and warming scenarios (P5.4). Report time to critical infection intensity, if reached.
Done when: a table by site and scenario exists. Caveat to document: Maine strains and hosts may differ from Chesapeake.
*Status (2026-10-13):* screen done (`scripts/p4_3_dermo_maine.jl`, outputs/p4_3/dermo_maine.md), not a forecast because P4.2 failed. Published model, parasite assumed present, warming +0/0.5/1/2 degC on 2016-2019 forcing. dam_lower: water above the 20 degC division threshold on about 2 days/yr, no parasite multiplication at present or +0.5 degC; multiplication starts at +1 degC (16 days/yr, 4x) and grows at +2 degC (38 days/yr, 18x). dam_upper: 58 days/yr already (35x), 90 days/yr at +2 degC. Hog Island (MD) for scale: 146 days/yr. Critical intensity (MD DNR high range >= 3.0, params/dermo_thresholds_mddnr.toml) is reached on the transmission day at every site, because the Bushek conversion turns the 1000 cells/g input alone into intensity ~5; the time-to-critical column is therefore uninformative. Usable conclusion: the upper estuary is warm enough for Dermo to multiply now; the lower estuary is near the threshold and crosses it with about 1 degC of warming. Caveat: Maine strains and hosts may differ.

**P4.4 Closure model.** From closure history (P1.5), estimate the probability that a growing area is closed on a given day by month and closure type. For rainfall closures, add a model linking daily precipitation (NOAA GHCN-Daily, fetched by a new Python task) to closure events using each area's published rainfall threshold.
Done when: simulated closure-day counts per year match the historical distribution per area (tested on held-out years). *Status (2026-10-14):* not built, because there is nothing to fit at the Maine sites. Both are in Approved water 2015-2025 with no rainfall-conditional rule, and the only dated closure found 2009-2025 is a 4-day emergency in 2021 that dam_lower was outside and dam_upper was at the edge of. Closure risk is therefore negligible next to the fixed Vibrio control window (P4.5) on present evidence; the 2025 PSP event is the open gap.

**P4.5 Vibrio rules.** Encode Maine's Vibrio control plan as constraints: the statewide oyster control period (June 1 to October 15 under the 2024 rules; verify current dates) and the time-to-temperature handling requirements, represented as a harvest-day labor cost or capacity limit. Optional extension: a Vibrio abundance index from water temperature using the FDA risk assessment relationship.
Done when: harvests inside the control period carry the configured cost or limit; tests cover the boundary dates.
*Status (2026-10-13):* done. `src/vibrio.jl` + `params/vibrio_maine.toml`: Maine DMR Ch. 115 (as amended 2024-03-19) control period June 1 - October 15 inclusive, statewide for oysters, calendar based (verified against the rule text). Harvests in the period are flagged; `vibrio_capacity` reduces weekly harvest capacity and `vibrio_handling_cost` adds a per-oyster cost. Both magnitudes are scenario inputs: no published compliance-cost estimate was found (only grower comments in the 2024 rulemaking record, e.g. ice over $150/week for ~2,000 lb/week). Boundary dates are tested. The rule is applied to all simulated years as the current regime (it took effect 2016 and became statewide in 2024).

**P4.6 Stochastic driver.** Monte Carlo wrapper that samples closure sequences, disease outcomes and parameter uncertainty, and runs the farm simulator many times. Uses a seeded random number generator for reproducibility.
Done when: running the same seed twice gives identical results; outputs are distributions of yield, mortality and lost harvest days.

---

## Phase 5: Decision layer (Julia)

**P5.1 Economics.** Revenue from harvest count × price (from P1.6, with a scenario range). Costs per oyster seeded and per harvest day (seed, gear, labor); source these from extension enterprise budgets and record citations.
Done when: net present value per lease is computed for a reference farm and checked by hand.
*Status (2026-10-13):* done as functions; reference-farm run pending P5.2. `src/economics.jl` (`EconParams`, `farm_cashflows`, `npv`) with `params/economics_reference_farm.toml`: price 0.7785 USD/oyster (Maine DMR 2024 farm value / count); costs from Engle & van Senten, Economic Analysis of Oyster Production in Maryland (NOAA IR noaa/37917, Tables 11-12, 2018 dollars): seed 0.009 USD, non-seed variable cost 0.3119 USD per harvested oyster, fixed 4,800 USD/yr, capital 155,375 USD, discount rate 10 %. Hand-calculated NPV example is a test. Caveats: Maryland 2018 costs, not inflation-adjusted; no published Maine enterprise budget with these items was found (Maine Aquaculture Association benchmarking gives cost of goods sold per oyster but its source table was not verified). *Revision (2026-10-14):* Maine economics in `params/economics_maine_maa.toml` from the Maine Aquaculture Association benchmarking study (Engle et al. 2020, 2019 survey, established suspended farms): price 0.72 USD, size tiers petite 63.5-76.2 mm 0.641 and legal >= 76.2 mm 0.740 USD (Massachusetts 2013 price ratio and mix, UMass Dartmouth 2015), costs 0.531 USD per harvested oyster (Table 27), capital 40,395 USD per planted acre at 214,985 seeds per acre (Table 24), straight-line salvage over a 10-year life (Maryland Table 12). `src/economics.jl` gained `price_at` and `capital_outlay`. The Maryland file is kept.

**P5.2 Decision variables and objective.** Decision variables: seeding month, seeding density, harvest size threshold. Objective: expected NPV (with an option for a risk-averse objective such as a lower quantile). Secondary output: nitrogen removed.
Done when: objective function runs end to end on one site.
*Status (2026-10-13):* done. `src/decision.jl`: `FarmDecision` (seed month, seeds per cohort, market size), `evaluate_decision` gives expected NPV, a lower-quantile NPV, expected nitrogen removed and harvest over weighted parameter scenarios (X_K posterior quantiles 0.1-0.9, 20 % each, exported by P2.6). Runs end to end on dam_upper (`scripts/p5_3_grid_search.jl`).

**P5.3 Optimizer.** Start with a grid search over the decision variables (simple and transparent). Add a derivative-free optimizer from Optimization.jl only if the grid is too coarse.
Done when: the optimal decision per site is reported with its expected NPV and uncertainty.
*Status (2026-10-13):* grid search done (`grid_search`). First dam_upper run: best expected NPV at 1.6 million seeds per cohort, 65 mm, seeding in October/April/May (about 0.19-0.20 million USD over the 2015-2020 block). Not a result yet: it is driven by scenario inputs without sources - weekly harvest capacity (10,000 oysters, binding), one price for all sizes (favours the smallest market size), a 4.3-year horizon with no value for unharvested stock, and Maryland 2018 costs. Needs size-specific prices and a Maine capacity/labour source before interpretation. *Revision (2026-10-14):* with Maine economics and an MAA-derived weekly capacity of 10,670 oysters (Table 13 employee hours x Table 42 oysters per labor hour), the optimum is interior: May seeding, 400,000 seeds per cohort, 76 mm, expected NPV about +0.9k USD over the 2015-2020 block (near breakeven, consistent with MAA net margins of 0.08-0.12 USD per oyster). Remaining caveats: horizon truncation and capacity derived from survey averages. *Revision 2 (2026-10-14), 10-year horizon:* `cycle_years` (src/forcing.jl) repeats the four complete measured years (2016-2019) in order to a 10-year farm equal to the capital life (no salvage credit), with 8 annual cohorts. dam_upper optimum: May seeding, 800,000 seeds per cohort, 76 mm; expected NPV 224k USD, 10th-percentile 98k USD (1.6 million seeds is worse, so the optimum is interior). Assumes the measured years are representative weather. *Economic sensitivity (2026-10-14, `scripts/p5_3b_economic_sensitivity.jl`, outputs/p5_3b/sensitivity.md):* price 0.72 +/- 0.12 USD (MAA Table 33) x labor +/- 0.21 USD per oyster (MAA Table 27), 9 variants per site. The decision is not robust: the margin per oyster drives both profitability and the best scale. With +1 SD labor cost every decision loses money (best -100k to -210k USD); with -1 SD price at mean labor both sites are near or below zero; with wider margins the best decision moves to 1.6 million seeds at 65 mm (grid edge). Keeping the baseline decision (May, 800,000, 76 mm) costs up to about 250k USD of NPV when the margin differs. A farm's own price and labor cost should set the scale before this model is used for advice.

**P5.4 Climate scenarios.** Delta-change warming: add +0.5, +1.0 and +2.0 °C to each site's temperature series (configurable). Document that this ignores changes in food and salinity.
Done when: scenarios run through Phases 2–5 from one script.
*Status (2026-10-13):* done with one gap (`scripts/p5_4_climate_scenarios.jl`, outputs/p5_4/scenarios.md; warming values are CLI arguments). +0/0.5/1/2 degC on the 2016-2020 blocks. Warming speeds growth: May-seeded time to market at dam_upper falls from 2.0 to 1.3 yr and at dam_lower from 2.35 to 2.0 yr. Best expected NPV rises (dam_upper 87k to 190k USD, dam_lower -51k to 63k USD; P5.3 caveats apply, compare scenarios only). Dermo division days rise (dam_upper 58 to 90/yr, dam_lower 2 to 38/yr). Gap: Dermo is not coupled into farm mortality because the model failed validation (P4.2), and closures (P4.4) are blocked, so the NPV gains are upper bounds on the net effect of warming.

**P5.5 Headline analysis.** For each site and scenario: change in time-to-market, change in risk (mortality, lost harvest days), change in optimal decisions and NPV. This is the core "does warming help or hurt" result.
Done when: `scripts/p5_5_headline.jl` produces the tables and figures for Chapter 3.
*Status (2026-10-13):* done with gaps (`julia --project=scripts/figures scripts/p5_5_headline.jl`; outputs/p5_5/headline.md and headline.png). `scripts/figures` is a separate environment holding CairoMakie so the package does not depend on it. Result: warming helps growth and economics at both Damariscotta sites (May-seeded time to market -0.3 yr at dam_lower and -0.7 yr at dam_upper by +2 degC; best expected NPV +104k to +114k USD; same best decision, October seeding at the highest density, at every level) and raises Dermo risk (division days 2 to 38/yr at dam_lower, 58 to 90/yr at dam_upper). The Vibrio control period is a fixed calendar window (137 days), so it does not change. Gaps: Dermo cost is not in the NPV (P4.2), closures are not modelled (P4.4 blocked), and the P5.3 economic caveats apply, so NPV gains are upper bounds. *Revision (2026-10-14):* rerun with Maine economics. dam_lower loses money at every warming level (best -6.9k to -1.6k USD, 100,000 seeds); dam_upper goes from +0.9k to +17.5k USD at +2 degC (400,000 seeds, May). Growth and Dermo conclusions unchanged. *Revision 2 (2026-10-14):* 10-year horizon. Best decision at both sites and all warming levels: April-May seeding, 800,000 seeds per cohort, 76 mm. Expected NPV dam_lower 102k -> 194k USD (+0 -> +2 degC; 10th percentile -71k -> +3k), dam_upper 224k -> 308k USD (10th percentile 98k -> 206k). Disease: the statistical model (P4.2c) cannot be applied to Maine; Damariscotta summers stay at or below the coldest Maryland bar-year even at +2 degC.

---

## Phase 6: Farm-data readiness

**P6.1 Intake templates.** Simple spreadsheet templates for farms: growth checks (date, cohort, shell height sample), mortality counts, seeding and harvest records, and an optional temperature logger file.
Done when: templates exist in `docs/farm_templates/` with a one-page explanation for farmers.
*Status (2026-10-13):* done. `docs/farm_templates/`: farm_info, growth_checks (one row per oyster), mortality_counts, seeding, harvest and temperature_logger CSV templates, with a one-page README for farmers.

**P6.2 Farm data ingest.** Python converter from those templates into contract files, with anonymized site ids.
Done when: a filled-in example converts and validates.
*Status (2026-10-13):* done. `oyster-data farm-ingest FOLDER OUT` (`data-pipeline/src/oyster_data/farm_ingest.py`) writes sites, growth_obs (cohort-date means with n) and mortality_obs contract tables, plus daily logger temperature (days with >= 12 hours of readings), seeding and harvest CSVs. Site id is farm_ + 8 hex digits of SHA-256(farm name, site label); names are not written. The synthetic example in data-pipeline/tests/fixtures/farm_example_synthetic converts and validates (Python and Julia). Logger temperature is not a forcing_daily table because that contract requires salinity; P6.3 merges it into public forcing.

**P6.3 Site recalibration.** Rerun P2.6 calibration with farm data for that site and compare against the public-data version.
Done when: a script produces the comparison for any farm site.
*Status (2026-10-13):* done (`scripts/p6_3_site_recalibration.jl FARM_DIR PUBLIC_SITE [--triploid]`). Fits X_K per farm cohort on public forcing with farm logger temperature substituted, and compares shell-height RMSE under the public X_K, the cohort's own fit and the other cohorts' median fit. Pipeline check with SYNTHETIC data (`--synthetic-demo`, true X_K 1.5 at dam_upper, noise 1.5 mm): fits 1.41 (1.28-1.55) and 1.44 (1.21-1.66); out-of-sample RMSE 1.1-1.8 mm vs 5.2-5.5 mm with the public X_K. No real farm data yet.

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
11. **Dermo model does not reproduce Maryland bar intensities (P4.2).** Options: (a) obtain MD DNR per-oyster intensity scores (Cooperative Oxford Laboratory; not in the published report) to compare distributions instead of bar means; (b) refit transmission inputs as well as division, which needs a source or an accepted calibration design; (c) use a different published Dermo model (e.g. Powell, Klinck & Hofmann 1996). Option (c) was implemented (P4.2b) and also failed; Dermo results remain a screen only. Remaining options: (a) per-oyster data, (b) fitting transmission terms with an accepted calibration design, or a statistical model of prevalence on temperature and salinity fitted to the MD DNR bars.

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
