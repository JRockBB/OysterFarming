# Tech stack: decisions and reasoning

**Decision:** Julia for the model, Python for data acquisition, files as the boundary between them. No Rust. No Oceananigans in the core plan.

## What the project actually needs

Most of the computation is small and repeated many times:

- A DEB model is 4–5 coupled ODEs per oyster cohort, run daily for 2–5 simulated years.
- Calibration means running that model thousands to millions of times (Bayesian or least-squares fitting).
- The decision layer wraps the whole simulator in an optimizer and Monte Carlo loops over stochastic risk and climate scenarios.
- The data side means pulling time series and satellite imagery from many web services and cleaning them.

So the model side needs: fast ODE solving, parameter estimation, optimization, uncertainty quantification. The data side needs: mature clients for ocean data servers, satellite catalogs and GIS services.

## Ecosystem comparison

| Need | Julia | Python | Rust |
|---|---|---|---|
| ODE solving | Excellent: DifferentialEquations.jl / OrdinaryDiffEq.jl | Good: SciPy `solve_ivp`, diffrax (JAX) | Usable: `diffsol` (BDF, Runge-Kutta, adjoints), but small ecosystem |
| DEB-specific code | No maintained animal DEB package; write it yourself (small) | `pydeb` (Bruggeman): simulation plus Bayesian inference; reads Add-my-Pet | None |
| Calibration / Bayesian | Excellent: Turing.jl, Optimization.jl, autodiff through solvers | Good: PyMC, emcee, NumPyro | Weak |
| Optimization | Excellent: JuMP.jl, Optimization.jl, Optim.jl | Good: SciPy, Pyomo | Weak |
| Speed of many small model runs | Native speed, no vectorization tricks needed | Needs NumPy vectorization, Numba or JAX to be fast | Native speed |
| Ocean data servers (ERDDAP, THREDDS) | Thin; NCDatasets.jl reads OPeNDAP | Excellent: `erddapy` (maintained by IOOS), xarray | None |
| Satellite (STAC, Landsat, Sentinel-2) | Exists: STAC.jl, Rasters.jl, ArchGDAL.jl; less documented | Excellent: pystac-client, odc-stac/stackstac, rioxarray, earthaccess | Minimal |
| GIS (ArcGIS REST, shapefiles) | GeoDataFrames.jl, ArchGDAL.jl | Excellent: geopandas, shapely, requests | Minimal |
| Plotting | Makie.jl (very good) | matplotlib, plotly | Minimal |
| Coding-agent reliability | Decent; more API mistakes, especially with older package versions | Best; most training data | Good for generic code, weak for this domain |

## Is Julia the right choice? Mostly yes, with two pushbacks

**Why Julia fits the model core.** There is no off-the-shelf Julia DEB package, but that matters less than it sounds: the standard DEB model is small, and you should write it yourself anyway so you understand and control every equation for the thesis. Where Julia clearly wins is what surrounds the model: calibration, sensitivity analysis, optimization and Monte Carlo over many scenarios, all at native speed with automatic differentiation through the ODE solver. That is the bulk of thesis chapters 2 and 3.

**Pushback 1: do not do the data side in Julia.** Julia's ocean-data and satellite clients exist but are thinner and less documented than Python's. Fighting them would cost weeks for no scientific gain. Python owns downloading and cleaning; Julia reads clean files.

**Pushback 2: coding agents write Python more reliably than Julia.** This is real but manageable. Mitigations: pin package versions, keep the Julia code plain (functions and structs, minimal metaprogramming), and require tests for every function (see AGENTS.md).

**When Python-only would be better instead:** if you want a single language, if your advisor or collaborators work in Python, or if you plan to reuse `pydeb`'s Bayesian machinery rather than write your own. That is a legitimate path; the plan in PLAN.md works in either language. Decide before Phase 2 starts and do not switch mid-thesis.

## Oceananigans.jl: not for the core plan

Oceananigans is a GPU fluid-dynamics model for idealized to high-resolution ocean simulations. The companion package OceanBioME.jl adds biogeochemistry and individual-based "active particles", and can run as a box model. It is excellent software, but it is the wrong tool here:

1. **Scale mismatch.** Growth, risk and decisions are site-level (zero-dimensional) problems. No fluid dynamics is needed.
2. **Maine estuaries are hard to simulate.** They are narrow, tidal and have intertidal flats. UMaine's own carrying-capacity work used FVCOM (an unstructured-grid coastal model that handles wetting and drying), and even that ran into open-boundary problems. Support for realistic regional setups in Oceananigans is still evolving; open-boundary schemes are active research as of 2026.
3. **Better substitutes exist.** For physical forcing beyond buoys and satellites, use existing model output: NECOFS (FVCOM for the Gulf of Maine) and NOAA's GoMOFS, both served over THREDDS.

**Where it could come back:** an optional, late extension on bay-scale carrying capacity, such as how food depletion spreads across many farms in one estuary. Even then, consider reading FVCOM or GoMOFS output first. If you do use it, OceanBioME's box-model mode is the natural entry point.

## Rust: no

`diffsol` is a credible ODE solver, but Rust has no ecosystem for ecological modeling, calibration or ocean data. Nothing in this project needs Rust's strengths (memory safety, systems programming, deployment).

## Chosen packages

**Julia (model, the `OysterFarming` package at the repository root):**
- OrdinaryDiffEq.jl: ODE solving
- Turing.jl: Bayesian calibration
- Optimization.jl (and Optim.jl): fitting and decision optimization
- DataFrames.jl, Arrow.jl, CSV.jl: reading the data contract files
- NCDatasets.jl: NetCDF and OPeNDAP if needed
- TOML (stdlib): parameter files
- Makie.jl (CairoMakie): figures
- Test (stdlib), Aqua.jl: testing

**Python (data, `data-pipeline/`):**
- uv: environment and dependency management
- pandas or polars, pyarrow: tables and Arrow files
- xarray, netCDF4: gridded data
- erddapy: buoy and ERDDAP data
- pystac-client, odc-stac or stackstac, rioxarray: Landsat and Sentinel-2
- geopandas, requests: Maine DMR ArcGIS services
- pdfplumber: extracting tables from Maryland survey PDFs
- pytest: testing

**File formats at the boundary:** Arrow IPC (`.arrow`) for tables, NetCDF (`.nc`) for gridded fields, TOML for parameters. Arrow is read natively by both pyarrow and Arrow.jl, which avoids format bugs between the two languages.

## References

- pydeb: https://github.com/jornbr/pydeb
- Add-my-Pet DEB database and DEBtool: https://www.bio.vu.nl/thb/deb/deblab/
- OceanBioME.jl (JOSS 2023): https://joss.theoj.org/papers/10.21105/joss.05669
- Oceananigans regional open-boundary thesis (2026): https://iris.sissa.it/handle/20.500.11767/151772
- UMaine FVCOM carrying-capacity project: https://umaine.edu/aquaculture/?p=442
- GoMOFS: https://catalog-beta.data.gov/dataset/regional-hydrodynamic-model-outputs-of-the-noaa-gulf-of-maine-operational-forecast-system-
- erddapy: https://github.com/ioos/erddapy
- STAC.jl: https://github.com/JuliaClimate/STAC.jl
- Rasters.jl: https://github.com/rafaqz/Rasters.jl
- diffsol: https://github.com/martinjrobins/diffsol
