# OysterFarming

A farm-scale simulation and decision model for eastern oyster (*Crassostrea virginica*) aquaculture in Maine, built as thesis work and designed so that real farm data can replace public data later.

## The question

The Gulf of Maine is warming faster than almost any other ocean region. For Maine oyster farms this cuts both ways:

- **Benefit:** warmer water speeds growth. Today most of coastal Maine takes 2–5 years to grow an oyster to market size, so shorter grow-out is a large economic gain.
- **Cost:** warmer water also raises the risk of *Vibrio* (now regulated statewide for oysters each summer), parasitic disease (Dermo and MSX parasites are present in Maine waters, and MSX caused die-offs in the 1990s and 2010), and possibly more harvest closures.

**Core research question:** *Where and when along the Maine coast does the growth benefit of warming outweigh the added risk, and how should farms adjust seeding, density and harvest timing in response?*

## What the model does

One simulator, built in layers that share the same environmental drivers (temperature, salinity, food):

| Layer | Answers | Approach |
|---|---|---|
| Growth core | How fast do oysters reach market size at a site? | Dynamic Energy Budget (DEB) model, forced by temperature, salinity and food |
| Farm layer | How does stocking density limit growth? | Cohorts on a lease plus a simple food-depletion box model |
| Risk layer | How many oysters die, and how many harvest days are lost? | Dermo within-host model; closure model (rainfall, biotoxin, wastewater); Vibrio rules as harvest and handling constraints |
| Decision layer | What should a farm do? | Optimize seeding date, density and harvest schedule for profit, under present and warmer climates |
| Outputs | What else does the farm produce? | Nitrogen removal accounting on harvested biomass |

Ocean acidification is out of scope (it is a hatchery and larval-stage problem, not grow-out).

## Thesis structure

1. **Growth:** build, calibrate and validate the DEB growth core; run it with Maine forcing data.
2. **Risk:** add disease mortality and lost harvest days, treated stochastically.
3. **Decisions:** optimize farm management under climate scenarios; profit as the main objective, nitrogen removal as a second.

## Repository layout

The repository root is the Julia package `OysterFarming`; the Python pipeline and data live beside it.

```
OysterFarming/
├── Project.toml         ← Julia package (the model itself)
├── src/                 ← Julia source
├── test/                ← Julia tests (one file per source file)
├── data-pipeline/       ← Python: download and clean public data
├── data/
│   ├── raw/             ← downloaded files (gitignored)
│   ├── processed/       ← clean files that follow the data contract
│   └── fixtures/        ← small synthetic contract files shared by both test suites
├── params/              ← parameter files, every value with a citation
├── scripts/             ← runnable analyses (one per figure or table)
└── docs/
    ├── README.md        ← this file
    ├── PLAN.md          ← phased work plan with agent-ready tasks
    ├── TECH_STACK.md    ← language and library choices, with reasoning
    └── AGENTS.md        ← rules for coding agents working in this repo
```

**Why two languages:** the model is written in Julia; downloading and cleaning public ocean, satellite and GIS data is done in Python, where those tools are much stronger. The two halves only talk through files in `data/processed/` that follow a fixed schema (see PLAN.md, "Data contract"). See TECH_STACK.md for the full reasoning.

## Study region and data

- **Primary region:** midcoast and southern Maine (Damariscotta River and neighboring estuaries).
- **Disease calibration site:** Chesapeake Bay, which has long-term disease records that Maine lacks. The disease model is calibrated there, then run with Maine conditions.

Main public data sources: Maine DMR (closures GIS, aquaculture landings), NERRS (Wells Reserve water quality), regional buoys, Landsat and Sentinel-2 satellite products, the Chesapeake Bay Program DataHub, Maryland DNR Fall Oyster Survey, and the Add-my-Pet DEB parameter database. Details and links are in PLAN.md.

## Status

Phase 0 (scaffolding, data contract validators, CI) and the start of Phase 2 (parameter loader, standard DEB model) are implemented. Published parameter sets are incomplete and marked `TODO(source)`.

## Related work to build on, not duplicate

- **UMaine (Kiffney, Brady et al.):** an eastern oyster DEB model for Maine coupled to Landsat and Sentinel-2, validated across four estuaries, with a time-to-market web tool. This project's contribution is the risk, farm-management and climate-decision layers on top of growth, not growth alone. Contact this group early.
- **FARM model (Ferreira et al.):** farm-scale production, carrying capacity, economics and nitrogen balance.
- **ABC platform (Ferreira et al.):** couples growth models to disease (SEIR) models.
- **Lavaud et al. 2017:** eastern oyster DEB with salinity effects.
