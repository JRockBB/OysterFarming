"""
    OysterFarming

Farm-scale simulation and decision model for eastern oyster aquaculture in Maine.
See `docs/README.md` for the research question and `docs/PLAN.md` for the work plan.
"""
module OysterFarming

using Arrow
using DataFrames
using Dates
using OrdinaryDiffEqTsit5
using TOML
using Tables

include("contract.jl")
include("params.jl")
include("deb.jl")

export ContractError, ColumnSpec, TableSpec, CONTRACT_SPECS
export validate_table, validate_sites, validate_forcing_daily, validate_closures,
       validate_growth_obs, validate_disease_obs, validate_economics, read_contract
export Param, ParamError, load_params, param_value
export DEBParams, arrhenius_factor, functional_response, salinity_factor, deb_fluxes,
       shell_height_mm, dry_tissue_g, gonado_somatic_index, spawns, initial_state,
       simulate_individual

end
