"""
    OysterFarming

Farm-scale simulation and decision model for eastern oyster aquaculture in Maine.
See `docs/README.md` for the research question and `docs/PLAN.md` for the work plan.
"""
module OysterFarming

using Arrow
using DataFrames
using Dates
using OrdinaryDiffEqTsit5: ODEProblem, Tsit5, init, solve, step!, derivative_discontinuity!, set_u!
using TOML
using Tables
using ForwardDiff
using LinearAlgebra: I, Symmetric, diag
using Optimization: Optimization, OptimizationFunction, OptimizationProblem, AutoForwardDiff
using OptimizationOptimJL: LBFGS
using Random: AbstractRNG, Xoshiro
using Statistics: mean, quantile
using Turing: @model, NUTS, MvNormal, Normal, Uniform, arraydist, sample, truncated

const SciMLBase = Optimization.SciMLBase

include("contract.jl")
include("params.jl")
include("deb.jl")
include("forcing.jl")
include("calibration.jl")
include("farm.jl")
include("nitrogen.jl")
include("dermo.jl")
include("stochastic.jl")
include("vibrio.jl")

export ContractError, ColumnSpec, TableSpec, CONTRACT_SPECS
export validate_table, validate_sites, validate_forcing_daily, validate_closures,
       validate_growth_obs, validate_disease_obs, validate_economics, read_contract
export Param, ParamError, load_params, param_value
export DEBParams, with_values, arrhenius_factor, functional_response, salinity_factor, deb_fluxes,
       shell_height_mm, dry_tissue_g, gonado_somatic_index, spawns, initial_state,
       simulate_individual, IndividualStepper, advance_day!
export forcing_for_site, continuous_blocks, longest_block, time_to_market
export DermoParams, dermo_division_rate, dermo_mortality_rate, dermo_transmission_cells, dermo_step,
       simulate_dermo, rank_to_log10_cells
export run_ensemble, farm_outcome
export VibrioRule, in_vibrio_control, vibrio_capacity, vibrio_handling_cost
export NitrogenParams, nitrogen_per_oyster_g, nitrogen_removed_g
export FarmParams, Cohort, HarvestRule, max_ingestion_ug_d, depleted_food, simulate_farm
export CalibrationProblem, predict_shell_height, sum_squared_error, fit_least_squares, fit_bayesian

end
