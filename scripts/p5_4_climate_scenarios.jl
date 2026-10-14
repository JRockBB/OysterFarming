# P5.4: delta-change warming scenarios run through growth (P2.7), farm economics (P5.2/P5.3) and the
# Dermo screen (P4.3) from one script.
#
# Warming: add dT to every day's temperature (default +0, +0.5, +1, +2 degC; override with CLI arguments).
# Food (chlorophyll) and salinity are unchanged: this ignores any change in food or salinity with warming.
# Inputs are the same files as the phase scripts:
#   growth: params/deb_c_virginica.toml with calibrated X_K (params/deb_c_virginica_xk.toml), 20 mm seed, 76 mm market
#   economics: X_K posterior-quantile scenarios and the scenario inputs of scripts/p5_3_grid_search.jl
#   Dermo: published Ragone Calvo et al. (2001) model (screen only, see PLAN.md P4.2/P4.3)
# The farm table carries the P5.3 caveats (MAA-derived harvest capacity, horizon truncation with no value for
# unharvested stock, 2019 Maine costs); compare scenarios with each other, not as absolute values.
#
# Usage: julia --project scripts/p5_4_climate_scenarios.jl [dT1 dT2 ...]
# Writes outputs/p5_4/scenarios.md and .arrow.

using OysterFarming
using Arrow, DataFrames, Dates, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "outputs", "p5_4")
const WARMING_C = isempty(ARGS) ? [0.0, 0.5, 1.0, 2.0] : parse.(Float64, ARGS)
const SITES = ["dam_lower", "dam_upper"]
const SEED_MM = 20.0
const MARKET_MM = 76.0
const INIT_FOOD_DAYS = 14
# Farm scenario inputs, as in scripts/p5_3_grid_search.jl
# Weekly harvest capacity, derived from the MAA 2020 benchmarking Tables: employee hours per week
# (Table 13: 2.8 full-time x 38 h + 2.5 part-time x 23.1 h = 164.2 h) x 65 singles harvested per hour of labor
# (Table 42, established suspended farms) = 10,670 oysters/week. Table 13 covers all Maine shellfish and seaweed farms.
const WEEKLY_CAPACITY = 10_670.0
const N_COHORTS = 2
const RISK_Q = 0.1
const SEED_MONTHS = [4, 5, 7, 10]
const SEEDS = [100_000.0, 200_000.0, 400_000.0, 800_000.0]

P(file) = load_params(joinpath(ROOT, "params", file))
vir = P("deb_c_virginica.toml")
xk = P("deb_c_virginica_xk.toml")
deb = DEBParams(vir; X_K=xk["X_K"].value)
qs = sort([k for k in keys(xk) if startswith(k, "X_K_posterior_q")])
length(qs) == 5 || error("rerun scripts/p2_6_virginica_calibration.jl to export X_K posterior quantiles")
scenarios = [(0.2, DEBParams(vir; X_K=xk[k].value)) for k in qs]
farm = FarmParams(merge(P("farm_damariscotta.toml"), P("deb_c_virginica_clearance.toml")))
econ = EconParams(P("economics_maine_maa.toml"))
nitro = NitrogenParams(P("nitrogen_chesapeake_bmp.toml"))
vib = VibrioRule(P("vibrio_maine.toml"))
dermo = DermoParams(P("dermo_ragone_calvo2001.toml"))
forcing = read_contract(joinpath(ROOT, "data", "processed", "forcing_daily.arrow"), "forcing_daily")

warm(f, dT) = transform(f, :temp_c => (t -> t .+ dT) => :temp_c)

"Median years to market over seedings on the 15th of `month` in the block (runs stop at the first chl gap)."
function median_years_to_market(block, month)
    ys = Float64[]; n_not = 0
    for y in year(first(block.date)):year(last(block.date))
        seed = Date(y, month, 15)
        (seed < first(block.date) || seed + Day(INIT_FOOD_DAYS) > last(block.date)) && continue
        run = block[block.date .>= seed, :]
        gap = findfirst(ismissing, run.chl_ug_l)
        gap === nothing || (run = run[1:gap-1, :])
        nrow(run) > INIT_FOOD_DAYS || continue
        u0 = initial_state(deb; shell_height_mm=SEED_MM,
                           e=functional_response(mean(run.chl_ug_l[1:INIT_FOOD_DAYS]), deb.X_K))
        d = time_to_market(deb, run, seed, u0; market_shell_height_mm=MARKET_MM)
        d === nothing ? (n_not += 1) : push!(ys, d / 365.25)
    end
    return (isempty(ys) ? missing : median(ys)), length(ys), n_not
end

rows = DataFrame(site=String[], warming_c=Float64[], mean_temp_c=Float64[],
                 ttm_may_yr=Union{Missing,Float64}[], ttm_oct_yr=Union{Missing,Float64}[], ttm_not_reached=Int[],
                 best_seed_month=Int[], best_seeds=Float64[], best_expected_npv=Float64[], best_lower_npv=Float64[],
                 best_harvested=Float64[], best_nitrogen_kg=Float64[], dermo_division_days_yr=Float64[])
for site in SITES
    block0 = longest_block(forcing_for_site(forcing, site))
    for dT in WARMING_C
        block = warm(block0, dT)
        may, _, nm = median_years_to_market(block, 5)
        oct, _, no = median_years_to_market(block, 10)
        setup = DecisionSetup(block, scenarios, farm, econ, nitro, vib, WEEKLY_CAPACITY, 1.0, 0.0,
                              SEED_MM, N_COHORTS, RISK_Q)
        g = grid_search(setup; seed_months=SEED_MONTHS, seeds_per_cohort=SEEDS, market_sizes_mm=[MARKET_MM])
        b = first(g)
        t = collect(skipmissing(block.temp_c))
        yrs = unique(year.(block.date))
        div_days = mean(count((year.(block.date) .== y) .& coalesce.(block.temp_c .> dermo.div_T_min_c, false)) for y in yrs
                        if Date(y, 1, 1) >= first(block.date) && Date(y, 12, 31) <= last(block.date))
        push!(rows, (site, dT, mean(t), may, oct, nm + no, b.seed_month, b.seeds_per_cohort, b.expected_npv,
                     b.lower_npv, b.expected_harvested, b.expected_nitrogen_kg, div_days))
        println("done $(site) +$(dT)")
    end
end

mkpath(OUT)
Arrow.write(joinpath(OUT, "scenarios.arrow"), rows)
fmt(x; d=2) = ismissing(x) ? "not reached" : string(round(x; digits=d))
open(joinpath(OUT, "scenarios.md"), "w") do io
    println(io, "# P5.4 warming scenarios (generated by scripts/p5_4_climate_scenarios.jl)\n")
    println(io, "Delta-change warming on each site's longest forcing block; food and salinity unchanged. ",
            "Time to market: $(SEED_MM) to $(MARKET_MM) mm at calibrated X_K = $(round(deb.X_K; digits=3)) ug/L, median over seeding years. ",
            "Farm: best of seed months $(SEED_MONTHS) x seeds/cohort $(Int.(SEEDS)) at $(MARKET_MM) mm by expected NPV over ",
            "X_K posterior quantiles (P5.3 scenario inputs and caveats apply). Dermo: days/yr above the ",
            "$(dermo.div_T_min_c) degC division threshold (P4.3 screen).\n")
    println(io, "| site | warming degC | mean temp degC | years to market, May seed | years to market, Oct seed | runs not reaching market | best seed month | best seeds/cohort | expected NPV USD | q$(RISK_Q) NPV USD | harvested | N removed kg | Dermo division days/yr |")
    println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in eachrow(rows)
        println(io, "| $(r.site) | $(r.warming_c) | $(fmt(r.mean_temp_c)) | $(fmt(r.ttm_may_yr)) | $(fmt(r.ttm_oct_yr)) | ",
                "$(r.ttm_not_reached) | $(r.best_seed_month) | $(Int(r.best_seeds)) | $(fmt(r.best_expected_npv; d=0)) | ",
                "$(fmt(r.best_lower_npv; d=0)) | $(fmt(r.best_harvested; d=0)) | $(fmt(r.best_nitrogen_kg; d=1)) | ",
                "$(fmt(r.dermo_division_days_yr; d=1)) |")
    end
end
println(read(joinpath(OUT, "scenarios.md"), String))
