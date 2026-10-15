# P4.2f Maine runs with the symbolic transmission formula (passed the held-out Great Bay test narrowly; see
# outputs/p4_2f/symbolic.md and docs/PLAN.md P4.2f for its caveats).
#
# Design (fixed before running):
# - Formula (from scripts/p4_2f_symbolic.jl, constants as found by symbolic regression):
#   g(T, S) = max(0, exp(A - exp(T - S + C)) + D + E * S).
# - Refit on all data (Maryland + Great Bay, arrival 2003.15) of the five mechanistic parameters (log beta, log mu,
#   logit host, log epsilon MD, log epsilon GB) by Nelder-Mead, as in the per-fold refits of the symbolic step.
# - Maine: Damariscotta sites, complete years 2016-2019 cycled to 10 years (as the farm runs), +0/0.5/1/2 degC; parasite
#   present at the start with true prevalence 0.36 at light intensity (level 10, below RFTM detection; Marquis et al.
#   2020 qPCR 34-39 %); background import epsilon_GB. Outputs: fall RFTM prevalence, true prevalence, Mackin intensity,
#   annual Dermo mortality; daily Dermo mortality files for the farm runs (P5.4 "with Dermo" scenario).
#
# Usage: julia --project scripts/p4_2f_maine.jl
# Writes outputs/p4_2f/{maine_formula.md, maine_formula.arrow, dermo_daily_<site>_<dT>.arrow} and params/dermo_ude_fit.toml.

using OysterFarming
using Arrow, DataFrames, Dates, Statistics
using Optimization, OptimizationOptimJL

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "outputs", "p4_2f")
include(joinpath(@__DIR__, "p4_2f_common.jl"))

const A, C, D, E = 1.5866469010762194, 10.625032545240591, -0.7141608793198079, 0.07903541132975125
gform(T, S) = max(0.0, exp(A - exp(T - S + C)) + D + E * S)
const WARMING = [0.0, 0.5, 1.0, 2.0]
const MAINE_PREV0 = 0.36
const MAINE_LEVEL0 = 10

θ̂ = [parse(Float64, split(l, ",")[2]) for l in readlines(joinpath(OUT, "ude_weights_full.csv"))[2:end]]
w = nn_weights(θ̂)

function objective(φ)
    θ = vcat(φ, w)
    md = 0.0; nmd = 0; gb = 0.0; ngb = 0
    for s in data
        r = simulate_establishment(params_from(θ, s.region), s.T, s.S, initial(s), 0.0;
                                   intro_day=intro_day(ARRIVAL["full"], s), background_mortality_d=BACKGROUND_D,
                                   record=s.rec, transmission_mod=gform)
        e = sum(abs2, r.prev .- s.obs)
        s.region == "gb" ? (gb += e; ngb += length(s.obs)) : (md += e; nmd += length(s.obs))
    end
    return md / nmd + gb / ngb
end
φ = solve(OptimizationProblem((x, _) -> objective(x), θ̂[1:5]), NelderMead(); maxiters=800).u
θ = vcat(φ, w)
println("refit objective ", round(objective(φ); digits=4), " params ", round.(φ; digits=3))

full_years(b) = [y for y in unique(year.(b.date)) if count(year.(b.date) .== y) == daysinyear(y)]
pM = params_from(θ, "gb")
rows = DataFrame(site=String[], warming_c=Float64[], year=Int[], fall_prev=Float64[], fall_true_prev=Float64[],
                 fall_mackin=Float64[], dermo_mortality_yr=Float64[])
mkpath(OUT)
for site in ("dam_lower", "dam_upper")
    b = longest_block(forcing_for_site(forcing, site))
    f = cycle_years(b, full_years(b), 10)
    dates = f.date; yrs = unique(year.(dates))
    for dT in WARMING
        T = Float64.(f.temp_c) .+ dT; S = Float64.(f.salinity_psu)
        N0 = zeros(EST_N_LEVELS); N0[1] = DENSITY0 * (1 - MAINE_PREV0); N0[MAINE_LEVEL0] = DENSITY0 * MAINE_PREV0
        rec = [findfirst(==(Date(y, 11, 1)), dates) for y in yrs]
        r = simulate_establishment(pM, T, S, N0, 0.0; background_mortality_d=BACKGROUND_D, record=rec, transmission_mod=gform)
        for (i, y) in enumerate(yrs)
            push!(rows, (site, dT, y, r.prev[i], r.true_prev[i], r.mackin[i], sum(r.deaths[year.(dates) .== y]) / DENSITY0))
        end
        Arrow.write(joinpath(OUT, "dermo_daily_$(site)_$(dT).arrow"), DataFrame(date=dates, dermo_mortality_d=r.deaths ./ DENSITY0))
    end
end
Arrow.write(joinpath(OUT, "maine_formula.arrow"), rows)

fmt(x; d=3) = string(round(x; digits=d))
open(joinpath(OUT, "maine_formula.md"), "w") do io
    println(io, "# P4.2f Maine runs with the symbolic transmission formula (scripts/p4_2f_maine.jl)\n")
    println(io, "Formula g(T, S) = max(0, exp($(round(A; digits=3)) - exp(T - S + $(round(C; digits=3)))) $(round(D; digits=3)) + $(round(E; digits=4)) S). ",
            "Refit objective $(fmt(objective(φ); d=4)). Parasite present at the start (36 % light prevalence).\n")
    println(io, "Summary by site and warming (mean over years 3-10, after the start-up years):\n")
    println(io, "| site | warming degC | fall RFTM prevalence | fall true prevalence | Mackin | Dermo mortality yr^-1 |\n|---|---|---|---|---|---|")
    for g in groupby(rows, [:site, :warming_c])
        r = g[3:end, :]
        println(io, "| $(g.site[1]) | $(g.warming_c[1]) | $(fmt(mean(r.fall_prev))) | $(fmt(mean(r.fall_true_prev))) | ",
                "$(fmt(mean(r.fall_mackin); d=2)) | $(fmt(mean(r.dermo_mortality_yr); d=4)) |")
    end
    println(io, "\nBy year:\n\n| site | warming degC | year | fall RFTM prevalence | fall true prevalence | Mackin | Dermo mortality yr^-1 |\n|---|---|---|---|---|---|---|")
    for r in eachrow(rows)
        println(io, "| $(r.site) | $(r.warming_c) | $(r.year) | $(fmt(r.fall_prev)) | $(fmt(r.fall_true_prev)) | $(fmt(r.fall_mackin; d=2)) | $(fmt(r.dermo_mortality_yr; d=4)) |")
    end
end

open(joinpath(ROOT, "params", "dermo_ude_fit.toml"), "w") do io
    println(io, "title = \"Dermo establishment model with the symbolic transmission formula (PLAN.md P4.2f)\"")
    println(io, "notes = \"Transmission multiplier g(T, S) = max(0, exp(g_a - exp(T - S + g_c)) + g_d + g_e * S) on the waterborne infection rate. Passed the held-out Great Bay test by about 1 % (0.207 vs 0.209 persistence); fails in Maryland; its temperature dependence is implausible (highest in cold water). Use only as a 'with Dermo' scenario.\"\n")
    srcf = "Symbolic regression (SymbolicRegression.jl) of the trained UDE network: scripts/p4_2f_symbolic.jl, outputs/p4_2f/symbolic.md"
    srcp = "Refit on Maryland DNR and NH Fish & Game prevalence with the formula: scripts/p4_2f_maine.jl, outputs/p4_2f/maine_formula.md"
    logistic_(x) = 1 / (1 + exp(-x))
    for (k, v, u, s) in (("g_a", A, "1", srcf), ("g_c", C, "degC or ppt", srcf), ("g_d", D, "1", srcf), ("g_e", E, "ppt^-1", srcf),
                         ("beta", exp(φ[1]), "m^2 cell^-1 d^-1", srcp), ("mu", exp(φ[2]), "d^-1", srcp),
                         ("host", logistic_(φ[3]), "1", srcp), ("epsilon_md", exp(φ[4]), "d^-1", srcp),
                         ("epsilon_gb", exp(φ[5]), "d^-1", srcp), ("gb_arrival_year", ARRIVAL["full"], "yr", "P4.2e full fit, scripts/p4_2e_establishment_fit.jl"))
        println(io, "[$k]\nvalue = $(round(v; sigdigits=6))\nunits = \"$u\"\nsource = \"$s\"\n")
    end
end
println(read(joinpath(OUT, "maine_formula.md"), String))
