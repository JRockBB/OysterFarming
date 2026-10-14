# P4.2b: test the Hofmann et al. (1995) / Powell et al. (1996) Dermo model on Maryland DNR disease bars.
#
# Design (fixed before running):
# - Continuous run per bar from 1996-10-01 to 2024-11-30 on the bar forcing (P4.2). Days with no forcing row
#   repeat the last available temperature and salinity; a survey year is scored only if fewer than 30 of the
#   365 days before it were repeated.
# - Population: market-size oysters only (size class 7, 76-88.9 mm), 20 oysters m^-2 with 50 % prevalence at the
#   start (Powell et al. 1996, p. 143; infected oysters start at one cell). Deaths are replaced by uninfected
#   oysters (constant density; a stand-in for recruitment, not from the papers). Non-Dermo mortality 0.0302 yr^-1
#   (intercept of the Maryland mortality fit, outputs/p4_dermo_mortality/fit.md).
# - Scored: fall survey years 1999-2024 (rows dated 1 November; MD DNR changed assays in 1999), so 1996-1998
#   is spin-up. Main metric: prevalence (thioglycollate-detectable, > 2^12 cells), which needs no intensity
#   conversion. Intensity: the model gives Mackin 0-5 and MD DNR reports a 0-7 scale with no published mapping,
#   so only rank correlation (Spearman) is reported.
# - Baseline: persistence (prevalence this fall = last fall's observed value).
# - Calibration: parasite halving time at 20 degC / 20 ppt chosen from the papers' own values {30, 60, 80, 120} h
#   by prevalence squared error on the calibration bars (same bars as P4.2); held-out bars report skill.
#
# Usage: julia --project scripts/p4_2b_powell_validation.jl
# Writes outputs/p4_2b/{validation.md, predictions.arrow}.

using OysterFarming
using Arrow, DataFrames, Dates, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "outputs", "p4_2b")
const CALIB = ["mddnr_hog_island", "mddnr_broome_island", "mddnr_swan_point", "mddnr_sandy_hill"]
const HOLDOUT = ["mddnr_hackett_point", "mddnr_oyster_shell_pt", "mddnr_sharkfin_shoal", "mddnr_goose_creek"]
const HALVING_H = [30.0, 60.0, 80.0, 120.0]
const START, STOP = Date(1996, 10, 1), Date(2024, 11, 30)
const FIRST_YEAR = 1999
const MAX_HELD_DAYS = 30
const MARKET_CLASS = 7
const DENSITY0 = 20.0
const PREV0 = 0.5
const BACKGROUND_YR = 0.0302

ptoml = load_params(joinpath(ROOT, "params", "dermo_powell1996.toml"))
forcing = read_contract(joinpath(ROOT, "data", "processed", "forcing_daily.arrow"), "forcing_daily")
disease = read_contract(joinpath(ROOT, "data", "processed", "disease_obs.arrow"), "disease_obs")
dermo = disease[disease.pathogen .== "P_marinus", :]

"Daily forcing over START:STOP with gaps filled by repeating the last value; `held` marks filled days."
function bar_forcing(site)
    f = forcing_for_site(forcing, site)
    have = Dict(r.date => (Float64(r.temp_c), Float64(r.salinity_psu)) for r in eachrow(f))
    dates = collect(START:Day(1):STOP)
    first_obs = minimum(keys(have))
    last = have[first_obs]
    T = Float64[]; S = Float64[]; held = Bool[]
    for d in dates
        if haskey(have, d)
            last = have[d]; push!(held, false)
        else
            push!(held, true)
        end
        push!(T, last[1]); push!(S, last[2])
    end
    return DataFrame(date=dates, temp_c=T, salinity_psu=S, held=held)
end

function initial_state(p)
    N = zeros(length(p.W), POWELL_N_LEVELS)
    N[MARKET_CLASS, 1] = DENSITY0 * (1 - PREV0)
    N[MARKET_CLASS, 2] = DENSITY0 * PREV0
    return N
end

forc = Dict(s => bar_forcing(s) for s in vcat(CALIB, HOLDOUT))

function run_bar(site, halving)
    p = PowellDermoParams(ptoml; halving_time_h=halving)
    f = forc[site]
    sim, _ = simulate_powell_dermo(p, f, initial_state(p); background_mortality_d=BACKGROUND_YR / 365)
    obs = sort(dermo[dermo.site_id .== site, :], :date)
    rows = DataFrame(site_id=String[], year=Int[], prev_obs=Float64[], prev_last=Union{Missing,Float64}[],
                     prev_model=Float64[], int_obs=Union{Missing,Float64}[], mackin_model=Float64[], held_days=Int[])
    for (i, r) in enumerate(eachrow(obs))
        y = year(r.date)
        y >= FIRST_YEAR || continue
        idx = findfirst(==(r.date), sim.date)
        idx === nothing && continue
        held = count(f.held[max(1, idx - 364):idx])
        held < MAX_HELD_DAYS || continue
        prev_last = i > 1 && year(obs.date[i-1]) == y - 1 ? Float64(obs.prevalence[i-1]) : missing
        push!(rows, (site, y, Float64(r.prevalence), prev_last, sim.prevalence[idx],
                     ismissing(r.intensity_mean) ? missing : Float64(r.intensity_mean), sim.mean_mackin[idx], held))
    end
    return rows
end

ranks(x) = invperm(sortperm(x))
spearman(a, b) = length(a) > 2 ? cor(ranks(a), ranks(b)) : NaN
rmse(a, b) = sqrt(mean((a .- b) .^ 2))

sse = Dict(h => sum(sum((r.prev_model .- r.prev_obs) .^ 2) for r in (run_bar(s, h) for s in CALIB)) for h in HALVING_H)
h_best = HALVING_H[argmin([sse[h] for h in HALVING_H])]
res = vcat([transform(run_bar(s, h_best), [] => (() -> s in CALIB ? "calibration" : "held out") => :role) for s in vcat(CALIB, HOLDOUT)]...)
std60 = vcat([run_bar(s, 60.0) for s in HOLDOUT]...)

mkpath(OUT)
Arrow.write(joinpath(OUT, "predictions.arrow"), res)
fmt(x) = (x isa Real && isnan(x)) ? "n/a" : string(round(x; digits=3))
open(joinpath(OUT, "validation.md"), "w") do io
    println(io, "# P4.2b Hofmann/Powell Dermo model on Maryland DNR bars (scripts/p4_2b_powell_validation.jl)\n")
    println(io, "Continuous runs 1996-10 to 2024-11, scored on fall surveys $(FIRST_YEAR)-2024. Halving time chosen on ",
            "calibration bars from $(Int.(HALVING_H)) h: **$(Int(h_best)) h** (prevalence SSE ",
            join(["$(Int(h)) h: $(fmt(sse[h]))" for h in HALVING_H], ", "), ").\n")
    println(io, "| bars | n | prevalence RMSE model | prevalence RMSE persistence | bias model | Pearson model | Spearman model | Spearman intensity (model Mackin vs MD 0-7) |")
    println(io, "|---|---|---|---|---|---|---|---|")
    for role in ("calibration", "held out")
        r = res[res.role .== role, :]
        pr = r[.!ismissing.(r.prev_last), :]
        ri = r[.!ismissing.(r.int_obs), :]
        println(io, "| $role | $(nrow(r)) | $(fmt(rmse(r.prev_model, r.prev_obs))) | ",
                "$(fmt(rmse(Float64.(pr.prev_last), pr.prev_obs))) (n=$(nrow(pr))) | $(fmt(mean(r.prev_model .- r.prev_obs))) | ",
                "$(fmt(cor(r.prev_model, r.prev_obs))) | $(fmt(spearman(r.prev_model, r.prev_obs))) | ",
                "$(fmt(spearman(ri.mackin_model, Float64.(ri.int_obs)))) |")
    end
    println(io, "\nHeld-out bars with the papers' standard 60 h halving time: prevalence RMSE ",
            "$(fmt(rmse(std60.prev_model, std60.prev_obs))), Pearson $(fmt(cor(std60.prev_model, std60.prev_obs))).\n")
    println(io, "Per bar (chosen halving time):\n")
    println(io, "| bar | role | n | mean observed prevalence | mean model prevalence | RMSE model | RMSE persistence | Pearson |")
    println(io, "|---|---|---|---|---|---|---|---|")
    for g in groupby(res, :site_id)
        pr = g[.!ismissing.(g.prev_last), :]
        println(io, "| $(g.site_id[1]) | $(g.role[1]) | $(nrow(g)) | $(fmt(mean(g.prev_obs))) | $(fmt(mean(g.prev_model))) | ",
                "$(fmt(rmse(g.prev_model, g.prev_obs))) | $(nrow(pr) > 0 ? fmt(rmse(Float64.(pr.prev_last), pr.prev_obs)) : "n/a") | ",
                "$(nrow(g) > 2 ? fmt(cor(g.prev_model, g.prev_obs)) : "n/a") |")
    end
    missing_bars = setdiff(vcat(CALIB, HOLDOUT), unique(res.site_id))
    isempty(missing_bars) || println(io, "\nNo scorable years (forcing gaps): ", join(missing_bars, ", "), ".")
end
println(read(joinpath(OUT, "validation.md"), String))
