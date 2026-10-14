# P4.2e: fit and test the Dermo establishment model (src/dermo_establishment.jl) on Maryland and Great Bay, NH.
#
# Design (fixed before running):
# - Data: fall RFTM prevalence. Maryland: 42 MD DNR bars, surveys dated 1 Nov 1999-2024. Great Bay: 5 beds, all
#   NH Fish & Game surveys 1996-2011 (sonde forcing ends 2011). Forcing: forcing_daily; days without a row repeat
#   the last value (Great Bay winters have no sonde data; below 6 degC the model's lab rate is constant).
# - Initial state: 20 market oysters m^-2 (Powell et al. 1996 p. 143). Maryland runs start 1996-10-01 with 50 %
#   prevalence (same source) at level 15 (just above RFTM detection); 1996-1998 is spin-up. Great Bay runs start
#   1995-01-01 uninfected (Dermo was first reported in Great Bay in fall 1995; NH F&G 2015 report p. 4); the
#   background import epsilon_GB starts at the fitted arrival day.
# - Fitted (Nelder-Mead on log/logit scales, 4 starts): beta, mu, host (shared), epsilon_MD, epsilon_GB, arrival
#   day for Great Bay (shared by its beds). Objective: mean squared prevalence error in Maryland + mean in Great Bay
#   (each region weighted equally).
# - Tests: (1) leave one Great Bay bed out (refit on the rest, predict the held-out bed); (2) leave one Maryland
#   region out (upper / middle / lower bay by latitude, as P4.2c). Baselines: persistence (last survey at the
#   same site) and the mean of the held-out data (generous: it uses the test data).
# - Maine (only reported if the model beats persistence on held-out Great Bay beds): Damariscotta forcing
#   (complete years 2016-2019 cycled to 10 years), +0/0.5/1/2 degC, parasite present at the start with true
#   prevalence 0.36 at light intensity (level 10, below RFTM detection; Marquis et al. 2020 qPCR, Damariscotta
#   sites 34-39 % in 2016-17), epsilon = fitted epsilon_GB. Outputs: fall RFTM prevalence, Mackin intensity,
#   annual Dermo mortality.
#
# Usage: julia --project scripts/p4_2e_establishment_fit.jl
# Writes outputs/p4_2e/{fit.md, maine.arrow, daily_mortality_<site>_<dT>.arrow} and params/dermo_establishment_fit.toml.

using OysterFarming
using Arrow, DataFrames, Dates, Random, Statistics
using Optimization, OptimizationOptimJL

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "outputs", "p4_2e")
const MD_START, MD_STOP = Date(1996, 10, 1), Date(2024, 11, 30)
const GB_START, GB_STOP = Date(1995, 1, 1), Date(2011, 12, 31)
const FIRST_MD_YEAR = 1999
const DENSITY0 = 20.0
const BACKGROUND_D = 0.0302 / 365
const WARMING = [0.0, 0.5, 1.0, 2.0]
const MAINE_PREV0 = 0.36
const MAINE_LEVEL0 = 10

lab = load_params(joinpath(ROOT, "params", "dermo_establishment_lab.toml"))
pw = load_params(joinpath(ROOT, "params", "dermo_powell1996.toml"))
forcing = read_contract(joinpath(ROOT, "data", "processed", "forcing_daily.arrow"), "forcing_daily")
disease = read_contract(joinpath(ROOT, "data", "processed", "disease_obs.arrow"), "disease_obs")
sites = read_contract(joinpath(ROOT, "data", "processed", "sites.arrow"), "sites")
lat = Dict(r.site_id => r.lat for r in eachrow(sites))
dermo = disease[disease.pathogen .== "P_marinus", :]

"Daily T, S over [a, b], repeating the last available value on days without a forcing row."
function held(site, a, b)
    f = forcing_for_site(forcing, site)
    have = Dict(r.date => (Float64(r.temp_c), Float64(r.salinity_psu)) for r in eachrow(f))
    dates = collect(a:Day(1):b)
    k0 = minimum(keys(have))
    last = have[k0]
    T = Float64[]; S = Float64[]
    for d in dates
        haskey(have, d) && (last = have[d])
        push!(T, last[1]); push!(S, last[2])
    end
    return dates, T, S
end

struct SiteData
    id::String
    region::String          # "md_upper", "md_middle", "md_lower", "gb"
    T::Vector{Float64}
    S::Vector{Float64}
    start::Date
    rec::Vector{Int}        # day indices of surveys (sorted)
    obs::Vector{Float64}
    prev_obs::Vector{Union{Missing,Float64}}   # last survey at the site before this one (persistence)
end

function build(site, region, a, b, keep)
    dates, T, S = held(site, a, b)
    o = sort(dermo[(dermo.site_id .== site) .& keep.(dermo.date), :], :date)
    o = o[(o.date .>= a) .& (o.date .<= b), :]
    nrow(o) == 0 && return nothing
    rec = [Dates.value(d - a) + 1 for d in o.date]
    prevs = Union{Missing,Float64}[i == 1 ? missing : Float64(o.prevalence[i-1]) for i in 1:nrow(o)]
    return SiteData(site, region, T, S, a, rec, Float64.(o.prevalence), prevs)
end

mdregion(s) = lat[s] >= 38.8 ? "md_upper" : lat[s] >= 38.45 ? "md_middle" : "md_lower"
fsites = Set(forcing.site_id)
data = SiteData[]
for s in sort(unique(dermo.site_id))
    s in fsites || continue
    if startswith(s, "mddnr_")
        x = build(s, mdregion(s), MD_START, MD_STOP, d -> year(d) >= FIRST_MD_YEAR)
    elseif startswith(s, "nh_")
        x = build(s, "gb", GB_START, GB_STOP, d -> true)
    else
        continue
    end
    x === nothing || push!(data, x)
end
println("sites: ", length(data), " (", count(d -> d.region == "gb", data), " Great Bay); surveys: ", sum(length(d.obs) for d in data))

function initial(site::SiteData)
    N0 = zeros(EST_N_LEVELS)
    if site.region == "gb"
        N0[1] = DENSITY0
    else
        N0[1] = DENSITY0 * 0.5; N0[15] = DENSITY0 * 0.5
    end
    return N0
end

# θ = [log β, log μ, logit host, log ε_MD, log ε_GB, arrival year offset (years after 1995-01-01)]
logistic(x) = 1 / (1 + exp(-x))
function params_from(θ, region)
    eps = region == "gb" ? exp(θ[5]) : exp(θ[4])
    return EstablishmentParams(lab, pw; beta=exp(θ[1]), mu=exp(θ[2]), epsilon=eps, host=logistic(θ[3]))
end
intro_day(θ, site) = site.region == "gb" ? clamp(round(Int, θ[6] * 365.25) + 1, 1, length(site.T)) : 1

function predict(θ, site::SiteData)
    p = params_from(θ, site.region)
    r = simulate_establishment(p, site.T, site.S, initial(site), 0.0; intro_day=intro_day(θ, site),
                               background_mortality_d=BACKGROUND_D, record=site.rec)
    return r.prev
end

function objective(θ, sites)
    md = 0.0; nmd = 0; gb = 0.0; ngb = 0
    for s in sites
        e = sum(abs2, predict(θ, s) .- s.obs)
        if s.region == "gb"
            gb += e; ngb += length(s.obs)
        else
            md += e; nmd += length(s.obs)
        end
    end
    return (nmd > 0 ? md / nmd : 0.0) + (ngb > 0 ? gb / ngb : 0.0)
end

const LB = [log(1e-14), log(1e-3), -6.0, log(1e-7), log(1e-8), 0.0]
const UB = [log(1e-3), log(1.0), 6.0, log(1e-1), log(1e-1), 12.0]
function fit(sites; starts=4, seed=20261014, θinit=nothing)
    rng = Xoshiro(seed)
    best = (θ=Float64[], f=Inf)
    inits = θinit === nothing ? [LB .+ rand(rng, 6) .* (UB .- LB) for _ in 1:starts] : [θinit]
    for θ0 in inits
        prob = OptimizationProblem((θ, _) -> objective(θ, sites), θ0; lb=LB, ub=UB)
        sol = solve(prob, NelderMead(); maxiters=1500)
        f = objective(sol.u, sites)
        f < best.f && (best = (θ=copy(sol.u), f=f))
    end
    return best
end

println("fitting all data ...")
full = fit(data; starts=6)
θ̂ = full.θ
println("fit: ", round.(θ̂; digits=3), " objective ", round(full.f; digits=4))

rmse(a, b) = sqrt(mean((a .- b) .^ 2))
function score(θ, test)
    pred = reduce(vcat, [predict(θ, s) for s in test]); obs = reduce(vcat, [s.obs for s in test])
    pers = reduce(vcat, [s.prev_obs for s in test])
    ok = .!ismissing.(pers)
    return (n=length(obs), rmse=rmse(pred, obs), cor=length(obs) > 2 ? cor(pred, obs) : NaN,
            rmse_pers=rmse(Float64.(pers[ok]), obs[ok]), rmse_pers_model=rmse(pred[ok], obs[ok]), n_pers=count(ok),
            rmse_mean=rmse(fill(mean(obs), length(obs)), obs), bias=mean(pred .- obs))
end

println("leave one Great Bay bed out ...")
loo = DataFrame(heldout=String[], n=Int[], rmse=Float64[], cor=Float64[], rmse_model_pers_rows=Float64[],
                rmse_persistence=Float64[], n_pers=Int[], rmse_mean=Float64[], bias=Float64[], arrival_year=Float64[])
gbsites = [s for s in data if s.region == "gb"]
for s in gbsites
    tr = [x for x in data if x.id != s.id]
    f = fit(tr; θinit=θ̂)
    sc = score(f.θ, [s])
    push!(loo, (s.id, sc.n, sc.rmse, sc.cor, sc.rmse_pers_model, sc.rmse_pers, sc.n_pers, sc.rmse_mean, sc.bias, 1995 + f.θ[6]))
end
println("leave one Maryland region out ...")
lro = similar(loo, 0)
for reg in ("md_upper", "md_middle", "md_lower")
    tr = [x for x in data if x.region != reg]; te = [x for x in data if x.region == reg]
    f = fit(tr; θinit=θ̂)
    sc = score(f.θ, te)
    push!(lro, (reg, sc.n, sc.rmse, sc.cor, sc.rmse_pers_model, sc.rmse_pers, sc.n_pers, sc.rmse_mean, sc.bias, 1995 + f.θ[6]))
end
insample = Dict(reg => score(θ̂, [s for s in data if s.region == reg]) for reg in ("gb", "md_upper", "md_middle", "md_lower"))

gb_pass = sum(loo.rmse_model_pers_rows .* loo.n_pers) / sum(loo.n_pers) < sum(loo.rmse_persistence .* loo.n_pers) / sum(loo.n_pers)

# Maine
maine = DataFrame(site=String[], warming_c=Float64[], year=Int[], fall_prev=Float64[], fall_true_prev=Float64[],
                  fall_mackin=Float64[], dermo_mortality_yr=Float64[])
full_years(b) = [y for y in unique(year.(b.date)) if count(year.(b.date) .== y) == daysinyear(y)]
if gb_pass
    pM = params_from(θ̂, "gb")
    for site in ("dam_lower", "dam_upper")
        b = longest_block(forcing_for_site(forcing, site))
        f = cycle_years(b, full_years(b), 10)
        dates = f.date
        for dT in WARMING
            T = Float64.(f.temp_c) .+ dT; S = Float64.(f.salinity_psu)
            N0 = zeros(EST_N_LEVELS); N0[1] = DENSITY0 * (1 - MAINE_PREV0); N0[MAINE_LEVEL0] = DENSITY0 * MAINE_PREV0
            rec = [findfirst(==(Date(y, 11, 1)), dates) for y in unique(year.(dates))]
            r = simulate_establishment(pM, T, S, N0, 0.0; background_mortality_d=BACKGROUND_D, record=rec)
            for (i, y) in enumerate(unique(year.(dates)))
                dy = r.deaths[year.(dates) .== y]
                push!(maine, (site, dT, y, r.prev[i], r.true_prev[i], r.mackin[i], sum(dy) / DENSITY0))
            end
            mkpath(OUT)
            Arrow.write(joinpath(OUT, "daily_mortality_$(site)_$(dT).arrow"),
                        DataFrame(date=dates, dermo_mortality_d=r.deaths ./ DENSITY0))
        end
    end
end

mkpath(OUT)
Arrow.write(joinpath(OUT, "maine.arrow"), maine)
fmt(x; d=3) = (x isa Real && isnan(x)) ? "n/a" : string(round(x; digits=d))
names = ["beta (m^2 cell^-1 d^-1)", "mu (d^-1)", "host factor", "epsilon Maryland (d^-1)", "epsilon Great Bay (d^-1)", "Great Bay arrival"]
vals = [exp(θ̂[1]), exp(θ̂[2]), logistic(θ̂[3]), exp(θ̂[4]), exp(θ̂[5]), 1995 + θ̂[6]]
open(joinpath(OUT, "fit.md"), "w") do io
    println(io, "# P4.2e Dermo establishment model (generated by scripts/p4_2e_establishment_fit.jl)\n")
    println(io, "$(length(data)) sites ($(length(gbsites)) Great Bay), $(sum(length(d.obs) for d in data)) surveys. Objective (mean squared ",
            "error Maryland + Great Bay): $(fmt(full.f; d=4)).\n")
    println(io, "| parameter | fitted value |\n|---|---|")
    for (n, v) in zip(names, vals)
        println(io, "| $n | $(round(v; sigdigits=3)) |")
    end
    println(io, "\n## In-sample fit by region\n\n| region | n | RMSE | Pearson | RMSE persistence | RMSE region mean | bias |\n|---|---|---|---|---|---|---|")
    for (reg, sc) in sort(collect(insample))
        println(io, "| $reg | $(sc.n) | $(fmt(sc.rmse)) | $(fmt(sc.cor)) | $(fmt(sc.rmse_pers)) | $(fmt(sc.rmse_mean)) | $(fmt(sc.bias)) |")
    end
    for (title, t) in (("Leave one Great Bay bed out", loo), ("Leave one Maryland region out", lro))
        println(io, "\n## $title\n\n| held out | n | RMSE | Pearson | RMSE model on rows with a previous survey | RMSE persistence (same rows) | RMSE held-out mean | bias | fitted Great Bay arrival |")
        println(io, "|---|---|---|---|---|---|---|---|---|")
        for r in eachrow(t)
            println(io, "| $(r.heldout) | $(r.n) | $(fmt(r.rmse)) | $(fmt(r.cor)) | $(fmt(r.rmse_model_pers_rows)) | $(fmt(r.rmse_persistence)) | ",
                    "$(fmt(r.rmse_mean)) | $(fmt(r.bias)) | $(fmt(r.arrival_year; d=1)) |")
        end
    end
    println(io, "\nGreat Bay held-out test passed (model beats persistence, survey-weighted): **$(gb_pass)**.\n")
    if gb_pass
        println(io, "## Maine (Damariscotta, parasite present at the start; years are cycled 2016-2019 forcing)\n")
        println(io, "| site | warming degC | year | fall RFTM prevalence | fall true prevalence | Mackin intensity | Dermo mortality (fraction yr^-1) |")
        println(io, "|---|---|---|---|---|---|---|")
        for r in eachrow(maine)
            println(io, "| $(r.site) | $(r.warming_c) | $(r.year) | $(fmt(r.fall_prev)) | $(fmt(r.fall_true_prev)) | $(fmt(r.fall_mackin; d=2)) | $(fmt(r.dermo_mortality_yr; d=4)) |")
        end
    else
        println(io, "Maine runs not made: the model did not beat persistence on held-out Great Bay beds.")
    end
end

if gb_pass
    open(joinpath(ROOT, "params", "dermo_establishment_fit.toml"), "w") do io
        println(io, "title = \"Fitted parameters of the Dermo establishment model (PLAN.md P4.2e)\"\n")
        src = "Fitted in this project: scripts/p4_2e_establishment_fit.jl (Nelder-Mead, Maryland DNR bars 1999-2024 and NH Fish & Game Great Bay beds 1996-2011, RFTM prevalence); results and held-out tests in outputs/p4_2e/fit.md"
        for (key, v, u) in (("beta", vals[1], "m^2 cell^-1 d^-1"), ("mu", vals[2], "d^-1"), ("host", vals[3], "1"),
                            ("epsilon_md", vals[4], "d^-1"), ("epsilon_gb", vals[5], "d^-1"), ("gb_arrival_year", vals[6], "yr"))
            println(io, "[$key]\nvalue = $(round(v; sigdigits=5))\nunits = \"$u\"\nsource = \"$src\"\n")
        end
    end
end
println(read(joinpath(OUT, "fit.md"), String))
