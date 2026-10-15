# P4.2f: universal differential equation version of the Dermo establishment model. A small neural network
# g_θ(S, T) >= 0 multiplies the waterborne transmission rate (simulate_establishment's transmission_mod), trained
# together with the mechanistic parameters on Maryland DNR and Great Bay, NH prevalence.
#
# Design (fixed before running; same data, objective and pass rule as P4.2e):
# - Data, forcing, initial states and survey handling: as scripts/p4_2e_establishment_fit.jl (helpers copied).
# - Network: inputs (S/30, (T - 15)/10), 2 -> 6 -> 6 -> 1, tanh hidden layers, softplus output; initialised at
#   g = 1 (last-layer weights 0, bias log(e - 1)). 67 weights. L2 penalty 1e-5 on all weights except the last bias
#   (a first run with 1e-3 collapsed the network to a constant: at the start the data gradient on the network is
#   ~1e-3 and the 1e-3 penalty gradient ~2e-3).
# - Mechanistic parameters fitted with it: log beta, log mu, logit host, log epsilon_MD, log epsilon_GB, starting from
#   the P4.2e full fit. The Great Bay arrival year is not differentiable and is fixed per training fold at the P4.2e
#   value for the same fold (fitted without the held-out data): full 2003.15; leave-out Adams 2004.1, Nannie 2003.2,
#   Oyster River 2002.2, Squamscott 2002.5, Woodman 2003.2; leave-out MD upper 2003.1, middle 2003.2, lower 2002.9.
# - Optimizer: LBFGS with forward-mode automatic differentiation (Optimization.jl, AutoForwardDiff), 150 iterations.
# - Tests: leave one Great Bay bed out (retrained from scratch per fold) and leave one Maryland region out. Pass rule:
#   beat persistence on held-out Great Bay beds, survey-weighted, on rows with a previous survey.
# - If it passes: Maine runs as P4.2e (Damariscotta, cycled 2016-2019, +0/0.5/1/2 degC, parasite present at 36 %
#   light prevalence), daily Dermo mortality files for farm coupling.
#
# Usage: julia --project scripts/p4_2f_ude_fit.jl
# Writes outputs/p4_2f/{ude_fit.md, ude_weights_full.csv, maine.arrow, daily_mortality_<site>_<dT>.arrow, progress.log}.

using OysterFarming
using Arrow, DataFrames, Dates, ForwardDiff, Random, Statistics
using Optimization, OptimizationOptimJL

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "outputs", "p4_2f")
const MD_START, MD_STOP = Date(1996, 10, 1), Date(2024, 11, 30)
const GB_START, GB_STOP = Date(1995, 1, 1), Date(2011, 12, 31)
const FIRST_MD_YEAR = 1999
const DENSITY0 = 20.0
const BACKGROUND_D = 0.0302 / 365
const WARMING = [0.0, 0.5, 1.0, 2.0]
const MAINE_PREV0 = 0.36
const MAINE_LEVEL0 = 10
const L2 = 1e-5
const MAXITERS = 150
const THETA_P42E = [-21.295, -0.69, 6.0, -8.478, -6.832]          # P4.2e full fit (log beta, log mu, logit host, log eps MD, log eps GB)
const ARRIVAL = Dict("full" => 2003.15, "nh_adams_point" => 2004.1, "nh_nannie_island" => 2003.2,
                     "nh_oyster_river" => 2002.2, "nh_squamscott_river" => 2002.5, "nh_woodman_point" => 2003.2,
                     "md_upper" => 2003.1, "md_middle" => 2003.2, "md_lower" => 2002.9)

mkpath(OUT)
const LOG = open(joinpath(OUT, "progress.log"), "w")
logln(x...) = (println(LOG, now(), " ", x...); flush(LOG); println(x...); flush(stdout))

lab = load_params(joinpath(ROOT, "params", "dermo_establishment_lab.toml"))
pw = load_params(joinpath(ROOT, "params", "dermo_powell1996.toml"))
forcing = read_contract(joinpath(ROOT, "data", "processed", "forcing_daily.arrow"), "forcing_daily")
disease = read_contract(joinpath(ROOT, "data", "processed", "disease_obs.arrow"), "disease_obs")
sites = read_contract(joinpath(ROOT, "data", "processed", "sites.arrow"), "sites")
lat = Dict(r.site_id => r.lat for r in eachrow(sites))
dermo = disease[disease.pathogen .== "P_marinus", :]

# ---- data helpers (as P4.2e) ----
function held(site, a, b)
    f = forcing_for_site(forcing, site)
    have = Dict(r.date => (Float64(r.temp_c), Float64(r.salinity_psu)) for r in eachrow(f))
    dates = collect(a:Day(1):b)
    last = have[minimum(keys(have))]
    T = Float64[]; S = Float64[]
    for d in dates
        haskey(have, d) && (last = have[d])
        push!(T, last[1]); push!(S, last[2])
    end
    return dates, T, S
end

struct SiteData
    id::String
    region::String
    T::Vector{Float64}
    S::Vector{Float64}
    rec::Vector{Int}
    obs::Vector{Float64}
    prev_obs::Vector{Union{Missing,Float64}}
end

function build(site, region, a, b, keep)
    _, T, S = held(site, a, b)
    o = sort(dermo[(dermo.site_id .== site) .& keep.(dermo.date), :], :date)
    o = o[(o.date .>= a) .& (o.date .<= b), :]
    nrow(o) == 0 && return nothing
    rec = [Dates.value(d - a) + 1 for d in o.date]
    prevs = Union{Missing,Float64}[i == 1 ? missing : Float64(o.prevalence[i-1]) for i in 1:nrow(o)]
    return SiteData(site, region, T, S, rec, Float64.(o.prevalence), prevs)
end

mdregion(s) = lat[s] >= 38.8 ? "md_upper" : lat[s] >= 38.45 ? "md_middle" : "md_lower"
fsites = Set(forcing.site_id)
data = SiteData[]
for s in sort(unique(dermo.site_id))
    s in fsites || continue
    x = startswith(s, "mddnr_") ? build(s, mdregion(s), MD_START, MD_STOP, d -> year(d) >= FIRST_MD_YEAR) :
        startswith(s, "nh_") ? build(s, "gb", GB_START, GB_STOP, d -> true) : nothing
    x === nothing || push!(data, x)
end
logln("sites: ", length(data), "; surveys: ", sum(length(d.obs) for d in data))

function initial(site::SiteData)
    N0 = zeros(EST_N_LEVELS)
    site.region == "gb" ? (N0[1] = DENSITY0) : (N0[1] = DENSITY0 / 2; N0[15] = DENSITY0 / 2)
    return N0
end

# ---- network ----
const NH = 6
const NW = NH * 2 + NH + NH * NH + NH + NH + 1             # 67
softplus(x) = log1p(exp(-abs(x))) + max(x, zero(x))
function g_net(w, T, S)
    x1 = S / 30; x2 = (T - 15) / 10
    o = 0
    W1 = reshape(view(w, o+1:o+2NH), NH, 2); o += 2NH
    b1 = view(w, o+1:o+NH); o += NH
    W2 = reshape(view(w, o+1:o+NH*NH), NH, NH); o += NH * NH
    b2 = view(w, o+1:o+NH); o += NH
    W3 = view(w, o+1:o+NH); o += NH
    b3 = w[o+1]
    h1 = tanh.(W1[:, 1] .* x1 .+ W1[:, 2] .* x2 .+ b1)
    h2 = tanh.(W2 * h1 .+ b2)
    return softplus(sum(W3 .* h2) + b3)
end
function init_weights(rng)
    w = 0.3 .* randn(rng, NW)
    w[end-NH:end-1] .= 0.0                  # W3 = 0
    w[end] = log(exp(1.0) - 1)              # softplus(b3) = 1
    return w
end
nn_weights(θ) = view(θ, 6:length(θ))

logistic(x) = 1 / (1 + exp(-x))
function params_from(θ, region)
    eps = region == "gb" ? exp(θ[5]) : exp(θ[4])
    return EstablishmentParams(lab, pw; beta=exp(θ[1]), mu=exp(θ[2]), epsilon=eps, host=logistic(θ[3]))
end
intro_day(arrival, site) = site.region == "gb" ? clamp(round(Int, (arrival - 1995) * 365.25) + 1, 1, length(site.T)) : 1

function predict(θ, site::SiteData, arrival)
    p = params_from(θ, site.region)
    w = nn_weights(θ)
    r = simulate_establishment(p, site.T, site.S, initial(site), 0.0; intro_day=intro_day(arrival, site),
                               background_mortality_d=BACKGROUND_D, record=site.rec,
                               transmission_mod=(t, s) -> g_net(w, t, s))
    return r.prev
end

function objective(θ, sites, arrival)
    md = zero(eltype(θ)); nmd = 0; gb = zero(eltype(θ)); ngb = 0
    for s in sites
        e = sum(abs2, predict(θ, s, arrival) .- s.obs)
        s.region == "gb" ? (gb += e; ngb += length(s.obs)) : (md += e; nmd += length(s.obs))
    end
    w = nn_weights(θ)
    reg = L2 * (sum(abs2, w) - w[end]^2)
    return (nmd > 0 ? md / nmd : 0.0) + (ngb > 0 ? gb / ngb : 0.0) + reg
end

function train(sites, arrival; seed=20261014)
    θ0 = vcat(THETA_P42E, init_weights(Xoshiro(seed)))
    f = OptimizationFunction((θ, _) -> objective(θ, sites, arrival), Optimization.AutoForwardDiff())
    t = @elapsed sol = solve(OptimizationProblem(f, θ0), LBFGS(); maxiters=MAXITERS)
    logln("  trained in ", round(t / 60; digits=1), " min, objective ", round(objective(sol.u, sites, arrival); digits=4))
    return sol.u
end

rmse(a, b) = sqrt(mean((a .- b) .^ 2))
function score(θ, test, arrival)
    pred = reduce(vcat, [predict(θ, s, arrival) for s in test]); obs = reduce(vcat, [s.obs for s in test])
    pers = reduce(vcat, [s.prev_obs for s in test]); ok = .!ismissing.(pers)
    return (n=length(obs), rmse=rmse(pred, obs), cor=length(obs) > 2 ? cor(pred, obs) : NaN,
            rmse_pers=rmse(Float64.(pers[ok]), obs[ok]), rmse_pers_model=rmse(pred[ok], obs[ok]), n_pers=count(ok),
            rmse_mean=rmse(fill(mean(obs), length(obs)), obs), bias=mean(pred .- obs))
end

# timing check
tg = @elapsed ForwardDiff.gradient(θ -> objective(θ, data, ARRIVAL["full"]), vcat(THETA_P42E, init_weights(Xoshiro(1))))
logln("one gradient: ", round(tg; digits=1), " s (includes compilation)")

logln("training on all data ...")
θ̂ = train(data, ARRIVAL["full"])
open(joinpath(OUT, "ude_weights_full.csv"), "w") do io
    println(io, "index,value"); for (i, v) in enumerate(θ̂); println(io, "$i,$v"); end
end

loo = DataFrame(heldout=String[], n=Int[], rmse=Float64[], cor=Float64[], rmse_model_pers_rows=Float64[],
                rmse_persistence=Float64[], n_pers=Int[], rmse_mean=Float64[], bias=Float64[])
for s in [x for x in data if x.region == "gb"]
    logln("leave out ", s.id)
    θf = train([x for x in data if x.id != s.id], ARRIVAL[s.id])
    sc = score(θf, [s], ARRIVAL[s.id])
    push!(loo, (s.id, sc.n, sc.rmse, sc.cor, sc.rmse_pers_model, sc.rmse_pers, sc.n_pers, sc.rmse_mean, sc.bias))
    logln("  held-out RMSE ", round(sc.rmse_pers_model; digits=3), " vs persistence ", round(sc.rmse_pers; digits=3))
end
lro = similar(loo, 0)
for reg in ("md_upper", "md_middle", "md_lower")
    logln("leave out ", reg)
    θf = train([x for x in data if x.region != reg], ARRIVAL[reg])
    sc = score(θf, [x for x in data if x.region == reg], ARRIVAL[reg])
    push!(lro, (reg, sc.n, sc.rmse, sc.cor, sc.rmse_pers_model, sc.rmse_pers, sc.n_pers, sc.rmse_mean, sc.bias))
    logln("  held-out RMSE ", round(sc.rmse_pers_model; digits=3), " vs persistence ", round(sc.rmse_pers; digits=3))
end
insample = Dict(reg => score(θ̂, [s for s in data if s.region == reg], ARRIVAL["full"]) for reg in ("gb", "md_upper", "md_middle", "md_lower"))
wavg(x, w) = sum(x .* w) / sum(w)
gb_model = wavg(loo.rmse_model_pers_rows, loo.n_pers); gb_pers = wavg(loo.rmse_persistence, loo.n_pers)
gb_pass = gb_model < gb_pers

maine = DataFrame(site=String[], warming_c=Float64[], year=Int[], fall_prev=Float64[], fall_true_prev=Float64[],
                  fall_mackin=Float64[], dermo_mortality_yr=Float64[])
full_years(b) = [y for y in unique(year.(b.date)) if count(year.(b.date) .== y) == daysinyear(y)]
if gb_pass
    pM = params_from(θ̂, "gb"); wM = nn_weights(θ̂)
    for site in ("dam_lower", "dam_upper")
        b = longest_block(forcing_for_site(forcing, site))
        f = cycle_years(b, full_years(b), 10)
        dates = f.date; yrs = unique(year.(dates))
        for dT in WARMING
            T = Float64.(f.temp_c) .+ dT; S = Float64.(f.salinity_psu)
            N0 = zeros(EST_N_LEVELS); N0[1] = DENSITY0 * (1 - MAINE_PREV0); N0[MAINE_LEVEL0] = DENSITY0 * MAINE_PREV0
            rec = [findfirst(==(Date(y, 11, 1)), dates) for y in yrs]
            r = simulate_establishment(pM, T, S, N0, 0.0; background_mortality_d=BACKGROUND_D, record=rec,
                                       transmission_mod=(t, s) -> g_net(wM, t, s))
            for (i, y) in enumerate(yrs)
                push!(maine, (site, dT, y, r.prev[i], r.true_prev[i], r.mackin[i], sum(r.deaths[year.(dates) .== y]) / DENSITY0))
            end
            Arrow.write(joinpath(OUT, "daily_mortality_$(site)_$(dT).arrow"), DataFrame(date=dates, dermo_mortality_d=r.deaths ./ DENSITY0))
        end
    end
end
Arrow.write(joinpath(OUT, "maine.arrow"), maine)

fmt(x; d=3) = (x isa Real && isnan(x)) ? "n/a" : string(round(x; digits=d))
open(joinpath(OUT, "ude_fit.md"), "w") do io
    println(io, "# P4.2f UDE Dermo model: learned transmission multiplier g(S, T) (scripts/p4_2f_ude_fit.jl)\n")
    println(io, "$(length(data)) sites, $(sum(length(d.obs) for d in data)) surveys; network $(NW) weights, L2 $(L2), LBFGS $(MAXITERS) iterations.\n")
    vals = [exp(θ̂[1]), exp(θ̂[2]), logistic(θ̂[3]), exp(θ̂[4]), exp(θ̂[5])]
    println(io, "| parameter | value |\n|---|---|")
    for (n, v) in zip(["beta", "mu (d^-1)", "host factor", "epsilon MD (d^-1)", "epsilon GB (d^-1)"], vals)
        println(io, "| $n | $(round(v; sigdigits=3)) |")
    end
    println(io, "\nLearned multiplier g(S, T) (full fit):\n")
    Ts = [5.0, 15.0, 22.0, 28.0]
    println(io, "| S (ppt) | " * join(["T = $(Int(t)) degC" for t in Ts], " | ") * " |")
    println(io, "|---|" * repeat("---|", length(Ts)))
    for S in (2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0)
        println(io, "| $(Int(S)) | " * join([fmt(g_net(nn_weights(θ̂), t, S)) for t in Ts], " | ") * " |")
    end
    println(io, "\n## In-sample\n\n| region | n | RMSE | Pearson | RMSE persistence | RMSE mean | bias |\n|---|---|---|---|---|---|---|")
    for (reg, sc) in sort(collect(insample))
        println(io, "| $reg | $(sc.n) | $(fmt(sc.rmse)) | $(fmt(sc.cor)) | $(fmt(sc.rmse_pers)) | $(fmt(sc.rmse_mean)) | $(fmt(sc.bias)) |")
    end
    for (title, t) in (("Leave one Great Bay bed out", loo), ("Leave one Maryland region out", lro))
        println(io, "\n## $title\n\n| held out | n | RMSE | Pearson | RMSE model (rows with previous survey) | RMSE persistence (same rows) | RMSE held-out mean | bias |")
        println(io, "|---|---|---|---|---|---|---|---|")
        for r in eachrow(t)
            println(io, "| $(r.heldout) | $(r.n) | $(fmt(r.rmse)) | $(fmt(r.cor)) | $(fmt(r.rmse_model_pers_rows)) | $(fmt(r.rmse_persistence)) | $(fmt(r.rmse_mean)) | $(fmt(r.bias)) |")
        end
    end
    println(io, "\nGreat Bay held-out, survey-weighted: model $(fmt(gb_model)) vs persistence $(fmt(gb_pers)). Pass: **$(gb_pass)**.\n")
    if gb_pass
        println(io, "## Maine (Damariscotta; parasite present at the start; cycled 2016-2019 forcing)\n")
        println(io, "| site | warming degC | year | fall RFTM prevalence | fall true prevalence | Mackin | Dermo mortality yr^-1 |\n|---|---|---|---|---|---|---|")
        for r in eachrow(maine)
            println(io, "| $(r.site) | $(r.warming_c) | $(r.year) | $(fmt(r.fall_prev)) | $(fmt(r.fall_true_prev)) | $(fmt(r.fall_mackin; d=2)) | $(fmt(r.dermo_mortality_yr; d=4)) |")
        end
    end
end
logln(read(joinpath(OUT, "ude_fit.md"), String))
close(LOG)
