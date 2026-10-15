# P4.2g: nested cross-validation of the Dermo establishment model with learned, shape-constrained transmission curves,
# trained on every system (Maryland DNR bars, Great Bay NH beds, Ford & Smolowitz 2007 New England sites, Maine qPCR).
#
# Design (fixed before running):
# - Sites and data:
#   Maryland: 40+ MD DNR bars, fall RFTM prevalence 1999-2024 (intensity scale 0-7 not used), forcing_daily from
#     1996-10-01 (spin-up), start 50 % infected at level 15 (Powell et al. 1996).
#   Great Bay: 5 NH F&G beds, fall RFTM prevalence and Mackin 0-5 weighted prevalence (WP) 1996-2011, NERRS forcing from
#     1995-01-01 (winter gaps repeat the last value), start uninfected.
#   New England: 8 Ford & Smolowitz (2007) sites, monthly RFTM prevalence and WP (data/literature/dermo/
#     ford2007_smolowitz_fig2_fig3_digitized.csv), daily temperature interpolated from the digitized Fig. 2 traces
#     (held after the last point), salinity constant at the midpoint of Table 3 ranges (Long Island Sound sites:
#     midpoint 27.5 of "20-35 ppt, mostly the high end", text p125; Edgartown: midpoint 15 of 6-24 ppt).
#     Start at the first sample: observed prevalence placed at level 15 (just above RFTM detection); that sample is
#     not scored.
#   Maine: Marquis et al. (2020) qPCR prevalence at the Damariscotta sites (data/literature/dermo/
#     marquis2020_maine_qpcr_monthly.csv, if present), forcing dam_lower / dam_upper; start at the first sample with the
#     observed prevalence at level 10 (light, below RFTM detection; qPCR detects it); compared with model true prevalence.
# - Model: simulate_establishment (smooth mode) with transmission multiplier g(T, S) = gS(S) * gT(T); gS and gT are
#   1 -> H -> 1 tanh networks with softplus output, initialised at 1. Penalties: monotone non-decreasing in S and in T
#   (squared negative steps on 1-unit grids, weight 1; Chu & Greene 1989: infective stages die within a day at 4 degC),
#   L2 on weights. Shared: log beta, log mu, logit host, gS, gT. Per system: log epsilon (background import). A held-out
#   site uses its system's epsilon fitted on the other sites of that system.
# - Loss: for each system in the training set, mean squared prevalence error + 0.5 * mean squared error of WP/5 vs model
#   Mackin/5 (where WP is observed); systems weighted equally.
# - Nested CV: outer folds hold out one site group (sites sharing a forcing station: Maryland by CBP station, Great Bay
#   by sonde, Long Island Sound Westport + Milford together, each other site alone). Inner: 3 folds of the remaining
#   site groups (round-robin within each system), hyperparameters H in {3, 6} x L2 in {1e-5, 1e-3}; the configuration
#   with the lowest mean inner validation loss is refit on all outer-training sites and scored once on the test group.
# - Training: LBFGS, AutoForwardDiff, 100 iterations, fixed initialisation. Folds run in parallel threads.
# - Scores per held-out site (scored rows only): prevalence RMSE vs persistence (previous observation at the same site)
#   and vs the training mean of the site's system. Pass rule per system: model RMSE below persistence, with the 95 %
#   site-bootstrap interval of (model - persistence) below 0.
#
# Usage: julia --project -t 12 scripts/p4_2g_nested.jl
# Writes outputs/p4_2g/{nested.md, predictions.arrow, folds.arrow, progress.log}.

using OysterFarming
using Arrow, DataFrames, Dates, ForwardDiff, LinearAlgebra, Random, Statistics
using Optimization, OptimizationOptimJL

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "outputs", get(ENV, "P42G_DRY", "0") == "1" ? "p4_2g_dry" : "p4_2g")
const LIT = joinpath(ROOT, "data", "literature", "dermo")
const DENSITY0 = 20.0
const BACKGROUND_D = 0.0302 / 365
const DRY = get(ENV, "P42G_DRY", "0") == "1"     # smoke test: 2 outer folds, 1 configuration, 3 iterations
const MAXITERS = DRY ? 3 : 100
const CONFIGS = DRY ? [(H=3, L2=1e-5)] : [(H=3, L2=1e-5), (H=3, L2=1e-3), (H=6, L2=1e-5), (H=6, L2=1e-3)]
const N_INNER = 3
const MONO_W = 1.0
const SYSTEMS = ["md", "gb", "ne", "me"]
mkpath(OUT)
const LOG = open(joinpath(OUT, "progress.log"), "a")
const LOCK = ReentrantLock()
logln(x...) = lock(LOCK) do
    println(LOG, now(), " ", x...); flush(LOG); println(x...); flush(stdout)
end

lab = load_params(joinpath(ROOT, "params", "dermo_establishment_lab.toml"))
pw = load_params(joinpath(ROOT, "params", "dermo_powell1996.toml"))
forcing = read_contract(joinpath(ROOT, "data", "processed", "forcing_daily.arrow"), "forcing_daily")
disease = read_contract(joinpath(ROOT, "data", "processed", "disease_obs.arrow"), "disease_obs")
sites_tbl = read_contract(joinpath(ROOT, "data", "processed", "sites.arrow"), "sites")
dermo = disease[disease.pathogen .== "P_marinus", :]

struct Site
    id::String
    system::String
    group::String
    T::Vector{Float64}
    S::Vector{Float64}
    N0::Vector{Float64}
    prev_rec::Vector{Int}
    prev_obs::Vector{Float64}
    pcr::Bool                       # true: compare with true prevalence (qPCR); false: RFTM-detectable prevalence
    wp_rec::Vector{Int}
    wp_obs::Vector{Float64}
    prev_persist::Vector{Union{Missing,Float64}}
end

"Daily T, S over [a, b] from forcing_daily, repeating the last value on days without a row."
function held(site, a, b)
    f = forcing_for_site(forcing, site)
    have = Dict(r.date => (Float64(r.temp_c), Float64(r.salinity_psu)) for r in eachrow(f))
    last = have[minimum(keys(have))]
    T = Float64[]; S = Float64[]
    for d in a:Day(1):b
        haskey(have, d) && (last = have[d])
        push!(T, last[1]); push!(S, last[2])
    end
    return T, S
end
dayidx(d, a) = Dates.value(d - a) + 1
persist(v) = Union{Missing,Float64}[i == 1 ? missing : v[i-1] for i in eachindex(v)]
function start_state(p0, level)
    N0 = zeros(EST_N_LEVELS); N0[1] = DENSITY0 * (1 - p0); N0[level] += DENSITY0 * p0
    return N0
end

sites = Site[]
fsites = Set(forcing.site_id)
src_of = Dict(r.site_id => r.source for r in eachrow(forcing) if startswith(r.site_id, "mddnr_"))
notes_of = Dict(r.site_id => coalesce(r.notes, "") for r in eachrow(sites_tbl))

# Maryland
let a = Date(1996, 10, 1), b = Date(2024, 11, 30)
    for s in sort(unique(dermo.site_id))
        (startswith(s, "mddnr_") && s in fsites) || continue
        o = sort(dermo[(dermo.site_id .== s) .& (year.(dermo.date) .>= 1999), :], :date)
        nrow(o) < 3 && continue
        T, S = held(s, a, b)
        m = match(r"cbp:([A-Z0-9.]+)", src_of[s])
        grp = "md_" * (m === nothing ? s : m[1])
        push!(sites, Site(s, "md", grp, T, S, start_state(0.5, 15), dayidx.(o.date, a), Float64.(o.prevalence), false,
                          Int[], Float64[], persist(Float64.(o.prevalence))))
    end
end
# Great Bay
let a = Date(1995, 1, 1), b = Date(2011, 12, 31)
    for s in sort(unique(dermo.site_id))
        (startswith(s, "nh_") && s in fsites) || continue
        o = sort(dermo[(dermo.site_id .== s) .& (dermo.date .<= b), :], :date)
        T, S = held(s, a, b)
        m = match(r"sonde (grb[a-z]+)", notes_of[s])
        push!(sites, Site(s, "gb", "gb_" * (m === nothing ? s : m[1]), T, S, start_state(0.0, 15), dayidx.(o.date, a),
                          Float64.(o.prevalence), false, dayidx.(o.date, a), Float64.(o.intensity_mean),
                          persist(Float64.(o.prevalence))))
    end
end
# New England (Ford & Smolowitz 2007)
"Split one CSV line, honouring double-quoted fields that contain commas."
function csvsplit(l)
    out = String[]; buf = IOBuffer(); q = false
    for c in l
        if c == '"'
            q = !q
        elseif c == ',' && !q
            push!(out, String(take!(buf)))
        else
            write(buf, c)
        end
    end
    push!(out, String(take!(buf)))
    return out
end
function read_csv_simple(f)
    lines = [l for l in readlines(joinpath(LIT, f)) if !startswith(l, "#") && !isempty(strip(l))]
    hdr = csvsplit(lines[1]); rows = [csvsplit(l) for l in lines[2:end]]
    return DataFrame([Symbol(h) => [String(r[i]) for r in rows] for (i, h) in enumerate(hdr)])
end
fs = read_csv_simple("ford2007_smolowitz_fig2_fig3_digitized.csv")
tr = read_csv_simple("ford2007_smolowitz_fig2_temperature_trace.csv")
const NE_SAL = Dict("Wellfleet MA" => 31.0, "Cotuit MA" => 26.0, "Oyster Bay NY" => 24.0, "Great Sound NJ" => 29.5,
                    "Delaware Bay NJ" => 20.0, "Martha's Vineyard MA (Butler's Neck)" => 15.0,
                    "Westport CT" => 27.5, "Milford CT" => 27.5)
for site in sort(unique(fs.site))
    haskey(NE_SAL, site) || continue
    t = tr[tr.site .== site, :]
    nrow(t) == 0 && continue
    tdates = Date.(t.date); tvals = parse.(Float64, t.temp_C)
    pv = fs[(fs.site .== site) .& (fs.variable .== "prevalence_pct"), :]
    wv = fs[(fs.site .== site) .& (fs.variable .== "weighted_prevalence_0to5"), :]
    pdates = Date.(pv.date); pvals = parse.(Float64, pv.value) ./ 100
    o = sortperm(pdates); pdates = pdates[o]; pvals = clamp.(pvals[o], 0.0, 1.0)
    a = pdates[1]; b = max(maximum(pdates), maximum(Date.(wv.date)))
    days = collect(a:Day(1):b)
    # daily temperature: mean of trace vertices per date, linear interpolation between dates, held at the ends
    ud = sort(unique(tdates)); uv = [mean(tvals[tdates .== d]) for d in ud]
    T = [d <= ud[1] ? uv[1] : d >= ud[end] ? uv[end] :
         (i = searchsortedlast(ud, d); uv[i] + (uv[i+1] - uv[i]) * Dates.value(d - ud[i]) / Dates.value(ud[i+1] - ud[i])) for d in days]
    S = fill(NE_SAL[site], length(days))
    wd = Date.(wv.date); wvals = parse.(Float64, wv.value); keepw = wd .> a
    grp = site in ("Westport CT", "Milford CT") ? "ne_LIS" : "ne_" * site
    sid = "ne_" * lowercase(replace(site, r"[^A-Za-z]+" => "_"))
    push!(sites, Site(sid, "ne", grp, T, S, start_state(pvals[1], 15), dayidx.(pdates[2:end], a), pvals[2:end], false,
                      dayidx.(wd[keepw], a), wvals[keepw], persist(pvals)[2:end]))
end
# Maine (Marquis et al. 2020), if digitized
const MAINE_FILE = "marquis2020_maine_qpcr_monthly.csv"
if isfile(joinpath(LIT, MAINE_FILE))
    me = read_csv_simple(MAINE_FILE)
    for (site, model_site) in (("Jack's Point", "dam_upper"), ("Prentiss Island", "dam_lower"))
        r = me[occursin.(lowercase(site), lowercase.(me.site)), :]
        nrow(r) == 0 && continue
        d = Date.(r.date); v = parse.(Float64, r.prevalence); o = sortperm(d); d = d[o]; v = v[o]
        v = any(v .> 1) ? v ./ 100 : v
        T, S = held(model_site, d[1], d[end])
        sid = "me_" * lowercase(replace(site, r"[^A-Za-z]+" => "_"))
        push!(sites, Site(sid, "me", "me_" * model_site, T, S, start_state(v[1], 10), dayidx.(d[2:end], d[1]), v[2:end], true,
                          Int[], Float64[], persist(v)[2:end]))
    end
end
logln("sites by system: ", Dict(sy => count(s -> s.system == sy, sites) for sy in SYSTEMS),
      "; groups: ", length(unique(s.group for s in sites)))

# ---- learned curves ----
nparams(H) = 3 * H + 1
softplus(x) = log1p(exp(-abs(x))) + max(x, zero(x))
function curve(w, x, H)
    W1 = view(w, 1:H); b1 = view(w, H+1:2H); W2 = view(w, 2H+1:3H); b2 = w[3H+1]
    return softplus(sum(W2 .* tanh.(W1 .* x .+ b1)) + b2)
end
function init_curve(rng, H)
    w = 0.5 .* randn(rng, nparams(H)); w[2H+1:3H] .= 0.0; w[end] = log(exp(1.0) - 1); return w
end
# θ layout: [log β, log μ, logit host, eps_md, eps_gb, eps_ne, eps_me, gS (3H+1), gT (3H+1)]
const NMECH = 7
eps_index(sy) = 3 + findfirst(==(sy), SYSTEMS)
gS_w(θ, H) = view(θ, NMECH+1:NMECH+nparams(H))
gT_w(θ, H) = view(θ, NMECH+nparams(H)+1:NMECH+2nparams(H))
logistic(x) = 1 / (1 + exp(-x))
function init_theta(H; seed=20261015)
    rng = Xoshiro(seed)
    return vcat([-21.0, -0.5, 4.0, log(2e-4), log(1e-3), log(1e-3), log(1e-3)], init_curve(rng, H), init_curve(rng, H))
end

function predict(θ, s::Site, H)
    p = EstablishmentParams(lab, pw; beta=exp(θ[1]), mu=exp(θ[2]), epsilon=exp(θ[eps_index(s.system)]), host=logistic(θ[3]))
    wS = gS_w(θ, H); wT = gT_w(θ, H)
    g(t, sal) = curve(wS, sal / 30, H) * curve(wT, (t - 15) / 10, H)
    rec = sort(unique(vcat(s.prev_rec, s.wp_rec)))
    r = simulate_establishment(p, s.T, s.S, s.N0, 0.0; background_mortality_d=BACKGROUND_D, record=rec,
                               transmission_mod=g, smooth=true)
    at = Dict(d => i for (i, d) in enumerate(rec))
    prev = [(s.pcr ? r.true_prev : r.prev)[at[d]] for d in s.prev_rec]
    wp = [r.mackin[at[d]] for d in s.wp_rec]
    return prev, wp
end

function penalty(θ, H, L2)
    wS = gS_w(θ, H); wT = gT_w(θ, H)
    mono = zero(eltype(θ))
    for k in 0:34
        d = curve(wS, (k + 1) / 30, H) - curve(wS, k / 30, H); mono += min(d, zero(d))^2
    end
    for k in 0:29
        d = curve(wT, (k + 1 - 15) / 10, H) - curve(wT, (k - 15) / 10, H); mono += min(d, zero(d))^2
    end
    reg = L2 * (sum(abs2, view(wS, 1:length(wS)-1)) + sum(abs2, view(wT, 1:length(wT)-1)))
    return MONO_W * mono + reg
end

function data_loss(θ, ss, H)
    tot = zero(eltype(θ)); nsys = 0
    for sy in SYSTEMS
        sub = [s for s in ss if s.system == sy]
        isempty(sub) && continue
        ep = zero(eltype(θ)); np = 0; ew = zero(eltype(θ)); nw = 0
        for s in sub
            prev, wp = predict(θ, s, H)
            ep += sum(abs2, prev .- s.prev_obs); np += length(prev)
            ew += sum(abs2, (wp .- s.wp_obs) ./ 5); nw += length(wp)
        end
        tot += ep / max(np, 1) + (nw > 0 ? 0.5 * ew / nw : 0.0); nsys += 1
    end
    return tot / nsys
end

function train(ss, cfg)
    θ0 = init_theta(cfg.H)
    f = OptimizationFunction((θ, _) -> data_loss(θ, ss, cfg.H) + penalty(θ, cfg.H, cfg.L2), Optimization.AutoForwardDiff())
    return solve(OptimizationProblem(f, θ0), LBFGS(); maxiters=MAXITERS).u
end

# ---- folds ----
groups = sort(unique(s.group for s in sites))
DRY && (groups = [groups[1], groups[end]])
"Inner folds: round-robin assignment of site groups to N_INNER folds within each system."
function inner_folds(train_sites)
    fold_of = Dict{String,Int}()
    for sy in SYSTEMS
        gs = sort(unique(s.group for s in train_sites if s.system == sy))
        for (i, g) in enumerate(gs)
            fold_of[g] = mod1(i, N_INNER)
        end
    end
    return fold_of
end

rows = Vector{NamedTuple}()
fold_info = Vector{NamedTuple}()
logln("outer folds: ", length(groups), "; trainings: ", length(groups) * (length(CONFIGS) * N_INNER + 1))
tstart = time()
Threads.@threads for gi in eachindex(groups)
    g = groups[gi]
    test = [s for s in sites if s.group == g]
    trainset = [s for s in sites if s.group != g]
    fo = inner_folds(trainset)
    vloss = zeros(length(CONFIGS))
    for (ci, cfg) in enumerate(CONFIGS), k in 1:N_INNER
        tr_ = [s for s in trainset if fo[s.group] != k]; va = [s for s in trainset if fo[s.group] == k]
        θ = train(tr_, cfg)
        vloss[ci] += data_loss(θ, va, cfg.H) / N_INNER
    end
    best = argmin(vloss); cfg = CONFIGS[best]
    θ = train(trainset, cfg)
    for s in test
        prev, _ = predict(θ, s, cfg.H)
        sysmean = mean(vcat([x.prev_obs for x in trainset if x.system == s.system]...))
        lock(LOCK) do
            for i in eachindex(prev)
                push!(rows, (site=s.id, system=s.system, group=g, day=s.prev_rec[i], obs=s.prev_obs[i], pred=prev[i],
                             persist=s.prev_persist[i], sysmean=sysmean))
            end
            push!(fold_info, (group=g, H=cfg.H, L2=cfg.L2, val_losses=join(round.(vloss; sigdigits=4), ";")))
        end
    end
    logln("fold ", g, " done (config H=", cfg.H, ", L2=", cfg.L2, "); elapsed ", round((time() - tstart) / 60; digits=1), " min")
end

pred = DataFrame(rows)
Arrow.write(joinpath(OUT, "predictions.arrow"), pred)
Arrow.write(joinpath(OUT, "folds.arrow"), DataFrame(fold_info))

rmse(a, b) = sqrt(mean((a .- b) .^ 2))
function boot_diff(d; B=2000, rng=Xoshiro(1))
    ss = unique(d.site)
    diffs = Float64[]
    for _ in 1:B
        pick = rand(rng, ss, length(ss))
        dd = vcat([d[d.site .== s, :] for s in pick]...)
        push!(diffs, rmse(dd.pred, dd.obs) - rmse(Float64.(dd.persist), dd.obs))
    end
    return quantile(diffs, [0.025, 0.975])
end
fmt(x; d=3) = string(round(x; digits=d))
open(joinpath(OUT, "nested.md"), "w") do io
    println(io, "# P4.2g nested cross-validation, establishment model with learned transmission curves (scripts/p4_2g_nested.jl)\n")
    println(io, "Sites: ", join(["$(sy) $(count(s -> s.system == sy, sites))" for sy in SYSTEMS], ", "), "; outer folds (site groups): $(length(groups)). ",
            "Configurations: ", join(["H=$(c.H), L2=$(c.L2)" for c in CONFIGS], "; "), ". Chosen per fold: see folds.arrow.\n")
    println(io, "| system | sites | scored rows | RMSE model | RMSE persistence | RMSE training-system mean | Pearson | 95 % CI of (model - persistence), site bootstrap | pass |")
    println(io, "|---|---|---|---|---|---|---|---|---|")
    for sy in SYSTEMS
        d = pred[pred.system .== sy, :]
        nrow(d) == 0 && continue
        dp = d[.!ismissing.(d.persist), :]
        ci = boot_diff(dp)
        pass = rmse(dp.pred, dp.obs) < rmse(Float64.(dp.persist), dp.obs) && ci[2] < 0
        println(io, "| $sy | $(length(unique(d.site))) | $(nrow(dp)) | $(fmt(rmse(dp.pred, dp.obs))) | $(fmt(rmse(Float64.(dp.persist), dp.obs))) | ",
                "$(fmt(rmse(dp.sysmean, dp.obs))) | $(fmt(cor(dp.pred, dp.obs))) | $(fmt(ci[1])) to $(fmt(ci[2])) | $(pass) |")
    end
    println(io, "\nPer site:\n\n| site | system | rows | RMSE model | RMSE persistence |\n|---|---|---|---|---|")
    for g in groupby(pred[.!ismissing.(pred.persist), :], :site)
        println(io, "| $(g.site[1]) | $(g.system[1]) | $(nrow(g)) | $(fmt(rmse(g.pred, g.obs))) | $(fmt(rmse(Float64.(g.persist), g.obs))) |")
    end
end
logln(read(joinpath(OUT, "nested.md"), String))
close(LOG)
