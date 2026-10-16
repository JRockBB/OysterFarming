# P4.2h step 5: nested cross-validation with the correction re-placed on the inputs selected by the residual diagnosis
# (scripts/p4_2h_residuals.jl, outputs/p4_2h/residuals.md).
#
# Design (fixed before running):
# - Data: scripts/p4_2g_data.jl (57 sites, 4 systems; Great Bay winter forcing from measured stations; Maine spin-up).
# - Mechanistic part: establishment model (smooth mode), lab in-host rates, waterborne pool, per-system import, Great Bay
#   import starting at the P4.2e arrival (2003.15, all folds). Fitted: log beta, log mu, logit host, 4 log epsilon.
# - Correction: a daily multiplier on the whole infection rate (waterborne + import), `rate_mod` in simulate_establishment:
#     g_d = gT(T_d) * gD(DD20_d) * gP(prev_d) * gQ(pressure_d) * gM(doy_d)
#   with inputs selected by the diagnosis (|Spearman| >= 0.2 in >= 2 systems):
#     T_d: same-day temperature; DD20_d: degree-days above 20 degC over the prior 150 days (DD15 dropped as collinear);
#     prev_d: the site's most recent observed prevalence before day d (held; before the first observation the pressure
#       value, else 0); pressure_d: previous-calendar-year mean prevalence of other sites in the system within 30 km
#       (system-wide where none; 0 where no data); doy_d: day of year as (sin, cos).
#   Each g is a small network with softplus output (1 -> H -> 1; gM is 2 -> H -> 1), initialised at 1. Shape constraints
#   (penalty): gP and gQ non-decreasing (more infection around -> more transmission). L2 on curve weights.
#   Salinity stays in the mechanistic in-host rates only (not selected by the diagnosis).
# - Maine observable: qPCR prevalence compared with model qPCR-detectable prevalence through a fitted smooth detection
#   threshold (log2 cells per oyster, one shared parameter); also scored as site-season means.
# - Loss, folds, inner selection, pass rule: as P4.2g. Outer folds = site groups sharing a forcing station; inner = 3
#   round-robin folds within each system over the training groups; configurations L2 in {1e-5, 1e-3} with H = 2 (gM: 3);
#   the configuration with the lowest mean inner validation loss is refit on all training sites and scored once on the
#   test group. Pass rule per system: held-out RMSE below persistence with the 95 % site-bootstrap interval of the
#   difference below zero. LBFGS + forward-mode AD, 80 iterations, fixed initialisation; per-training logging.
#
# Usage: julia --project -t 14 scripts/p4_2h_nested.jl        (P42H_DRY=1 for a smoke test: 2 folds, 1 config, 3 iterations)
# Writes outputs/p4_2h/{nested.md, predictions.arrow, folds.arrow, progress.log}.

using OysterFarming
using Arrow, DataFrames, Dates, ForwardDiff, LinearAlgebra, Random, Statistics
using Optimization, OptimizationOptimJL

const ROOT = normpath(joinpath(@__DIR__, ".."))
const DRY = get(ENV, "P42H_DRY", "0") == "1"
const OUT = joinpath(ROOT, "outputs", DRY ? "p4_2h_dry" : "p4_2h")
const MAXITERS = DRY ? 3 : 80
const H1, HM = 2, 3
const CONFIGS = DRY ? [(L2=1e-5,)] : [(L2=1e-5,), (L2=1e-3,)]
const N_INNER = 3
const MONO_W = 1.0
const NEAR_KM = 30.0
const WINDOW = 150
const GB_ARRIVAL = 2003.15
mkpath(OUT)
const LOG = open(joinpath(OUT, "progress.log"), "a")
const LOCK = ReentrantLock()
logln(x...) = lock(LOCK) do
    println(LOG, now(), " ", x...); flush(LOG); println(x...); flush(stdout)
end
include(joinpath(@__DIR__, "p4_2g_data.jl"))
logln("sites by system: ", Dict(sy => count(s -> s.system == sy, sites) for sy in SYSTEMS))

# ---- daily correction inputs per site ----
coords = Dict(r.site_id => (r.lat, r.lon) for r in eachrow(sites_tbl))
const NE_COORDS = Dict("ne_wellfleet_ma" => (41.93, -70.03), "ne_cotuit_ma" => (41.61, -70.44),
                       "ne_martha_s_vineyard_ma_butler_s_neck_" => (41.36, -70.55), "ne_oyster_bay_ny" => (40.88, -73.53),
                       "ne_milford_ct" => (41.21, -73.05), "ne_westport_ct" => (41.12, -73.37),
                       "ne_great_sound_nj" => (39.07, -74.78), "ne_delaware_bay_nj" => (39.12, -75.05))
site_coord(s::Site) = s.system == "ne" ? NE_COORDS[s.id] : s.system == "me" ? coords[replace(s.group, "me_" => "")] : coords[s.id]
function km(a, b)
    la1, lo1, la2, lo2 = deg2rad(a[1]), deg2rad(a[2]), deg2rad(b[1]), deg2rad(b[2])
    return 2 * 6371 * asin(sqrt(sin((la2 - la1) / 2)^2 + cos(la1) * cos(la2) * sin((lo2 - lo1) / 2)^2))
end
obsdate(s::Site, d) = s.start + Day(d - 1)
"Full observation history (date, prevalence) of a site: all survey years for Maryland and Great Bay; the Site's own
observations (plus the unscored start sample) for the literature sites."
function history(s::Site)
    if s.system in ("md", "gb")
        o = sort(dermo[dermo.site_id .== s.id, :], :date)
        return collect(zip(o.date, Float64.(o.prevalence)))
    end
    h = [(obsdate(s, d), v) for (d, v) in zip(s.prev_rec, s.prev_obs)]
    s.system == "ne" && pushfirst!(h, (s.start, 1 - s.N0[1] / DENSITY0))
    return h
end
hist = Dict(s.id => history(s) for s in sites)
function daily_inputs(s::Site)
    nd = length(s.T)
    c = cumsum(max.(s.T .- 20, 0.0))
    dd = [c[d] - (d > WINDOW ? c[d-WINDOW] : 0.0) for d in 1:nd]
    doy = [dayofyear(obsdate(s, d)) for d in 1:nd]
    others = [(o, km(site_coord(s), site_coord(o))) for o in sites if o.system == s.system && o.id != s.id]
    near = [o for (o, d) in others if d <= NEAR_KM]
    pool = isempty(near) ? [o for (o, _) in others] : near
    byyear = Dict{Int,Float64}()
    for y in year(s.start)-1:year(obsdate(s, nd))
        v = [p for o in pool for (dt, p) in hist[o.id] if year(dt) == y - 1]
        isempty(v) || (byyear[y] = mean(v))
    end
    pressure = [get(byyear, year(obsdate(s, d)), 0.0) for d in 1:nd]
    own = hist[s.id]
    prev = zeros(nd); j = 0
    for d in 1:nd
        dt = obsdate(s, d)
        while j < length(own) && own[j+1][1] < dt; j += 1; end
        prev[d] = j == 0 ? pressure[d] : own[j][2]
    end
    return (dd=dd ./ 300, sin=sin.(2pi .* doy ./ 365.25), cos=cos.(2pi .* doy ./ 365.25), prev=prev, pressure=pressure)
end
inputs = Dict(s.id => daily_inputs(s) for s in sites)

# ---- curves and parameters ----
np1(H) = 3H + 1                      # 1-input curve
npM(H) = 4H + 1                      # 2-input curve
softplus(x) = log1p(exp(-abs(x))) + max(x, zero(x))
curve1(w, x, H) = softplus(sum(view(w, 2H+1:3H) .* tanh.(view(w, 1:H) .* x .+ view(w, H+1:2H))) + w[3H+1])
curve2(w, x1, x2, H) = softplus(sum(view(w, 3H+1:4H) .* tanh.(view(w, 1:H) .* x1 .+ view(w, H+1:2H) .* x2 .+ view(w, 2H+1:3H))) + w[4H+1])
const NMECH = 8                      # log beta, log mu, logit host, 4 log epsilon, qPCR log2 threshold
const OFF = cumsum([NMECH, np1(H1), np1(H1), np1(H1), np1(H1)])   # ends of: mech, gT, gD, gP, gQ; gM follows
const NPAR = OFF[end] + npM(HM)
wT(θ) = view(θ, OFF[1]+1:OFF[2]); wD(θ) = view(θ, OFF[2]+1:OFF[3]); wP(θ) = view(θ, OFF[3]+1:OFF[4])
wQ(θ) = view(θ, OFF[4]+1:OFF[5]); wM(θ) = view(θ, OFF[5]+1:NPAR)
eps_index(sy) = 3 + findfirst(==(sy), SYSTEMS)
logistic(x) = 1 / (1 + exp(-x))
function init_theta(; seed=20261016)
    rng = Xoshiro(seed)
    c1() = (w = 0.5 .* randn(rng, np1(H1)); w[2H1+1:3H1] .= 0; w[end] = log(exp(1.0) - 1); w)
    cM() = (w = 0.5 .* randn(rng, npM(HM)); w[3HM+1:4HM] .= 0; w[end] = log(exp(1.0) - 1); w)
    return vcat([-21.0, -0.5, 4.0, log(2e-4), log(1e-3), log(1e-3), log(1e-3), 4.0], c1(), c1(), c1(), c1(), cM())
end

function predict(θ, s::Site)
    p = EstablishmentParams(lab, pw; beta=exp(θ[1]), mu=exp(θ[2]), epsilon=exp(θ[eps_index(s.system)]), host=logistic(θ[3]))
    x = inputs[s.id]
    g = [curve1(wT(θ), (s.T[d] - 15) / 10, H1) * curve1(wD(θ), x.dd[d], H1) * curve1(wP(θ), x.prev[d], H1) *
         curve1(wQ(θ), x.pressure[d], H1) * curve2(wM(θ), x.sin[d], x.cos[d], HM) for d in eachindex(s.T)]
    intro = s.system == "gb" ? clamp(round(Int, (GB_ARRIVAL - year(s.start)) * 365.25) + 1, 1, length(s.T)) : 1
    rec = sort(unique(vcat(s.prev_rec, s.wp_rec)))
    r = simulate_establishment(p, s.T, s.S, s.N0, 0.0; intro_day=intro, background_mortality_d=BACKGROUND_D, record=rec,
                               smooth=true, rate_mod=g, pcr_detection_log2=θ[8])
    at = Dict(d => i for (i, d) in enumerate(rec))
    prev = [(s.pcr ? r.pcr_prev : r.prev)[at[d]] for d in s.prev_rec]
    wp = [r.mackin[at[d]] for d in s.wp_rec]
    return prev, wp
end

function penalty(θ, L2)
    mono = zero(eltype(θ))
    for w in (wP(θ), wQ(θ)), k in 0:19
        d = curve1(w, (k + 1) / 20, H1) - curve1(w, k / 20, H1); mono += min(d, zero(d))^2
    end
    reg = zero(eltype(θ))
    for w in (wT(θ), wD(θ), wP(θ), wQ(θ), wM(θ))
        reg += sum(abs2, view(w, 1:length(w)-1))
    end
    return MONO_W * mono + L2 * reg
end

function data_loss(θ, ss)
    tot = zero(eltype(θ)); nsys = 0
    for sy in SYSTEMS
        sub = [s for s in ss if s.system == sy]
        isempty(sub) && continue
        ep = zero(eltype(θ)); np = 0; ew = zero(eltype(θ)); nw = 0
        for s in sub
            prev, wp = predict(θ, s)
            ep += sum(abs2, prev .- s.prev_obs); np += length(prev)
            ew += sum(abs2, (wp .- s.wp_obs) ./ 5); nw += length(wp)
        end
        tot += ep / max(np, 1) + (nw > 0 ? 0.5 * ew / nw : 0.0); nsys += 1
    end
    return tot / nsys
end

function train(ss, cfg, label)
    f = OptimizationFunction((θ, _) -> data_loss(θ, ss) + penalty(θ, cfg.L2), Optimization.AutoForwardDiff())
    t = @elapsed sol = solve(OptimizationProblem(f, init_theta()), LBFGS(); maxiters=MAXITERS)
    logln("  ", label, ": trained in ", round(t / 60; digits=1), " min, data loss ", round(data_loss(sol.u, ss); digits=4))
    return sol.u
end

groups = sort(unique(s.group for s in sites))
DRY && (groups = [groups[1], groups[end]])
function inner_folds(train_sites)
    fold_of = Dict{String,Int}()
    for sy in SYSTEMS
        for (i, g) in enumerate(sort(unique(s.group for s in train_sites if s.system == sy)))
            fold_of[g] = mod1(i, N_INNER)
        end
    end
    return fold_of
end

rows = Vector{NamedTuple}(); fold_info = Vector{NamedTuple}()
logln("parameters: ", NPAR, "; outer folds: ", length(groups), "; trainings: ", length(groups) * (length(CONFIGS) * N_INNER + 1))
tstart = time()
Threads.@threads for gi in eachindex(groups)
    g = groups[gi]
    test = [s for s in sites if s.group == g]; trainset = [s for s in sites if s.group != g]
    fo = inner_folds(trainset)
    vloss = zeros(length(CONFIGS))
    for (ci, cfg) in enumerate(CONFIGS), k in 1:N_INNER
        tr_ = [s for s in trainset if fo[s.group] != k]; va = [s for s in trainset if fo[s.group] == k]
        θ = train(tr_, cfg, "$g inner L2=$(cfg.L2) fold $k")
        vloss[ci] += data_loss(θ, va) / N_INNER
    end
    cfg = CONFIGS[argmin(vloss)]
    θ = train(trainset, cfg, "$g final L2=$(cfg.L2)")
    lock(LOCK) do
        for s in test
            prev, _ = predict(θ, s)
            sysmean = mean(vcat([x.prev_obs for x in trainset if x.system == s.system]...))
            for i in eachindex(prev)
                push!(rows, (site=s.id, system=s.system, group=g, date=obsdate(s, s.prev_rec[i]), obs=s.prev_obs[i],
                             pred=prev[i], persist=s.prev_persist[i], sysmean=sysmean))
            end
        end
        push!(fold_info, (group=g, L2=cfg.L2, val_losses=join(round.(vloss; sigdigits=4), ";"),
                          theta=join(round.(θ; sigdigits=5), ";")))
    end
    logln("fold ", g, " done (L2=", cfg.L2, "); elapsed ", round((time() - tstart) / 60; digits=1), " min")
end

pred = DataFrame(rows)
Arrow.write(joinpath(OUT, "predictions.arrow"), pred)
Arrow.write(joinpath(OUT, "folds.arrow"), DataFrame(fold_info))
rmse(a, b) = sqrt(mean((a .- b) .^ 2))
function boot_diff(d; B=2000, rng=Xoshiro(1))
    ss = unique(d.site); diffs = Float64[]
    for _ in 1:B
        dd = vcat([d[d.site .== s, :] for s in rand(rng, ss, length(ss))]...)
        push!(diffs, rmse(dd.pred, dd.obs) - rmse(Float64.(dd.persist), dd.obs))
    end
    return quantile(diffs, [0.025, 0.975])
end
fmt(x; d=3) = (x isa Real && isnan(x)) ? "n/a" : string(round(x; digits=d))
open(joinpath(OUT, "nested.md"), "w") do io
    println(io, "# P4.2h nested cross-validation, correction on diagnosed inputs (scripts/p4_2h_nested.jl)\n")
    println(io, "Sites: ", join(["$(sy) $(count(s -> s.system == sy, sites))" for sy in SYSTEMS], ", "), "; outer folds $(length(groups)); ",
            "parameters $(NPAR); configurations L2 in ", join([c.L2 for c in CONFIGS], ", "), ". Inputs: same-day T, DD20 (150 d), ",
            "site's previous observation, neighbours' previous-year prevalence, day of year.\n")
    println(io, "| system | sites | scored rows | RMSE model | RMSE persistence | RMSE training-system mean | Pearson | 95 % CI (model - persistence) | pass |")
    println(io, "|---|---|---|---|---|---|---|---|---|")
    for sy in SYSTEMS
        d = pred[pred.system .== sy, :]; nrow(d) == 0 && continue
        dp = d[.!ismissing.(d.persist), :]
        ci = boot_diff(dp)
        pass = rmse(dp.pred, dp.obs) < rmse(Float64.(dp.persist), dp.obs) && ci[2] < 0
        println(io, "| $sy | $(length(unique(d.site))) | $(nrow(dp)) | $(fmt(rmse(dp.pred, dp.obs))) | $(fmt(rmse(Float64.(dp.persist), dp.obs))) | ",
                "$(fmt(rmse(dp.sysmean, dp.obs))) | $(fmt(cor(dp.pred, dp.obs))) | $(fmt(ci[1])) to $(fmt(ci[2])) | $(pass) |")
    end
    me = pred[pred.system .== "me", :]
    if nrow(me) > 0
        println(io, "\nMaine site-season means (qPCR vs model qPCR-detectable prevalence):\n\n| site | year | observed mean | model mean |\n|---|---|---|---|")
        for g in groupby(transform(me, :date => ByRow(year) => :year), [:site, :year])
            println(io, "| $(g.site[1]) | $(g.year[1]) | $(fmt(mean(g.obs))) | $(fmt(mean(g.pred))) |")
        end
    end
    println(io, "\nPer site:\n\n| site | system | rows | RMSE model | RMSE persistence |\n|---|---|---|---|---|")
    for g in groupby(pred[.!ismissing.(pred.persist), :], :site)
        println(io, "| $(g.site[1]) | $(g.system[1]) | $(nrow(g)) | $(fmt(rmse(g.pred, g.obs))) | $(fmt(rmse(Float64.(g.persist), g.obs))) |")
    end
end
logln(read(joinpath(OUT, "nested.md"), String))
close(LOG)
