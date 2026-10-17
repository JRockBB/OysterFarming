# P4.2i: staged compartmental Dermo model (SIS-type with a detection threshold), per-site transmission rates, what
# explains them, the self-sustaining threshold for the Damariscotta, and the same nested held-out test as P4.2g/h.
#
# A first cut with four unstaged states (S, L, D, dead) failed two ways and was replaced (2026-10-17): exponential
# delays smeared the seasonal timing (Great Bay in-sample 0.37 vs 0.19 for the staged establishment model), and with
# every infected oyster equally infectious and no release at death, warming raised the self-sustaining threshold
# (faster death removed sources faster than it added transmission). This version keeps the compartment idea but
# stages the delays and weights infectiousness by parasite load.
#
# States (fractions of a constant-density population; dead oysters are replaced by susceptibles): S susceptible;
# L1..L4 infected below RFTM detection (loads 1 to 2^12 cells per oyster, 3 doublings per stage); D1..D4 detectable
# (2^12 to 2^26 cells, 3.5 doublings per stage; the lethal load of size class 7, Powell et al. 1996). Daily transitions:
#   progression up one stage at NSUB * max(r, 0) / ln(span of the compartment), regression down one stage at
#     NSUB * max(-r, 0) / ln(span), with r(T, S) the lab-anchored net in-host rate (params/dermo_establishment_lab.toml,
#     host factor 1); regression below L1 is clearance to S; progression out of D4 is death (replaced by S);
#   clearance of light infections L -> S at c_system (fitted);
#   force of infection lambda = beta_site * fT(T) * fS(S) * I_eff + epsilon_system (beta, epsilon fitted), where
#     fT = 0 at <= 4 degC (zoospores die within a day at 4 degC, Chu & Greene 1989), linear to 1 at 20 degC;
#     fS = 0 at <= S0, linear to 1 at S1 (the (S0, S1) pair chosen from a grid by total likelihood; the learned term of
#       P4.2f and the Ragone Calvo et al. 2001 transmission term both gate transmission by salinity);
#     I_eff = sum over stages of (cells^0.667 / cells_top^0.667) x stage fraction + burst x deaths that day:
#       shedding scales with cells^0.667 (Bushek et al. 2002) and a dying oyster releases its whole load (2^26 cells),
#       equal to burst = 2^26 / (10^-0.752 x cells_top^0.667) days of shedding by a top-stage oyster (about 5,000).
# Observables: RFTM prevalence = sum of D; qPCR prevalence (Maine) = sum of L and D.
# Likelihood: binomial with the survey sample sizes (Maryland 30, MD DNR 2024 report p. 14; Great Bay as transcribed;
# New England 25; Maine 24); counts = round(prevalence x n). Parameters are bounded through a logistic map.
# Stage A: per system, joint fit of per-site (beta, I0) and system (c, epsilon) by LBFGS + forward-mode AD, repeated
#   over the salinity-ramp grid; Great Bay import starts at the P4.2e arrival (2003.15).
# Stage B: OLS of log beta_site on site covariates (mean Jun-Oct salinity, mean annual degree-days above 20 degC, mean
#   observed prevalence of other sites within 30 km or the system), with and without system dummies; R2 and
#   leave-one-site-out R2.
# Stage C: Damariscotta threshold: critical beta at which a vanishing introduction neither grows nor decays over 10
#   cycled years in the linearised system (epsilon = 0, cool-system median clearance), at +0/0.5/1/2 degC.
# Stage D: nested test, same outer folds (site groups); per fold the salinity ramp and all fits use training sites
#   only; inner leave-one-site-out chooses between the two covariate models; a held-out site gets beta from the
#   covariate model, c and epsilon from its system's training fit, I0 = its first observation (also given to the
#   persistence baseline); scored on later observations; pass rule as P4.2g.
#
# Usage: julia --project -t 14 scripts/p4_2i_compartment.jl [stages, default ABCD]
# Writes outputs/p4_2i/{sites.arrow, compartment.md, nested.md, predictions.arrow}.

using OysterFarming
using Arrow, DataFrames, Dates, ForwardDiff, LinearAlgebra, Random, Statistics
using Optimization, OptimizationOptimJL

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "outputs", "p4_2i")
const STAGES = isempty(ARGS) ? "ABCD" : ARGS[1]
mkpath(OUT)
include(joinpath(@__DIR__, "p4_2g_data.jl"))

const PLAB = EstablishmentParams(lab, pw; beta=0.0, mu=0.0, epsilon=0.0, host=1.0)
const NSUB = 4
const LN_DET = log(2.0^12)
const LN_LETHAL = log(2.0^26 / 2.0^12)
const T_TR0, T_TR1 = 4.0, 20.0
const GB_ARRIVAL = 2003.15
const NEAR_KM = 30.0
const SGRID = [(0.0, 4.0), (3.0, 10.0), (6.0, 12.0), (8.0, 16.0), (10.0, 20.0), (12.0, 25.0)]
# load-weighted infectiousness (Bushek et al. 2002: shedding ∝ cells^0.667), normalised to the top D stage
const LOG2_L = [3k - 1.5 for k in 1:NSUB]                 # mid-stage log2 cells, 0 to 12
const LOG2_D = [12 + 3.5k - 1.75 for k in 1:NSUB]         # 12 to 26
const WL = (2.0 .^ LOG2_L) .^ 0.667 ./ (2.0^LOG2_D[end])^0.667
const WD = (2.0 .^ LOG2_D) .^ 0.667 ./ (2.0^LOG2_D[end])^0.667
const BURST = 2.0^26 / (10.0^-0.752 * (2.0^LOG2_D[end])^0.667)
ftr(T) = clamp((T - T_TR0) / (T_TR1 - T_TR0), 0.0, 1.0)
fsal(S, sal) = clamp((S - sal[1]) / (sal[2] - sal[1]), 0.0, 1.0)
logistic(x) = 1 / (1 + exp(-x))

struct Drivers
    f::Vector{Float64}        # transmission temperature factor
    pL::Vector{Float64}       # L stage progression fraction per day
    rL::Vector{Float64}       # L stage regression
    pD::Vector{Float64}       # D stage progression (out of D4 = death)
    rD::Vector{Float64}       # D stage regression (out of D1 -> L4)
    intro::Int
end
function drivers(s::Site)
    r = [establishment_rate(t, sal, PLAB) for (t, sal) in zip(s.T, s.S)]
    g = max.(r, 0.0); d = max.(-r, 0.0)
    intro = s.system == "gb" ? clamp(round(Int, (GB_ARRIVAL - year(s.start)) * 365.25) + 1, 1, length(s.T)) : 1
    return Drivers(ftr.(s.T), 1 .- exp.(-NSUB .* g ./ LN_DET), 1 .- exp.(-NSUB .* d ./ LN_DET),
                   1 .- exp.(-NSUB .* g ./ LN_LETHAL), 1 .- exp.(-NSUB .* d ./ LN_LETHAL), intro)
end
const DRV = Dict(s.id => drivers(s) for s in sites)

"Simulate one site; returns the observable (sum of D, or L + D for qPCR sites) at the recorded days."
function simulate(s::Site, beta, I0, c, eps, sal; rec=s.prev_rec, dv=DRV[s.id])
    R = promote_type(typeof(beta), typeof(I0), typeof(c), typeof(eps))
    L = zeros(R, NSUB); D = zeros(R, NSUB)
    S = one(R) - I0; D[1] = R(I0)
    cl = 1 - exp(-c)
    deaths = zero(R)
    newL = zeros(R, NSUB); newD = zeros(R, NSUB)
    out = zeros(R, length(rec)); ri = 1
    for d in eachindex(s.T)
        ieff = dot(WL, L) + dot(WD, D) + BURST * deaths
        lam = beta * dv.f[d] * fsal(s.S[d], sal) * ieff + (d >= dv.intro ? eps : zero(R))
        inS = (1 - exp(-lam)) * S
        pL = dv.pL[d]; rL = dv.rL[d]; pD = dv.pD[d]; rD = dv.rD[d]
        toS = zero(R)
        fill!(newL, zero(R)); fill!(newD, zero(R))
        for k in 1:NSUB
            up = pL * L[k]; down = rL * L[k]; clr = cl * L[k]
            newL[k] += L[k] - up - down - clr
            k < NSUB ? (newL[k+1] += up) : (newD[1] += up)
            k > 1 ? (newL[k-1] += down) : (toS += down)        # regression below one cell is clearance
            toS += clr
        end
        died = zero(R)
        for k in 1:NSUB
            up = pD * D[k]; down = rD * D[k]
            newD[k] += D[k] - up - down
            k < NSUB ? (newD[k+1] += up) : (died += up)
            k > 1 ? (newD[k-1] += down) : (newL[NSUB] += down)
        end
        newL[1] += inS
        S += -inS + toS + died
        L .= newL; D .= newD; deaths = died
        while ri <= length(rec) && rec[ri] == d
            out[ri] = s.pcr ? sum(L) + sum(D) : sum(D); ri += 1
        end
    end
    return out
end

function nll(s::Site, pred)
    tot = zero(eltype(pred))
    for (p, obs, n) in zip(pred, s.prev_obs, s.prev_n)
        k = round(obs * n); q = clamp(p, 1e-6, 1 - 1e-6)
        tot -= k * log(q) + (n - k) * log(1 - q)
    end
    return tot
end

# ---- Stage A: per-system joint fits ----
# Bounded parameters through a logistic map (an unbounded log scale let beta run to 1e70 where the likelihood is flat).
const BOX = (beta=(1e-5, 10.0), I0=(1e-3, 0.999), c=(1e-5, 0.2), eps=(1e-7, 1e-2))
boxmap(x, b) = b[1] * (b[2] / b[1])^logistic(x)
boxinv(v, b) = (q = log(v / b[1]) / log(b[2] / b[1]); log(q / (1 - q)))
function unpack(θ, m)
    return (beta=[boxmap(θ[i], BOX.beta) for i in 1:m], I0=[boxmap(θ[m+i], BOX.I0) for i in 1:m],
            c=boxmap(θ[2m+1], BOX.c), eps=boxmap(θ[2m+2], BOX.eps))
end
function fit_system(ss::Vector{Site}, sal; iters=400)
    m = length(ss)
    function f(θ, _)
        u = unpack(θ, m)
        return sum(nll(ss[i], simulate(ss[i], u.beta[i], u.I0[i], u.c, u.eps, sal)) for i in 1:m)
    end
    θ0 = vcat(fill(boxinv(0.01, BOX.beta), m), [boxinv(clamp(s.prev_obs[1], 0.02, 0.98), BOX.I0) for s in ss],
              [boxinv(0.005, BOX.c), boxinv(1e-4, BOX.eps)])
    sol = solve(OptimizationProblem(OptimizationFunction(f, Optimization.AutoForwardDiff()), θ0), LBFGS(); maxiters=iters)
    u = unpack(sol.u, m)
    return (beta=u.beta, I0=u.I0, c=u.c, eps=u.eps, nll=f(sol.u, nothing))
end
"Fit every system over the salinity grid; returns the best grid point, its fits and the nll profile."
function fit_all(ss::Vector{Site})
    best = nothing; profile = Float64[]
    for sal in SGRID
        fits = Dict(sy => fit_system([s for s in ss if s.system == sy], sal) for sy in SYSTEMS if any(s.system == sy for s in ss))
        tot = sum(f.nll for f in values(fits)); push!(profile, tot)
        (best === nothing || tot < best.nll) && (best = (sal=sal, fits=fits, nll=tot))
    end
    return best, profile
end

rmse(a, b) = sqrt(mean((a .- b) .^ 2))

# ---- Stage B: site covariates ----
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
function covariates(s::Site)
    dates = [s.start + Day(d - 1) for d in eachindex(s.T)]
    summer = [month(d) in 6:10 for d in dates]
    yrs = [y for y in unique(year.(dates)) if count(==(y), year.(dates)) >= 365]
    dd = isempty(yrs) ? sum(max.(s.T .- 20, 0.0)) * 365 / length(s.T) : mean(sum(max.(s.T[year.(dates) .== y] .- 20, 0.0)) for y in yrs)
    others = [(o, km(site_coord(s), site_coord(o))) for o in sites if o.system == s.system && o.id != s.id]
    near = [o for (o, d) in others if d <= NEAR_KM]
    pool = isempty(near) ? [o for (o, _) in others] : near
    neigh = isempty(pool) ? NaN : mean(vcat([o.prev_obs for o in pool]...))
    return (S_summer=mean(s.S[summer]), DD20_yr=dd, neigh_prev=neigh, local_neigh=!isempty(near))
end
const COV = Dict(s.id => covariates(s) for s in sites)
const SYSOF = Dict(s.id => s.system for s in sites)
ols(X, y) = (β = X \ y; (β=β, r2=1 - sum(abs2, y .- X * β) / sum(abs2, y .- mean(y))))
function design(ids, with_system)
    X = hcat(ones(length(ids)), [COV[i].S_summer for i in ids], [COV[i].DD20_yr for i in ids] ./ 100, [COV[i].neigh_prev for i in ids])
    if with_system
        for sy in SYSTEMS[2:end]
            X = hcat(X, [SYSOF[i] == sy ? 1.0 : 0.0 for i in ids])
        end
    end
    return X
end
function loo_r2(ids, y, with_system)
    pred = similar(y)
    for i in eachindex(ids)
        tr = setdiff(eachindex(ids), i)
        β = design(ids[tr], with_system) \ y[tr]
        pred[i] = (design(ids[i:i], with_system) * β)[1]
    end
    return 1 - sum(abs2, y .- pred) / sum(abs2, y .- mean(y)), pred
end

# ---- Stage C: Damariscotta threshold (linearised growth of a vanishing introduction) ----
full_years(b) = [y for y in unique(year.(b.date)) if count(year.(b.date) .== y) == daysinyear(y)]
function maine_site(site, dT)
    b = longest_block(forcing_for_site(forcing, site))
    f = cycle_years(b, full_years(b), 10)
    T = Float64.(f.temp_c) .+ dT; S = Float64.(f.salinity_psu)
    return Site("thr_" * site, "me", "me", f.date[1], T, S, zeros(EST_N_LEVELS), [length(T)], [0.0], [1.0], false,
                Int[], Float64[], Union{Missing,Float64}[missing])
end
function growth_rate(s::Site, beta, c, sal)
    dv = drivers(s); cl = 1 - exp(-c)
    L = zeros(NSUB); D = zeros(NSUB); L[1] = 1e-6; deaths = 0.0
    newL = zeros(NSUB); newD = zeros(NSUB)
    for d in eachindex(s.T)
        ieff = dot(WL, L) + dot(WD, D) + BURST * deaths
        inS = 1 - exp(-beta * dv.f[d] * fsal(s.S[d], sal) * ieff)       # S = 1
        pL = dv.pL[d]; rL = dv.rL[d]; pD = dv.pD[d]; rD = dv.rD[d]
        fill!(newL, 0.0); fill!(newD, 0.0)
        for k in 1:NSUB
            up = pL * L[k]; down = rL * L[k]; clr = cl * L[k]
            newL[k] += L[k] - up - down - clr
            k < NSUB ? (newL[k+1] += up) : (newD[1] += up)
            k > 1 && (newL[k-1] += down)
        end
        died = 0.0
        for k in 1:NSUB
            up = pD * D[k]; down = rD * D[k]
            newD[k] += D[k] - up - down
            k < NSUB ? (newD[k+1] += up) : (died += up)
            k > 1 ? (newD[k-1] += down) : (newL[NSUB] += down)
        end
        newL[1] += inS
        L .= newL; D .= newD; deaths = died
        tot = sum(L) + sum(D)
        tot > 1e3 && return 1e3            # unmistakable growth; stop early (numbers are linear, scale is arbitrary)
    end
    return log(max(sum(L) + sum(D), 1e-300) / 1e-6)
end
function critical_beta(s::Site, c, sal)
    lo, hi = log(1e-5), log(100.0)
    growth_rate(s, exp(hi), c, sal) <= 0 && return Inf
    growth_rate(s, exp(lo), c, sal) >= 0 && return 0.0
    for _ in 1:40
        mid = (lo + hi) / 2
        growth_rate(s, exp(mid), c, sal) > 0 ? (hi = mid) : (lo = mid)
    end
    return exp((lo + hi) / 2)
end

# ---- Stage D: nested test ----
groups = sort(unique(s.group for s in sites))
function stage_D()
    rows = Vector{NamedTuple}(); lk = ReentrantLock()
    Threads.@threads for gi in eachindex(groups)
        g = groups[gi]
        test = [s for s in sites if s.group == g]; trainset = [s for s in sites if s.group != g]
        best, _ = fit_all(trainset)
        ids = String[]; lb = Float64[]
        for (sy, f) in best.fits, (i, s) in enumerate([x for x in trainset if x.system == sy])
            push!(ids, s.id); push!(lb, log(f.beta[i]))
        end
        r2_pool, _ = loo_r2(ids, lb, false); r2_sys, _ = loo_r2(ids, lb, true)
        with_sys = r2_sys > r2_pool
        β = design(ids, with_sys) \ lb
        for s in test
            haskey(best.fits, s.system) || continue
            beta_hat = exp((design([s.id], with_sys) * β)[1])
            f = best.fits[s.system]
            pred = simulate(s, beta_hat, clamp(s.prev_obs[1], 0.005, 0.995), f.c, f.eps, best.sal)
            sysmean = mean(vcat([x.prev_obs for x in trainset if x.system == s.system]...))
            lock(lk) do
                for i in 2:length(pred)
                    push!(rows, (site=s.id, system=s.system, group=g, obs=s.prev_obs[i], pred=pred[i], persist=s.prev_obs[i-1],
                                 sysmean=sysmean, beta_hat=beta_hat, with_system=with_sys, sal0=best.sal[1], sal1=best.sal[2]))
                end
            end
        end
        println("fold $g done (sal ", best.sal, ", with_system ", with_sys, ")"); flush(stdout)
    end
    return DataFrame(rows)
end
function boot_diff(d; B=2000, rng=Xoshiro(1))
    ss = unique(d.site); diffs = Float64[]
    for _ in 1:B
        dd = vcat([d[d.site .== s, :] for s in rand(rng, ss, length(ss))]...)
        push!(diffs, rmse(dd.pred, dd.obs) - rmse(dd.persist, dd.obs))
    end
    return quantile(diffs, [0.025, 0.975])
end

fmt(x; d=3) = (x isa Real && (isnan(x) || isinf(x))) ? string(x) : string(round(x; digits=d))
if occursin("A", STAGES)
    t = @elapsed (best, profile) = fit_all(sites)
    println("stage A: salinity ramp ", best.sal, " nll ", round(best.nll; digits=1), "; profile ", round.(profile; digits=0), "; ", round(t; digits=0), " s")
    fits = DataFrame(site=String[], system=String[], beta=Float64[], I0=Float64[], c=Float64[], eps=Float64[],
                     n_obs=Int[], rmse_model=Float64[], rmse_persist=Float64[])
    for sy in SYSTEMS
        haskey(best.fits, sy) || continue
        ss = [s for s in sites if s.system == sy]; f = best.fits[sy]
        for (i, s) in enumerate(ss)
            pred = simulate(s, f.beta[i], f.I0[i], f.c, f.eps, best.sal)
            ok = .!ismissing.(s.prev_persist)
            push!(fits, (s.id, sy, f.beta[i], f.I0[i], f.c, f.eps, length(pred), rmse(pred, s.prev_obs),
                         count(ok) > 0 ? rmse(Float64.(s.prev_persist[ok]), s.prev_obs[ok]) : NaN))
        end
    end
    for k in (:S_summer, :DD20_yr, :neigh_prev)
        fits[!, k] = [getfield(COV[i], k) for i in fits.site]
    end
    Arrow.write(joinpath(OUT, "sites.arrow"), fits)
    open(joinpath(OUT, "compartment.md"), "w") do io
        println(io, "# P4.2i staged compartmental Dermo model (scripts/p4_2i_compartment.jl)\n")
        println(io, "## Stage A: per-system fits\n\nSalinity ramp on transmission chosen by total likelihood: 0 at $(best.sal[1]) ppt, 1 at $(best.sal[2]) ppt ",
                "(grid nll: ", join(["$(g[1])-$(g[2]): $(round(p; digits=0))" for (g, p) in zip(SGRID, profile)], "; "), ").\n")
        println(io, "| system | sites | c (d^-1) | epsilon (d^-1) | beta q10 / median / q90 (d^-1) | in-sample RMSE model | RMSE persistence | RMSE system mean |\n|---|---|---|---|---|---|---|---|")
        for g in groupby(fits, :system)
            q = quantile(g.beta, [0.1, 0.5, 0.9]); ok = .!isnan.(g.rmse_persist)
            allobs = vcat([s.prev_obs for s in sites if s.system == g.system[1]]...)
            println(io, "| $(g.system[1]) | $(nrow(g)) | $(round(g.c[1]; sigdigits=3)) | $(round(g.eps[1]; sigdigits=3)) | ",
                    "$(round(q[1]; sigdigits=2)) / $(round(q[2]; sigdigits=2)) / $(round(q[3]; sigdigits=2)) | ",
                    "$(fmt(sqrt(sum(g.rmse_model .^ 2 .* g.n_obs) / sum(g.n_obs)))) | $(fmt(sqrt(sum(g.rmse_persist[ok] .^ 2 .* g.n_obs[ok]) / sum(g.n_obs[ok])))) | $(fmt(std(allobs; corrected=false))) |")
        end
        if occursin("B", STAGES)
            println(io, "\n## Stage B: what explains log beta across sites\n")
            ids = fits.site; lb = log.(fits.beta)
            for (label, ws) in (("covariates only", false), ("covariates + system dummies", true))
                o = ols(design(ids, ws), lb); r2l, _ = loo_r2(ids, lb, ws)
                names = vcat(["intercept", "S_summer (ppt)", "DD20 per 100", "neighbour prevalence"], ws ? ["sys " * sy for sy in SYSTEMS[2:end]] : String[])
                println(io, "$label: R2 $(fmt(o.r2)), leave-one-site-out R2 $(fmt(r2l)); coefficients ",
                        join(["$(n) $(round(b; sigdigits=3))" for (n, b) in zip(names, o.β)], "; "), "\n")
            end
            println(io, "| system | corr(log beta, S_summer) | corr(log beta, DD20) | corr(log beta, neighbour prevalence) |\n|---|---|---|---|")
            for g in groupby(fits, :system)
                nrow(g) < 4 && continue
                lbg = log.(g.beta)
                println(io, "| $(g.system[1]) | $(fmt(cor(lbg, g.S_summer))) | $(fmt(cor(lbg, g.DD20_yr))) | $(fmt(cor(lbg, g.neigh_prev))) |")
            end
        end
        if occursin("C", STAGES)
            println(io, "\n## Stage C: self-sustaining threshold at the Damariscotta (epsilon = 0, linearised growth over 10 cycled years)\n")
            ccool = median([fits.c[fits.system .== sy][1] for sy in ("gb", "ne", "me") if any(fits.system .== sy)])
            println(io, "Clearance c = $(round(ccool; sigdigits=3)) d^-1 (median of the cool systems); salinity ramp as fitted. Critical beta (d^-1): the ",
                    "parasite spreads if the site's beta exceeds it.\n")
            println(io, "| site | +0 degC | +0.5 | +1 | +2 |\n|---|---|---|---|---|")
            for site in ("dam_lower", "dam_upper")
                vals = [critical_beta(maine_site(site, dT), ccool, best.sal) for dT in (0.0, 0.5, 1.0, 2.0)]
                println(io, "| $site | " * join([v == Inf ? "never" : string(round(v; sigdigits=2)) for v in vals], " | ") * " |")
            end
            println(io, "\nFitted beta for comparison (median by system): ", join(["$(sy) $(round(median(fits.beta[fits.system .== sy]); sigdigits=2))" for sy in SYSTEMS if any(fits.system .== sy)], "; "),
                    "; cool sites: ", join(["$(r.site) $(round(r.beta; sigdigits=2))" for r in eachrow(fits) if r.site in ("nh_nannie_island", "nh_adams_point", "ne_wellfleet_ma", "ne_cotuit_ma", "me_jack_s_point", "me_prentiss_island")], "; "), ".")
        end
    end
    println(read(joinpath(OUT, "compartment.md"), String))
end
if occursin("D", STAGES)
    pred = stage_D()
    Arrow.write(joinpath(OUT, "predictions.arrow"), pred)
    open(joinpath(OUT, "nested.md"), "w") do io
        println(io, "# P4.2i nested held-out test of the compartmental model with covariate-predicted beta\n")
        println(io, "| system | sites | rows | RMSE model | RMSE persistence | RMSE training-system mean | Pearson | 95 % CI (model - persistence) | pass |\n|---|---|---|---|---|---|---|---|---|")
        for sy in SYSTEMS
            d = pred[pred.system .== sy, :]; nrow(d) == 0 && continue
            ci = boot_diff(d)
            pass = rmse(d.pred, d.obs) < rmse(d.persist, d.obs) && ci[2] < 0
            println(io, "| $sy | $(length(unique(d.site))) | $(nrow(d)) | $(fmt(rmse(d.pred, d.obs))) | $(fmt(rmse(d.persist, d.obs))) | $(fmt(rmse(d.sysmean, d.obs))) | ",
                    "$(fmt(cor(d.pred, d.obs))) | $(fmt(ci[1])) to $(fmt(ci[2])) | $pass |")
        end
        folds = unique(pred[!, [:group, :with_system, :sal0, :sal1]])
        println(io, "\nPer fold: system dummies chosen in $(count(folds.with_system)) of $(nrow(folds)) folds; salinity ramps chosen: ",
                join(["$(k[1])-$(k[2]) ($(v))" for (k, v) in sort(collect(Dict((r.sal0, r.sal1) => count((folds.sal0 .== r.sal0) .& (folds.sal1 .== r.sal1)) for r in eachrow(folds))))], ", "), ".\n")
        println(io, "Per site:\n\n| site | system | rows | RMSE model | RMSE persistence | predicted beta |\n|---|---|---|---|---|---|")
        for g in groupby(pred, :site)
            println(io, "| $(g.site[1]) | $(g.system[1]) | $(nrow(g)) | $(fmt(rmse(g.pred, g.obs))) | $(fmt(rmse(g.persist, g.obs))) | $(round(g.beta_hat[1]; sigdigits=2)) |")
        end
    end
    println(read(joinpath(OUT, "nested.md"), String))
end
