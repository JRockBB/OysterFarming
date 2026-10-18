# P4.2j: fit the staged establishment model to every observation type at once (per-oyster class counts, monthly
# intensities, bar means, qPCR prevalence and copy numbers), instead of prevalence alone.
#
# Observation types and likelihoods (counts come from the survey sample sizes):
#   classes  Great Bay (NH F&G class counts 1996-2014; later years when transcribed) and Mussel Watch Chesapeake sites
#            (per-oyster Mackin scores 1988-2009, nominal date 1 Jan of the fiscal year): multinomial over the classes
#            0, 0.5, 1, 2, 3, 4, 5. Model class probabilities: the level distribution at the survey day (28 parasite
#            load levels, `record_levels`) mapped level -> Mackin rank (Hofmann et al. 1995 eq. 23; 0 at or below
#            RFTM detection) -> class by the rank bins [0.75, 1.5, 2.5, 3.5, 4.5); Mussel Watch scores in thirds
#            are rounded to the nearest class. A 2 % uniform mixture guards against zero probabilities.
#   prev     every site: binomial on detectable prevalence (qPCR-detectable for Maine, fitted threshold).
#   wp       Ford & Smolowitz monthly weighted prevalence (0-5, zeros included): Gaussian on the model mean rank,
#            sd fitted per system (New England; Great Bay's WP is implied by its class counts and not used twice).
#   md_int   Maryland bar mean intensity (0-7 scale, zeros included): Gaussian on a + b x model mean rank, with a, b
#            and the sd fitted (no published 0-7 to 0-5 mapping).
#   copies   Maine log10 qPCR copies per ug DNA: Gaussian on a_q + b_q x model mean log2 cells among qPCR-positive
#            oysters, with a_q, b_q and the sd fitted.
# Model: simulate_establishment, smooth mode, allow_clearance = true (infections below one cell are lost), host factor
#   1; per site beta and I0 (bounded), per system pool decay mu and import epsilon; Great Bay import from 2003.15.
# Fit: per system, LBFGS + forward-mode AD on the summed negative log-likelihood of all its observation types.
# Report: in-sample prevalence RMSE vs persistence and system mean; within-site year-to-year skill (correlation of
#   site-demeaned predicted and observed prevalence; the quantity every earlier variant missed); class-distribution
#   fit; intensity fits; parameters at bounds.
#
# Usage: julia --project -t 14 scripts/p4_2j_fit.jl
# Writes outputs/p4_2j/{fit.md, sites.arrow, obs.arrow}.

using OysterFarming
using Arrow, DataFrames, Dates, ForwardDiff, LinearAlgebra, Random, Statistics
using Optimization, OptimizationOptimJL

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "outputs", "p4_2j")
mkpath(OUT)
include(joinpath(@__DIR__, "p4_2g_data.jl"))

const CLASSES = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0]
const EDGES = [0.75, 1.5, 2.5, 3.5, 4.5]
const MIX = 0.02
const GB_ARRIVAL = 2003.15
const PLAB = EstablishmentParams(lab, pw; beta=0.0, mu=0.0, epsilon=0.0, host=1.0)

# ---- extra observations ----
struct Obs
    day::Int
    kind::Symbol               # :prev, :classes, :wp, :md_int, :copies
    n::Float64
    value::Float64             # prevalence / wp / md_int / log10 copies (unused for :classes)
    counts::Vector{Float64}    # class counts (:classes)
end
toclass(v) = CLASSES[argmin(abs.(CLASSES .- v))]
obs = Dict{String,Vector{Obs}}()
for s in sites
    o = Obs[]
    for (d, v, n) in zip(s.prev_rec, s.prev_obs, s.prev_n)
        push!(o, Obs(d, :prev, n, v, Float64[]))
    end
    for (d, v) in zip(s.wp_rec, s.wp_obs)
        s.system == "ne" && push!(o, Obs(d, :wp, 25.0, v, Float64[]))
    end
    obs[s.id] = o
end
# Great Bay class counts from the pipeline source (transcribed NH F&G Table 2) and the later reports' CSV if present
gb_rows = Tuple{Date,String,Int,Vector{Int}}[]
for l in readlines(joinpath(ROOT, "data-pipeline", "src", "oyster_data", "sources", "nh_greatbay.py"))
    m = match(r"\(\"(\d{4}-\d{2}-\d{2})\", \"(nh_[a-z_]+)\", (\d+), \((\d+), (\d+), (\d+), (\d+), (\d+), (\d+)\)\)", l)
    m === nothing || push!(gb_rows, (Date(m[1]), m[2], parse(Int, m[3]), parse.(Int, [m[i] for i in 4:9])))
end
newer = joinpath(LIT, "nhfg_greatbay_dermo_2015_2019.csv")
if isfile(newer)
    t = read_csv_simple("nhfg_greatbay_dermo_2015_2019.csv")
    known = Set((d, s) for (d, s, _, _) in gb_rows)
    name_map = Dict("Adams Point" => "nh_adams_point", "Adams Pt." => "nh_adams_point", "Nannie Island" => "nh_nannie_island",
                    "Woodman Point" => "nh_woodman_point", "Oyster River" => "nh_oyster_river", "Squamscott River" => "nh_squamscott_river")
    # the CSV holds two full transcriptions (articles 1416 and 1441); keep the longest single series (1441, 1996-2019),
    # skip rows whose class counts do not add up to the printed prevalence (listed in the CSV header), and skip the
    # aquaculture lots (Table 3)
    bad = Set([("2004-11-18", "Nannie Island"), ("2008-10-22", "Nannie Island"), ("2008-10-27", "Squamscott River"),
               ("2016-10-07", "Oyster River"), ("2018-10-24", "Woodman Point")])
    for r in eachrow(t)
        (r.source_article == "1441" && haskey(name_map, r.site) && !((r.date, r.site) in bad)) || continue
        d = Date(r.date); sid = name_map[r.site]
        (d, sid) in known && continue
        push!(gb_rows, (d, sid, parse(Int, r.n), parse.(Int, [r.c0_5, r.c1, r.c2, r.c3, r.c4, r.c5])))
    end
end
sitemap = Dict(s.id => s for s in sites)
for (d, sid, n, c) in gb_rows
    haskey(sitemap, sid) || continue
    s = sitemap[sid]
    day = dayidx(d, s.start)
    1 <= day <= length(s.T) || continue
    push!(obs[sid], Obs(day, :classes, n, 0.0, vcat(n - sum(c), Float64.(c))))
end
# Maryland bar mean intensity (0-7)
for s in sites
    s.system == "md" || continue
    o = sort(dermo[(dermo.site_id .== s.id) .& (year.(dermo.date) .>= 1999) .& .!ismissing.(dermo.intensity_mean), :], :date)
    for r in eachrow(o)
        push!(obs[s.id], Obs(dayidx(r.date, s.start), :md_int, 30.0, Float64(r.intensity_mean), Float64[]))
    end
end
# Maine copies
if isfile(joinpath(LIT, MAINE_FILE))
    me = read_csv_simple(MAINE_FILE)
    for (site, sid) in (("Jack's Point", "me_jack_s_point"), ("Prentiss Island", "me_prentiss_island"))
        haskey(sitemap, sid) || continue
        r = me[occursin.(lowercase(site), lowercase.(me.site)), :]
        for x in eachrow(r)
            c = tryparse(Float64, x.intensity_copies_per_ug)
            (c === nothing || c <= 0) && continue
            push!(obs[sid], Obs(dayidx(Date(x.date), sitemap[sid].start), :copies, 24.0, log10(c), Float64[]))
        end
    end
end
# Mussel Watch sites (per-oyster scores): class counts + prevalence, nominal 1 January of the fiscal year
mw = read_csv_simple("mussel_watch_chesapeake_oysters.csv")
fsites = Set(forcing.site_id)
for code in sort(unique(mw.nst_site))
    sid = "mw_" * lowercase(code)
    sid in fsites || continue
    r = mw[mw.nst_site .== code, :]
    fys = sort(unique(parse.(Int, r.fiscal_year)))
    a = Date(1987, 1, 1)
    T, S = held(sid, a, Date(fys[end], 12, 31))
    o = Obs[]; prevs = Float64[]; days = Int[]; ns = Float64[]
    for fy in fys
        vals = parse.(Float64, r.dermo_value[parse.(Int, r.fiscal_year) .== fy])
        d = dayidx(Date(fy, 1, 1), a)
        counts = [count(==(cl), toclass.(vals)) for cl in CLASSES]
        push!(o, Obs(d, :classes, length(vals), 0.0, Float64.(counts)))
        push!(o, Obs(d, :prev, length(vals), count(>(0), vals) / length(vals), Float64[]))
        push!(prevs, count(>(0), vals) / length(vals)); push!(days, d); push!(ns, length(vals))
    end
    push!(sites, Site(sid, "mw", "mw_" * code, a, T, S, start_state(0.5, 15), days, prevs, ns, false, Int[], Float64[], persist(prevs)))
    obs[sid] = o
end
const SYS = ["md", "mw", "gb", "ne", "me"]
println("sites: ", Dict(sy => count(s -> s.system == sy, sites) for sy in SYS), "; observations: ",
        Dict(k => sum(count(o -> o.kind == k, v) for v in values(obs)) for k in (:prev, :classes, :wp, :md_int, :copies)))

# ---- level -> class mapping (constant) ----
const RANK = [est_mackin(k, PLAB) for k in 1:EST_N_LEVELS]
const CLASS_OF_LEVEL = [k == 1 || RANK[k] == 0 ? 1 : 2 + count(e -> RANK[k] >= e, EDGES) for k in 1:EST_N_LEVELS]   # index into CLASSES
const LOG2_CELLS = [k == 1 ? 0.0 : Float64(k - 2) for k in 1:EST_N_LEVELS]

# ---- parameters ----
const BOX = (beta=(1e-12, 1e-6), I0=(1e-3, 0.999), mu=(1e-3, 3.0), eps=(1e-7, 1e-2), sd=(0.05, 3.0), thr=(0.0, 14.0))
logistic(x) = 1 / (1 + exp(-x))
boxmap(x, b) = b[1] * (b[2] / b[1])^logistic(x)
boxinv(v, b) = (q = clamp(log(v / b[1]) / log(b[2] / b[1]), 1e-6, 1 - 1e-6); log(q / (1 - q)))
# per-system theta: [beta_1..m, I0_1..m, mu, eps, obs params...]; obs params: md: a, b, sd; ne: sd_wp; me: thr, a_q, b_q, sd_q
nobs(sy) = sy == "md" ? 3 : sy == "ne" ? 1 : sy == "me" ? 4 : 0
function unpack(θ, m, sy)
    o = θ[2m+3:end]
    (beta=[boxmap(θ[i], BOX.beta) for i in 1:m], I0=[boxmap(θ[m+i], BOX.I0) for i in 1:m], mu=boxmap(θ[2m+1], BOX.mu),
     eps=boxmap(θ[2m+2], BOX.eps),
     a=sy == "md" ? o[1] : sy == "me" ? o[2] : 0.0, b=sy == "md" ? o[2] : sy == "me" ? o[3] : 1.0,
     sd=sy == "md" ? boxmap(o[3], BOX.sd) : sy == "ne" ? boxmap(o[1], BOX.sd) : sy == "me" ? boxmap(o[4], BOX.sd) : 1.0,
     thr=sy == "me" ? BOX.thr[2] * logistic(o[1]) : -Inf)      # linear box 0-14 (log2 cells)
end
function init_theta(ss, sy)
    m = length(ss)
    θ = vcat(fill(boxinv(1e-9, BOX.beta), m), [boxinv(clamp(s.prev_obs[1], 0.02, 0.98), BOX.I0) for s in ss],
             [boxinv(0.5, BOX.mu), boxinv(1e-4, BOX.eps)])
    sy == "md" && (θ = vcat(θ, [0.0, 1.0, boxinv(0.7, BOX.sd)]))
    sy == "ne" && (θ = vcat(θ, [boxinv(0.5, BOX.sd)]))
    sy == "me" && (θ = vcat(θ, [log(6.0 / 8.0), 0.0, 0.3, boxinv(0.5, BOX.sd)]))   # thr 6 = 14 * logistic(log(6/8))
    return θ
end

function site_nll(s::Site, o::Vector{Obs}, u)
    p = EstablishmentParams(lab, pw; beta=u.beta, mu=u.mu, epsilon=u.eps, host=1.0)
    rec = sort(unique(x.day for x in o))
    intro = s.system == "gb" ? clamp(round(Int, (GB_ARRIVAL - year(s.start)) * 365.25) + 1, 1, length(s.T)) : 1
    r = simulate_establishment(p, s.T, s.S, s.N0, 0.0; intro_day=intro, background_mortality_d=BACKGROUND_D, record=rec,
                               smooth=true, allow_clearance=true, record_levels=true, pcr_detection_log2=u.thr)
    at = Dict(d => i for (i, d) in enumerate(rec))
    R = eltype(r.prev)
    tot = zero(R)
    for x in o
        i = at[x.day]
        if x.kind == :prev
            q = clamp(s.pcr ? r.pcr_prev[i] : r.prev[i], 1e-6, 1 - 1e-6)
            k = round(x.value * x.n)
            tot -= k * log(q) + (x.n - k) * log(1 - q)
        elseif x.kind == :classes
            pc = zeros(R, length(CLASSES))
            for k in 1:EST_N_LEVELS
                pc[CLASS_OF_LEVEL[k]] += r.levels[k, i]
            end
            pc = (1 - MIX) .* pc .+ MIX / length(CLASSES)
            tot -= sum(x.counts .* log.(pc))
        elseif x.kind == :wp
            tot += 0.5 * ((x.value - r.mackin[i]) / u.sd)^2 + log(u.sd)
        elseif x.kind == :md_int
            tot += 0.5 * ((x.value - (u.a + u.b * r.mackin[i])) / u.sd)^2 + log(u.sd)
        elseif x.kind == :copies
            pos = sum(r.levels[k, i] for k in 2:EST_N_LEVELS)
            mean_log2 = pos > 1e-9 ? sum(r.levels[k, i] * LOG2_CELLS[k] for k in 2:EST_N_LEVELS) / pos : zero(R)
            tot += 0.5 * ((x.value - (u.a + u.b * mean_log2)) / u.sd)^2 + log(u.sd)
        end
    end
    return tot
end

const ITERS = parse(Int, get(ENV, "P42J_ITERS", "300"))     # P42J_ITERS=2 for a smoke test
const PRIOR_W = 0.01
function fit_system(sy; iters=ITERS)
    ss = [s for s in sites if s.system == sy]
    m = length(ss)
    # weak quadratic penalty on the unbounded scale keeps the logistic box maps from saturating (a parameter at a
    # bound has zero gradient and cannot move back); PRIOR_W = 0.01 is about 1 nat at 10 box-units from centre
    function f(θ, _)
        u = unpack(θ, m, sy)
        return sum(site_nll(ss[i], obs[ss[i].id], (beta=u.beta[i], I0=u.I0[i], mu=u.mu, eps=u.eps, a=u.a, b=u.b, sd=u.sd, thr=u.thr))
                   for i in 1:m) + PRIOR_W * sum(abs2, θ)
    end
    θ0 = init_theta(ss, sy)
    t = @elapsed sol = solve(OptimizationProblem(OptimizationFunction(f, Optimization.AutoForwardDiff()), θ0), LBFGS(); maxiters=iters)
    u = unpack(sol.u, m, sy)
    println("fit $sy: $m sites, nll $(round(f(sol.u, nothing); digits=1)), mu $(round(u.mu; sigdigits=3)), eps $(round(u.eps; sigdigits=3)), $(round(t / 60; digits=1)) min"); flush(stdout)
    return (sites=ss, u=u, nll=f(sol.u, nothing))
end

function predictions(fit)
    rows = Vector{NamedTuple}()
    for (i, s) in enumerate(fit.sites)
        u = fit.u; o = obs[s.id]
        p = EstablishmentParams(lab, pw; beta=u.beta[i], mu=u.mu, epsilon=u.eps, host=1.0)
        rec = sort(unique(x.day for x in o))
        intro = s.system == "gb" ? clamp(round(Int, (GB_ARRIVAL - year(s.start)) * 365.25) + 1, 1, length(s.T)) : 1
        r = simulate_establishment(p, s.T, s.S, s.N0, 0.0; intro_day=intro, background_mortality_d=BACKGROUND_D, record=rec,
                                   smooth=true, allow_clearance=true, record_levels=true, pcr_detection_log2=u.thr)
        at = Dict(d => k for (k, d) in enumerate(rec))
        for x in o
            k = at[x.day]
            pc = zeros(length(CLASSES)); for l in 1:EST_N_LEVELS; pc[CLASS_OF_LEVEL[l]] += r.levels[l, k]; end
            pred = x.kind == :prev ? (s.pcr ? r.pcr_prev[k] : r.prev[k]) : x.kind == :wp ? r.mackin[k] :
                   x.kind == :md_int ? u.a + u.b * r.mackin[k] : x.kind == :classes ? sum(pc .* CLASSES) :
                   (pos = sum(r.levels[l, k] for l in 2:EST_N_LEVELS); u.a + u.b * (pos > 1e-9 ? sum(r.levels[l, k] * LOG2_CELLS[l] for l in 2:EST_N_LEVELS) / pos : 0.0))
            observed = x.kind == :classes ? sum(x.counts .* CLASSES) / sum(x.counts) : x.value
            push!(rows, (site=s.id, system=s.system, day=x.day, date=s.start + Day(x.day - 1), kind=string(x.kind), n=x.n,
                         obs=observed, pred=pred, beta=u.beta[i], I0=u.I0[i],
                         pc=join(round.(pc; digits=3), ";"), oc=x.kind == :classes ? join(x.counts, ";") : ""))
        end
    end
    return rows
end

rmse(a, b) = sqrt(mean((a .- b) .^ 2))
fits = Dict{String,Any}(); fitlock = ReentrantLock()
syslist = [sy for sy in SYS if any(s.system == sy for s in sites)]
Threads.@threads for sy in syslist
    f = fit_system(sy)
    lock(fitlock) do; fits[sy] = f; end
end
pred = DataFrame(vcat([predictions(f) for f in values(fits)]...))
Arrow.write(joinpath(OUT, "obs.arrow"), pred)
sitetab = DataFrame(site=String[], system=String[], beta=Float64[], I0=Float64[])
for (sy, f) in fits, (i, s) in enumerate(f.sites)
    push!(sitetab, (s.id, sy, f.u.beta[i], f.u.I0[i]))
end
Arrow.write(joinpath(OUT, "sites.arrow"), sitetab)

fmt(x; d=3) = (x isa Real && isnan(x)) ? "n/a" : string(round(x; digits=d))
function demeaned_cor(d)
    g = groupby(d, :site)
    o = vcat([x.obs .- mean(x.obs) for x in g]...); p = vcat([x.pred .- mean(x.pred) for x in g]...)
    return length(o) > 2 ? cor(o, p) : NaN
end
open(joinpath(OUT, "fit.md"), "w") do io
    println(io, "# P4.2j all-observation fit of the staged establishment model (scripts/p4_2j_fit.jl)\n")
    println(io, "Observations used: ", join(["$(k) $(count(==(string(k)), pred.kind))" for k in (:prev, :classes, :wp, :md_int, :copies)], ", "), ".\n")
    println(io, "## Prevalence, in-sample\n\n| system | sites | n | RMSE model | RMSE persistence | RMSE system mean | within-site year-to-year correlation |\n|---|---|---|---|---|---|---|")
    for sy in SYS
        d = pred[(pred.system .== sy) .& (pred.kind .== "prev"), :]; nrow(d) == 0 && continue
        pers = Float64[]; po = Float64[]
        for s in fits[sy].sites, (i, v) in enumerate(s.prev_persist)
            ismissing(v) || (push!(pers, v); push!(po, s.prev_obs[i]))
        end
        println(io, "| $sy | $(length(unique(d.site))) | $(nrow(d)) | $(fmt(rmse(d.pred, d.obs))) | $(isempty(pers) ? "n/a" : fmt(rmse(pers, po))) | $(fmt(rmse(fill(mean(d.obs), nrow(d)), d.obs))) | $(fmt(demeaned_cor(d))) |")
    end
    println(io, "\n## Intensity and class fits, in-sample\n\n| system | kind | n | RMSE | correlation | within-site correlation |\n|---|---|---|---|---|---|")
    for sy in SYS, k in ("classes", "wp", "md_int", "copies")
        d = pred[(pred.system .== sy) .& (pred.kind .== k), :]; nrow(d) == 0 && continue
        println(io, "| $sy | $(k == "classes" ? "mean class (0-5)" : k) | $(nrow(d)) | $(fmt(rmse(d.pred, d.obs))) | $(fmt(cor(d.pred, d.obs))) | $(fmt(demeaned_cor(d))) |")
    end
    println(io, "\n## Class distributions (pooled observed vs predicted shares)\n\n| system | class | observed | predicted |\n|---|---|---|---|")
    for sy in ("gb", "mw")
        d = pred[(pred.system .== sy) .& (pred.kind .== "classes"), :]; nrow(d) == 0 && continue
        oc = sum(parse.(Float64, split(x, ";")) for x in d.oc); pc = sum(parse.(Float64, split(x, ";")) .* parse.(Float64, split(y, ";")) for (x, y) in zip(d.pc, d.oc) for _ in 1:1) # weighted by n
        pcw = sum(parse.(Float64, split(x, ";")) .* sum(parse.(Float64, split(y, ";"))) for (x, y) in zip(d.pc, d.oc))
        for (j, cl) in enumerate(CLASSES)
            println(io, "| $sy | $cl | $(fmt(oc[j] / sum(oc))) | $(fmt(pcw[j] / sum(pcw))) |")
        end
    end
    println(io, "\n## Parameters\n\n| system | mu (d^-1) | epsilon (d^-1) | beta q10 / median / q90 | beta at bounds | observation params |\n|---|---|---|---|---|---|")
    for sy in SYS
        haskey(fits, sy) || continue
        u = fits[sy].u; q = quantile(u.beta, [0.1, 0.5, 0.9])
        atb = count(b -> b <= 1.01 * BOX.beta[1] || b >= 0.99 * BOX.beta[2], u.beta)
        extra = sy == "md" ? "0-7 = $(round(u.a; digits=2)) + $(round(u.b; digits=2)) x rank, sd $(round(u.sd; digits=2))" :
                sy == "ne" ? "WP sd $(round(u.sd; digits=2))" :
                sy == "me" ? "qPCR threshold 2^$(round(u.thr; digits=1)) cells; log10 copies = $(round(u.a; digits=2)) + $(round(u.b; digits=3)) x log2 cells, sd $(round(u.sd; digits=2))" : ""
        println(io, "| $sy | $(round(u.mu; sigdigits=3)) | $(round(u.eps; sigdigits=3)) | $(round(q[1]; sigdigits=2)) / $(round(q[2]; sigdigits=2)) / $(round(q[3]; sigdigits=2)) | $atb of $(length(u.beta)) | $extra |")
    end
end
println(read(joinpath(OUT, "fit.md"), String))
