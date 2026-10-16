# P4.2h step 1: where do the mechanistic model's errors sit? Residual diagnosis before placing any learned correction.
#
# Design (fixed before running):
# - Model: establishment model with NO learned term, P4.2e full-fit parameters (log beta -21.295, log mu -0.69, logit host
#   6.0, log epsilon MD -8.478, log epsilon GB -6.832); New England and Maine use the Great Bay import (cool, salty
#   systems; no fitted value exists for them). Great Bay import starts at the P4.2e arrival 2003.15. Smooth mode.
# - Sites, forcing, start states and observations: scripts/p4_2g_data.jl (57 sites, 4 systems).
# - Residual = observed - predicted prevalence per observation (RFTM-detectable prevalence; model true prevalence for
#   the Maine qPCR sites).
# - Candidate inputs per observation, all computable at prediction time:
#     T_day, S_day: same-day temperature and salinity;
#     DD20_150, DD15_150: degree-days above 20 / 15 degC over the prior 150 days;
#     winter_min_T: minimum daily T over the most recent complete 1 Dec-31 Mar window before the observation;
#     prev_obs: previous observation at the site;
#     pressure: mean prevalence in the previous calendar year of OTHER sites in the same system within 30 km
#       (system-wide previous-year mean where no site is within 30 km; column pressure_local says which);
#     years_since_first_positive: years since the site's first positive observation (NaN before it);
#     month: calendar month of the observation;
#     fill_frac: fraction of the prior 150 days with an observed (qc_flag 0) forcing row (NaN for literature sites).
# - Report: Pearson r and Spearman rho of the residual with each candidate, by system and pooled; residual mean and sd;
#   base-model RMSE vs persistence. Selection rule: a candidate enters the correction if |Spearman rho| >= 0.2 in at
#   least two systems.
#
# Usage: julia --project scripts/p4_2h_residuals.jl
# Writes outputs/p4_2h/{residuals.md, residuals.arrow}.

using OysterFarming
using Arrow, DataFrames, Dates, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "outputs", "p4_2h")
include(joinpath(@__DIR__, "p4_2g_data.jl"))

const THETA = (beta=exp(-21.295), mu=exp(-0.69), host=1 / (1 + exp(-6.0)), eps_md=exp(-8.478), eps_gb=exp(-6.832))
const GB_ARRIVAL = 2003.15
const WINDOW = 150
const NEAR_KM = 30.0

# site coordinates: Maryland / Great Bay / Maine from sites.arrow; New England from the paper's site descriptions
# (Ford & Smolowitz 2007 Table 1), approximate to a few km
coords = Dict(r.site_id => (r.lat, r.lon) for r in eachrow(sites_tbl))
const NE_COORDS = Dict("ne_wellfleet_ma" => (41.93, -70.03), "ne_cotuit_ma" => (41.61, -70.44),
                       "ne_martha_s_vineyard_ma_butler_s_neck_" => (41.36, -70.55), "ne_oyster_bay_ny" => (40.88, -73.53),
                       "ne_milford_ct" => (41.21, -73.05), "ne_westport_ct" => (41.12, -73.37),
                       "ne_great_sound_nj" => (39.07, -74.78), "ne_delaware_bay_nj" => (39.12, -75.05))
function site_coord(s::Site)
    s.system == "ne" && return NE_COORDS[s.id]
    s.system == "me" && return coords[replace(s.group, "me_" => "")]
    return coords[s.id]
end
function km(a, b)
    la1, lo1, la2, lo2 = deg2rad(a[1]), deg2rad(a[2]), deg2rad(b[1]), deg2rad(b[2])
    h = sin((la2 - la1) / 2)^2 + cos(la1) * cos(la2) * sin((lo2 - lo1) / 2)^2
    return 2 * 6371 * asin(sqrt(h))
end
forcing_id(s::Site) = s.system == "ne" ? nothing : s.system == "me" ? replace(s.group, "me_" => "") : s.id

function predict_base(s::Site)
    eps = s.system == "md" ? THETA.eps_md : THETA.eps_gb
    p = EstablishmentParams(lab, pw; beta=THETA.beta, mu=THETA.mu, epsilon=eps, host=THETA.host)
    intro = s.system == "gb" ? clamp(round(Int, (GB_ARRIVAL - year(s.start)) * 365.25) + 1, 1, length(s.T)) : 1
    r = simulate_establishment(p, s.T, s.S, s.N0, 0.0; intro_day=intro, background_mortality_d=BACKGROUND_D,
                               record=s.prev_rec, smooth=true)
    return s.pcr ? r.true_prev : r.prev
end

# observation table with dates
obsdate(s::Site, d) = s.start + Day(d - 1)
sysobs = Dict(sy => [(s, obsdate(s, d), v) for s in sites if s.system == sy for (d, v) in zip(s.prev_rec, s.prev_obs)] for sy in SYSTEMS)
qc0 = Dict{String,Set{Date}}()
for r in eachrow(forcing)
    r.qc_flag == 0 && push!(get!(qc0, r.site_id, Set{Date}()), r.date)
end

rows = DataFrame(site=String[], system=String[], date=Date[], obs=Float64[], pred=Float64[], resid=Float64[],
                 T_day=Float64[], S_day=Float64[], DD20_150=Float64[], DD15_150=Float64[], winter_min_T=Float64[],
                 prev_obs=Float64[], pressure=Float64[], pressure_local=Bool[], years_since_first_positive=Float64[],
                 month=Float64[], fill_frac=Float64[])
for s in sites
    pred = predict_base(s)
    firstpos = findfirst(>(0), s.prev_obs)
    firstpos_date = firstpos === nothing ? nothing : obsdate(s, s.prev_rec[firstpos])
    here = site_coord(s)
    others = [(o, km(here, site_coord(o))) for o in sites if o.system == s.system && o.id != s.id]
    near = [o for (o, d) in others if d <= NEAR_KM]
    fid = forcing_id(s)
    for (i, d) in enumerate(s.prev_rec)
        dt = obsdate(s, d)
        lo = max(1, d - WINDOW + 1)
        Tw = view(s.T, lo:d)
        # most recent complete Dec-Mar window before dt
        wy = month(dt) >= 4 ? year(dt) : year(dt) - 1
        w1, w2 = Date(wy - 1, 12, 1), Date(wy, 3, 31)
        wmin = (w1 >= s.start && w2 < dt) ? minimum(view(s.T, dayidx(w1, s.start):dayidx(w2, s.start))) : NaN
        pool = isempty(near) ? [o for (o, _) in others] : near
        py = year(dt) - 1
        pv = [v for o in pool for (dd, v) in zip(o.prev_rec, o.prev_obs) if year(obsdate(o, dd)) == py]
        pressure = isempty(pv) ? NaN : mean(pv)
        ysfp = firstpos_date === nothing || dt < firstpos_date ? NaN : Dates.value(dt - firstpos_date) / 365.25
        ff = fid === nothing ? NaN : count(x -> x in qc0[fid], (dt - Day(WINDOW - 1)):Day(1):dt) / WINDOW
        push!(rows, (s.id, s.system, dt, s.prev_obs[i], pred[i], s.prev_obs[i] - pred[i], s.T[d], s.S[d],
                     sum(max.(Tw .- 20, 0)), sum(max.(Tw .- 15, 0)), wmin,
                     ismissing(s.prev_persist[i]) ? NaN : s.prev_persist[i], pressure, !isempty(near), ysfp, Float64(month(dt)), ff))
    end
end

function tiedrank(x)
    o = sortperm(x); r = zeros(length(x)); i = 1
    while i <= length(x)
        j = i
        while j < length(x) && x[o[j+1]] == x[o[i]]; j += 1; end
        r[o[i:j]] .= (i + j) / 2; i = j + 1
    end
    return r
end
spearman(a, b) = cor(tiedrank(a), tiedrank(b))
rmse(a, b) = sqrt(mean((a .- b) .^ 2))
const CANDS = [:T_day, :S_day, :DD20_150, :DD15_150, :winter_min_T, :prev_obs, :pressure, :years_since_first_positive, :month, :fill_frac]

mkpath(OUT)
Arrow.write(joinpath(OUT, "residuals.arrow"), rows)
fmt(x) = isnan(x) ? "n/a" : string(round(x; digits=2))
open(joinpath(OUT, "residuals.md"), "w") do io
    println(io, "# P4.2h residual diagnosis of the mechanistic Dermo model (scripts/p4_2h_residuals.jl)\n")
    println(io, "Model: establishment model, no learned term, P4.2e parameters (New England and Maine use the Great Bay import). ",
            "Residual = observed - predicted prevalence.\n")
    println(io, "| system | sites | n | mean residual | sd residual | RMSE model | RMSE persistence | RMSE system mean |\n|---|---|---|---|---|---|---|---|")
    for sy in SYSTEMS
        d = rows[rows.system .== sy, :]; nrow(d) == 0 && continue
        dp = d[.!isnan.(d.prev_obs), :]
        println(io, "| $sy | $(length(unique(d.site))) | $(nrow(d)) | $(fmt(mean(d.resid))) | $(fmt(std(d.resid))) | $(fmt(rmse(d.pred, d.obs))) | ",
                "$(fmt(rmse(dp.prev_obs, dp.obs))) | $(fmt(rmse(fill(mean(d.obs), nrow(d)), d.obs))) |")
    end
    println(io, "\n## Correlation of the residual with each candidate input (Pearson r / Spearman rho; n in brackets)\n")
    println(io, "| candidate | " * join(SYSTEMS, " | ") * " | pooled | systems with abs(rho) >= 0.2 |")
    println(io, "|---|" * repeat("---|", length(SYSTEMS) + 2))
    selected = Symbol[]
    for c in CANDS
        cells = String[]; hits = 0
        for sy in vcat(SYSTEMS, ["pooled"])
            d = sy == "pooled" ? rows : rows[rows.system .== sy, :]
            x = d[!, c]; ok = .!isnan.(x)
            if count(ok) < 5 || std(x[ok]) == 0
                push!(cells, "n/a"); continue
            end
            r = cor(x[ok], d.resid[ok]); rho = spearman(x[ok], d.resid[ok])
            sy != "pooled" && abs(rho) >= 0.2 && (hits += 1)
            push!(cells, "$(fmt(r)) / $(fmt(rho)) ($(count(ok)))")
        end
        hits >= 2 && push!(selected, c)
        println(io, "| $c | " * join(cells, " | ") * " | $hits |")
    end
    println(io, "\nSelected for the correction (|rho| >= 0.2 in at least two systems): ", isempty(selected) ? "none" : join(string.(selected), ", "), ".\n")
    println(io, "Notes: pressure uses sites within $(Int(NEAR_KM)) km where any exist (pressure_local), else the system-wide previous-year mean; ",
            "winter_min_T needs a complete Dec-Mar window inside the site's forcing; fill_frac is NaN for literature sites.")
end
println(read(joinpath(OUT, "residuals.md"), String))
