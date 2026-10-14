# P4.2c: statistical Dermo model fitted to Maryland DNR bars, tested on held-out bars, applied to Maine.
#
# Both published process models failed the Maryland test (P4.2, P4.2b). This fits the survey data directly.
# Design (fixed before running):
# - Bars: every MD DNR disease bar with bar forcing (42 of 43; P4.2c data expansion). Split, fixed and
#   deterministic: the 8 bars of P4.2 keep their roles; the others, ordered north to south by latitude, alternate
#   calibration / held out. Second test, leave one region out: bars grouped by latitude into upper (>= 38.8 N),
#   middle (38.45-38.8 N) and lower (< 38.45 N) bay; fit on two regions, score on the third.
# - Forcing gaps filled by repeating the last value; a bar-year is used only if fewer than 30 days of its window
#   were filled. Survey years 1999-2024 (assay change in 1999).
# - Response per bar-year y (fall survey, dated 1 Nov): mean intensity I (MD DNR 0-7 scale, zeros included) and
#   prevalence P (fitted on the logit scale, P clipped to [0.01, 0.99]).
# - Climate predictors from the 365 days ending 31 Oct of y: DD20 = sum of max(T - 20, 0) (degC days; 20 degC is
#   where P. marinus division starts in both published models), S = mean salinity 1 Jun-31 Oct, Tw = mean
#   temperature 1 Jan-31 Mar.
# - Model A (one year ahead): response_y = a + b response_{y-1} + c DD20 + d S + e Tw.
#   Model B (climate only):    response_y = a + c DD20 + d S + e Tw.
#   Ordinary least squares on the calibration bars; skill on held-out bars against persistence
#   (response_{y-1}) and against the calibration-bar mean.
# - Maine: each Damariscotta site's predictors (present and +0.5/1/2 degC delta warming) are compared with the
#   Maryland calibration range. Model B is applied only where all predictors are inside that range; otherwise the
#   prediction is reported as not applicable (extrapolation). Excess oyster mortality where applicable = 0.0902
#   yr^-1 per unit intensity (params/dermo_mortality_mddnr.toml).
#
# Usage: julia --project scripts/p4_2c_dermo_statistical.jl
# Writes outputs/p4_2c/{model.md, fits.arrow, maine.arrow}.

using OysterFarming
using Arrow, DataFrames, Dates, LinearAlgebra, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "outputs", "p4_2c")
const CALIB8 = ["mddnr_hog_island", "mddnr_broome_island", "mddnr_swan_point", "mddnr_sandy_hill"]
const HOLDOUT8 = ["mddnr_hackett_point", "mddnr_oyster_shell_pt", "mddnr_sharkfin_shoal", "mddnr_goose_creek"]
const FIRST_YEAR = 1999
const MAX_HELD_DAYS = 30
const WARMING = [0.0, 0.5, 1.0, 2.0]
const PCLIP = 0.01

forcing = read_contract(joinpath(ROOT, "data", "processed", "forcing_daily.arrow"), "forcing_daily")
disease = read_contract(joinpath(ROOT, "data", "processed", "disease_obs.arrow"), "disease_obs")
dermo = disease[disease.pathogen .== "P_marinus", :]
mort = load_params(joinpath(ROOT, "params", "dermo_mortality_mddnr.toml"))["dermo_mortality_per_intensity_yr"].value
sites = read_contract(joinpath(ROOT, "data", "processed", "sites.arrow"), "sites")
lat = Dict(r.site_id => r.lat for r in eachrow(sites))
bars = sort([b for b in unique(forcing.site_id) if startswith(b, "mddnr_")]; by=b -> -lat[b])
others = [b for b in bars if !(b in CALIB8 || b in HOLDOUT8)]
const CALIB = vcat(CALIB8, others[1:2:end])
const HOLDOUT = vcat(HOLDOUT8, others[2:2:end])
region(b) = lat[b] >= 38.8 ? "upper" : lat[b] >= 38.45 ? "middle" : "lower"

"Daily T and S over the site's whole range, gaps filled by repeating the last value; `held` marks filled days."
function filled(site)
    f = forcing_for_site(forcing, site)
    have = Dict(r.date => (Float64(r.temp_c), Float64(r.salinity_psu)) for r in eachrow(f)
                if !ismissing(r.temp_c) && !ismissing(r.salinity_psu))
    dates = collect(minimum(keys(have)):Day(1):maximum(keys(have)))
    last = have[dates[1]]
    T = Float64[]; S = Float64[]; held = Bool[]
    for d in dates
        haskey(have, d) ? (last = have[d]; push!(held, false)) : push!(held, true)
        push!(T, last[1]); push!(S, last[2])
    end
    return DataFrame(date=dates, temp_c=T, salinity_psu=S, held=held)
end

"Predictors for survey year y, or nothing if the window is not covered or too many days were filled."
function predictors(f, y; dT=0.0)
    w = (f.date .> Date(y - 1, 10, 31)) .& (f.date .<= Date(y, 10, 31))
    count(w) == Dates.value(Date(y, 10, 31) - Date(y - 1, 10, 31)) || return nothing
    count(f.held[w]) < MAX_HELD_DAYS || return nothing
    T = f.temp_c[w] .+ dT; d = f.date[w]; S = f.salinity_psu[w]
    summer = (month.(d) .>= 6) .& (month.(d) .<= 10)
    winter = month.(d) .<= 3
    return (DD20=sum(max.(T .- 20, 0.0)), S=mean(S[summer]), Tw=mean(T[winter]))
end

logit(p) = log(p / (1 - p))
expit(x) = 1 / (1 + exp(-x))

rows = DataFrame(site_id=String[], role=String[], year=Int[], I=Float64[], I_prev=Union{Missing,Float64}[],
                 P=Float64[], P_prev=Union{Missing,Float64}[], DD20=Float64[], S=Float64[], Tw=Float64[])
for s in vcat(CALIB, HOLDOUT)
    f = filled(s)
    obs = sort(dermo[(dermo.site_id .== s) .& .!ismissing.(dermo.intensity_mean), :], :date)
    byyear = Dict(year(r.date) => r for r in eachrow(obs))
    for (y, r) in byyear
        y >= FIRST_YEAR || continue
        x = predictors(f, y)
        x === nothing && continue
        prev = get(byyear, y - 1, nothing)
        push!(rows, (s, s in CALIB ? "calibration" : "held out", y, Float64(r.intensity_mean),
                     prev === nothing ? missing : Float64(prev.intensity_mean), Float64(r.prevalence),
                     prev === nothing ? missing : Float64(prev.prevalence), x.DD20, x.S, x.Tw))
    end
end
sort!(rows, [:site_id, :year])

"OLS fit: coefficients, standard errors, residual sd."
function ols(X, y)
    β = X \ y
    r = y .- X * β
    σ2 = sum(abs2, r) / (length(y) - size(X, 2))
    se = sqrt.(diag(σ2 * inv(X' * X)))
    return (β=β, se=se, sd=sqrt(σ2))
end

climX(d) = hcat(ones(nrow(d)), d.DD20, d.S, d.Tw)
tr_resp(v, kind) = kind == :I ? v : logit.(clamp.(v, PCLIP, 1 - PCLIP))
back(v, kind) = kind == :I ? clamp.(v, 0.0, 7.0) : expit.(v)
rmse(a, b) = sqrt(mean((a .- b) .^ 2))

results = Dict{Tuple{Symbol,Symbol},Any}()
cal = rows[rows.role .== "calibration", :]
hold = rows[rows.role .== "held out", :]
for kind in (:I, :P)
    resp, prevcol = kind == :I ? (:I, :I_prev) : (:P, :P_prev)
    # Model B: climate only
    fb = ols(climX(cal), tr_resp(cal[!, resp], kind))
    # Model A: lag + climate (rows with a previous year)
    ca = cal[.!ismissing.(cal[!, prevcol]), :]
    Xa(d) = hcat(climX(d)[:, 1], tr_resp(Float64.(d[!, prevcol]), kind), climX(d)[:, 2:end])
    fa = ols(Xa(ca), tr_resp(ca[!, resp], kind))
    ha = hold[.!ismissing.(hold[!, prevcol]), :]
    results[(kind, :B)] = (fit=fb, pred=back(climX(hold) * fb.β, kind), obs=hold[!, resp], n_cal=nrow(cal))
    results[(kind, :A)] = (fit=fa, pred=back(Xa(ha) * fa.β, kind), obs=ha[!, resp], pers=Float64.(ha[!, prevcol]),
                           n_cal=nrow(ca))
end

# Leave one region out (all bars)
regional = DataFrame(region=String[], response=String[], model=String[], n=Int[], rmse_model=Float64[],
                     rmse_persistence=Union{Missing,Float64}[], rmse_train_mean=Float64[], pearson=Float64[])
rows.region = region.(rows.site_id)
for reg in ("upper", "middle", "lower"), kind in (:I, :P)
    resp, prevcol = kind == :I ? (:I, :I_prev) : (:P, :P_prev)
    tr = rows[rows.region .!= reg, :]; te = rows[rows.region .== reg, :]
    fb = ols(climX(tr), tr_resp(tr[!, resp], kind))
    pb = back(climX(te) * fb.β, kind)
    push!(regional, (reg, string(kind), "B", nrow(te), rmse(pb, te[!, resp]), missing,
                     rmse(fill(mean(tr[!, resp]), nrow(te)), te[!, resp]), cor(pb, te[!, resp])))
    tra = tr[.!ismissing.(tr[!, prevcol]), :]; tea = te[.!ismissing.(te[!, prevcol]), :]
    Xr(d) = hcat(ones(nrow(d)), tr_resp(Float64.(d[!, prevcol]), kind), d.DD20, d.S, d.Tw)
    fa = ols(Xr(tra), tr_resp(tra[!, resp], kind))
    pa = back(Xr(tea) * fa.β, kind)
    push!(regional, (reg, string(kind), "A", nrow(tea), rmse(pa, tea[!, resp]), rmse(Float64.(tea[!, prevcol]), tea[!, resp]),
                     rmse(fill(mean(tra[!, resp]), nrow(tea)), tea[!, resp]), cor(pa, tea[!, resp])))
end

# Maine application, model B
ranges = Dict(v => extrema(cal[!, v]) for v in (:DD20, :S, :Tw))
maine = DataFrame(site=String[], warming_c=Float64[], years=Int[], DD20=Float64[], DD20_frac_md_min=Float64[],
                  S=Float64[], Tw=Float64[], outside_range=String[], intensity=Union{Missing,Float64}[],
                  prevalence=Union{Missing,Float64}[], excess_mortality_yr=Union{Missing,Float64}[])
for site in ("dam_lower", "dam_upper"), dT in WARMING
    f = filled(site)
    xs = [x for y in year(first(f.date))+1:year(last(f.date)) for x in (predictors(f, y; dT=dT),) if x !== nothing]
    isempty(xs) && continue
    x = (DD20=mean(getfield.(xs, :DD20)), S=mean(getfield.(xs, :S)), Tw=mean(getfield.(xs, :Tw)))
    out = join([string(v) for v in (:DD20, :S, :Tw) if !(ranges[v][1] <= getfield(x, v) <= ranges[v][2])], ", ")
    xv = [1.0, x.DD20, x.S, x.Tw]
    ok = isempty(out)
    I = ok ? back([dot(xv, results[(:I, :B)].fit.β)], :I)[1] : missing
    P = ok ? back([dot(xv, results[(:P, :B)].fit.β)], :P)[1] : missing
    push!(maine, (site, dT, length(xs), x.DD20, x.DD20 / ranges[:DD20][1], x.S, x.Tw, ok ? "none" : out, I, P,
                  ok ? mort * I : missing))
end

mkpath(OUT)
Arrow.write(joinpath(OUT, "fits.arrow"), rows)
Arrow.write(joinpath(OUT, "maine.arrow"), maine)
fmt(x; d=3) = string(round(x; digits=d))
names_B = ["intercept", "DD20", "S", "Tw"]
names_A = ["intercept", "lag", "DD20", "S", "Tw"]
open(joinpath(OUT, "model.md"), "w") do io
    println(io, "# P4.2c statistical Dermo model (generated by scripts/p4_2c_dermo_statistical.jl)\n")
    println(io, "Bars with usable years: $(length(unique(rows.site_id))). Calibration: $(length(unique(cal.site_id))) bars ($(nrow(cal)) bar-years). ",
            "Held out: $(length(unique(hold.site_id))) bars ($(nrow(hold)) bar-years). Survey years $(FIRST_YEAR)-2024.\n")
    println(io, "Calibration bars: $(join(sort(unique(cal.site_id)), ", ")).\n\nHeld-out bars: $(join(sort(unique(hold.site_id)), ", ")).\n")
    println(io, "## Held-out skill\n")
    println(io, "| response | model | n | RMSE model | RMSE persistence | RMSE calibration mean | Pearson model |")
    println(io, "|---|---|---|---|---|---|---|")
    for kind in (:I, :P), m in (:A, :B)
        r = results[(kind, m)]
        calmean = mean(cal[!, kind == :I ? :I : :P])
        pers = m == :A ? fmt(rmse(r.pers, r.obs)) : "n/a"
        println(io, "| $(kind == :I ? "intensity (0-7)" : "prevalence") | $(m == :A ? "A: lag + climate" : "B: climate only") | ",
                "$(length(r.obs)) | $(fmt(rmse(r.pred, r.obs))) | $pers | $(fmt(rmse(fill(calmean, length(r.obs)), r.obs))) | ",
                "$(fmt(cor(r.pred, r.obs))) |")
    end
    println(io, "\n## Leave one region out (fit on the other two regions)\n")
    println(io, "| held-out region | response | model | n | RMSE model | RMSE persistence | RMSE training mean | Pearson |")
    println(io, "|---|---|---|---|---|---|---|---|")
    for r in eachrow(regional)
        println(io, "| $(r.region) | $(r.response == "I" ? "intensity" : "prevalence") | $(r.model) | $(r.n) | $(fmt(r.rmse_model)) | ",
                "$(ismissing(r.rmse_persistence) ? "n/a" : fmt(r.rmse_persistence)) | $(fmt(r.rmse_train_mean)) | $(fmt(r.pearson)) |")
    end
    println(io, "\n## Coefficients (calibration bars; prevalence on the logit scale)\n")
    for kind in (:I, :P), m in (:A, :B)
        r = results[(kind, m)]
        nm = m == :A ? names_A : names_B
        println(io, "$(kind == :I ? "Intensity" : "Prevalence"), model $(m) (n = $(r.n_cal), residual sd $(fmt(r.fit.sd))): ",
                join(["$(nm[i]) $(round(r.fit.β[i]; sigdigits=3)) (se $(round(r.fit.se[i]; sigdigits=2)))" for i in eachindex(nm)], "; "), "\n")
    end
    println(io, "Calibration predictor ranges: ", join(["$(v) $(fmt(ranges[v][1]; d=1))-$(fmt(ranges[v][2]; d=1))" for v in (:DD20, :S, :Tw)], "; "), ".\n")
    println(io, "## Maine (mean over each site's forcing years)\n")
    println(io, "| site | warming degC | years | DD20 | DD20 / Maryland minimum | summer S | winter T | outside MD range | intensity (0-7) | prevalence | excess mortality yr^-1 |")
    println(io, "|---|---|---|---|---|---|---|---|---|---|---|")
    na(x; d=3) = ismissing(x) ? "n/a (extrapolation)" : fmt(x; d=d)
    for r in eachrow(maine)
        println(io, "| $(r.site) | $(r.warming_c) | $(r.years) | $(fmt(r.DD20; d=1)) | $(fmt(r.DD20_frac_md_min; d=2)) | ",
                "$(fmt(r.S; d=1)) | $(fmt(r.Tw; d=1)) | $(r.outside_range) | $(na(r.intensity)) | $(na(r.prevalence)) | ",
                "$(na(r.excess_mortality_yr; d=4)) |")
    end
    println(io, "\nEvery Maine case is outside the Maryland data: summer salinity about 30 (Maryland at most ",
            "$(fmt(ranges[:S][2]; d=1))) and far fewer degree-days above 20 degC. The fitted model cannot be used for ",
            "Maine. What the data do say: even with 2 degC warming, Damariscotta summers stay at or below the coldest ",
            "Maryland bar-year in the calibration set (DD20 / Maryland minimum <= 1).")
end
println(read(joinpath(OUT, "model.md"), String))
