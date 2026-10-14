# P4.2d: out-of-region test of Dermo models: fit on Maryland DNR bars, predict Great Bay, New Hampshire.
#
# Great Bay is the only data found with Maine-like conditions (Great Bay buoy: 37-169 degree-days above 20 degC,
# summer salinity 21-30). A model that cannot predict Great Bay should not be used for the Damariscotta.
# Design (fixed before running):
# - Response: fall Dermo prevalence (RFTM in both programs; intensity scales differ, 0-7 vs 0-5, so not used).
#   Maryland: fall survey (dated 1 Nov), years 1999-2024. Great Bay: last survey of each bed-year from August on.
# - Predictors (summer only, because Great Bay sondes are removed in winter): DD20 = sum of max(T - 20, 0) over
#   1 Jun-31 Oct; S = mean salinity 1 Jun-31 Oct. A bed-year is used only if at least 90 % of those 153 days
#   have forcing. Maryland gaps are as filled in forcing_daily; Great Bay forcing has measured fills only.
# - Candidates, all fitted on Maryland only:
#   A  logit P_y = a + b logit P_{y-1} + c DD20 + d S   (needs last year's survey at the same bed)
#   B  logit P_y = a + c DD20 + d S                      (climate only; the form usable for Maine)
#   PK Hofmann/Powell infection-class model (P4.2b, 60 h halving time) run continuously on each Great Bay bed
#      from 1996-01-01 with the P4.2b settings; winter days without sonde data repeat the last value (in that
#      model every temperature below 10 degC gives the same parasite mortality, and division is negligible).
# - Baselines: persistence (last year's observed prevalence at the bed) and the Great Bay mean (which uses the test
#   data, so it is a generous baseline).
#
# Usage: julia --project scripts/p4_2d_greatbay_test.jl
# Writes outputs/p4_2d/greatbay.md and .arrow.

using OysterFarming
using Arrow, DataFrames, Dates, LinearAlgebra, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "outputs", "p4_2d")
const FIRST_MD_YEAR = 1999
const MIN_COVER = 0.9
const PCLIP = 0.01
const BACKGROUND_YR = 0.0302          # as P4.2b

forcing = read_contract(joinpath(ROOT, "data", "processed", "forcing_daily.arrow"), "forcing_daily")
disease = read_contract(joinpath(ROOT, "data", "processed", "disease_obs.arrow"), "disease_obs")
dermo = disease[disease.pathogen .== "P_marinus", :]

logit(p) = log(p / (1 - p))
expit(x) = 1 / (1 + exp(-x))
lg(p) = logit(clamp(p, PCLIP, 1 - PCLIP))
rmse(a, b) = sqrt(mean((a .- b) .^ 2))

"Summer predictors for year y from a site's forcing rows, or nothing if coverage is too low."
function summer(f, y)
    w = (f.date .>= Date(y, 6, 1)) .& (f.date .<= Date(y, 10, 31))
    count(w) >= MIN_COVER * 153 || return nothing
    return (DD20=sum(max.(Float64.(f.temp_c[w]) .- 20, 0.0)) * 153 / count(w), S=mean(Float64.(f.salinity_psu[w])))
end

function site_years(prefix, first_year; month_from=1)
    rows = DataFrame(site_id=String[], year=Int[], P=Float64[], P_prev=Union{Missing,Float64}[], DD20=Float64[], S=Float64[])
    have_forcing = Set(forcing.site_id)
    for s in unique(dermo.site_id[startswith.(dermo.site_id, prefix)])
        s in have_forcing || continue
        f = forcing_for_site(forcing, s)
        obs = sort(dermo[(dermo.site_id .== s) .& (month.(dermo.date) .>= month_from), :], :date)
        last = Dict{Int,Float64}()
        for r in eachrow(obs)
            last[year(r.date)] = Float64(r.prevalence)       # last survey of the year wins
        end
        for (y, p) in sort(collect(last))
            y >= first_year || continue
            x = summer(f, y)
            x === nothing && continue
            push!(rows, (s, y, p, get(last, y - 1, missing), x.DD20, x.S))
        end
    end
    return rows
end

md = site_years("mddnr_", FIRST_MD_YEAR)
gb = site_years("nh_", 1996; month_from=8)

XB(d) = hcat(ones(nrow(d)), d.DD20, d.S)
XA(d) = hcat(ones(nrow(d)), lg.(Float64.(d.P_prev)), d.DD20, d.S)
mdA = md[.!ismissing.(md.P_prev), :]
βB = XB(md) \ lg.(md.P)
βA = XA(mdA) \ lg.(mdA.P)
gbA = gb[.!ismissing.(gb.P_prev), :]
predB = expit.(XB(gb) * βB)
predA = expit.(XA(gbA) * βA)

# Hofmann/Powell model, continuous run per Great Bay bed
ptoml = load_params(joinpath(ROOT, "params", "dermo_powell1996.toml"))
pk = PowellDermoParams(ptoml)
function powell_bed(s)
    f = sort(forcing_for_site(forcing, s), :date)
    have = Dict(r.date => (Float64(r.temp_c), Float64(r.salinity_psu)) for r in eachrow(f))
    dates = collect(Date(1996, 1, 1):Day(1):Date(2011, 12, 31))
    firstv = have[minimum(keys(have))]
    last = firstv; T = Float64[]; S = Float64[]
    for d in dates
        haskey(have, d) && (last = have[d])
        push!(T, last[1]); push!(S, last[2])
    end
    N0 = zeros(length(pk.W), POWELL_N_LEVELS); N0[7, 1] = 10.0; N0[7, 2] = 10.0
    sim, _ = simulate_powell_dermo(pk, DataFrame(date=dates, temp_c=T, salinity_psu=S), N0;
                                   background_mortality_d=BACKGROUND_YR / 365)
    return Dict(d => p for (d, p) in zip(sim.date, sim.prevalence))
end
pkpred = Float64[]; pkobs = Float64[]
for s in unique(gb.site_id)
    sim = powell_bed(s)
    obs = dermo[(dermo.site_id .== s) .& (month.(dermo.date) .>= 8) .& (year.(dermo.date) .<= 2011) .& (year.(dermo.date) .>= 1998), :]
    for r in eachrow(obs)
        push!(pkpred, sim[r.date]); push!(pkobs, Float64(r.prevalence))
    end
end

mkpath(OUT)
res = transform(gb, [:DD20, :S] => ByRow((a, b) -> 0.0) => :dummy)
res.pred_B = predB
Arrow.write(joinpath(OUT, "greatbay.arrow"), select(res, Not(:dummy)))
fmt(x) = string(round(x; digits=3))
open(joinpath(OUT, "greatbay.md"), "w") do io
    println(io, "# P4.2d Dermo models fitted on Maryland, tested on Great Bay NH (scripts/p4_2d_greatbay_test.jl)\n")
    println(io, "Maryland: $(nrow(md)) bar-years ($(length(unique(md.site_id))) bars). Great Bay: $(nrow(gb)) bed-years ",
            "($(length(unique(gb.site_id))) beds), years $(minimum(gb.year))-$(maximum(gb.year)).\n")
    println(io, "Summer predictor ranges: Maryland DD20 $(fmt(minimum(md.DD20)))-$(fmt(maximum(md.DD20))), S $(fmt(minimum(md.S)))-$(fmt(maximum(md.S))); ",
            "Great Bay DD20 $(fmt(minimum(gb.DD20)))-$(fmt(maximum(gb.DD20))), S $(fmt(minimum(gb.S)))-$(fmt(maximum(gb.S))).\n")
    println(io, "| model | n | RMSE | bias (pred - obs) | Pearson | RMSE persistence (same rows) | RMSE Great Bay mean (same rows) |")
    println(io, "|---|---|---|---|---|---|---|")
    gbm(o) = rmse(fill(mean(gb.P), length(o)), o)
    println(io, "| A: lag + DD20 + S | $(nrow(gbA)) | $(fmt(rmse(predA, gbA.P))) | $(fmt(mean(predA .- gbA.P))) | $(fmt(cor(predA, gbA.P))) | ",
            "$(fmt(rmse(Float64.(gbA.P_prev), gbA.P))) | $(fmt(gbm(gbA.P))) |")
    println(io, "| B: DD20 + S (climate only) | $(nrow(gb)) | $(fmt(rmse(predB, gb.P))) | $(fmt(mean(predB .- gb.P))) | $(fmt(cor(predB, gb.P))) | n/a | $(fmt(gbm(gb.P))) |")
    println(io, "| Hofmann/Powell (60 h) | $(length(pkobs)) | $(fmt(rmse(pkpred, pkobs))) | $(fmt(mean(pkpred .- pkobs))) | ",
            "$(length(pkobs) > 2 ? fmt(cor(pkpred, pkobs)) : "n/a") | n/a | $(fmt(rmse(fill(mean(pkobs), length(pkobs)), pkobs))) |")
    println(io, "\nModel B coefficients (logit scale): intercept $(round(βB[1]; sigdigits=3)), DD20 $(round(βB[2]; sigdigits=3)), S $(round(βB[3]; sigdigits=3)).")
    println(io, "Model A coefficients: intercept $(round(βA[1]; sigdigits=3)), lag $(round(βA[2]; sigdigits=3)), DD20 $(round(βA[3]; sigdigits=3)), S $(round(βA[4]; sigdigits=3)).\n")
    println(io, "Great Bay bed-years (observed vs model B):\n")
    println(io, "| bed | year | DD20 | S | observed | last year | model B |")
    println(io, "|---|---|---|---|---|---|---|")
    for (i, r) in enumerate(eachrow(gb))
        println(io, "| $(r.site_id) | $(r.year) | $(round(r.DD20; digits=0)) | $(round(r.S; digits=1)) | $(fmt(r.P)) | ",
                "$(ismissing(r.P_prev) ? "" : fmt(r.P_prev)) | $(fmt(predB[i])) |")
    end
end
println(read(joinpath(OUT, "greatbay.md"), String))
