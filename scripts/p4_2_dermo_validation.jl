# P4.2: test and calibrate the Ragone Calvo et al. (2001) Dermo model on Maryland DNR disease bars.
#
# Design (fixed before running):
# - One-year-ahead test: for each bar and each year y with an observed fall mean intensity I_y and I_{y+1}
#   and continuous forcing between the two surveys, start from P0 = cells(I_y) (Bushek et al. 1994 via
#   rank_to_log10_cells) on 1 November of y, run 365 days, and compare cells_to_intensity(P) with I_{y+1}.
# - Baseline: persistence (predict I_{y+1} = I_y).
# - Calibration: one multiplier k on the parasite division rate (mu_d -> k * mu_d), chosen by grid search to
#   minimise squared error on the calibration bars; skill reported on held-out bars (whole bars held out).
#   Calibration bars: Hog Island, Broome Island, Swan Point, Sandy Hill.
#   Held-out bars: Hackett Point, Oyster Shell Pt, Sharkfin Shoal, Goose Creek.
# Fall survey rows are dated 1 November (year resolution only; see PLAN.md P1.7).
# Only survey years from 1999 are used: before 1999 MD DNR used less sensitive hemolymph assays instead of
# rectal-tissue RFTM (2024 Fall Survey report, p. 11), so earlier intensities are not comparable.
# MD DNR mean intensity includes zeros for uninfected oysters (report p. 11), i.e. the "weighted prevalence"
# that Ragone Calvo et al. (2001, p. 234) convert to mean parasite density.
# The fitted multiplier is written to params/ only if it is inside the grid (an edge value is not an estimate).
#
# Usage: julia --project scripts/p4_2_dermo_validation.jl
# Writes outputs/p4_2/{validation.md, predictions.arrow} and params/dermo_division_multiplier.toml.

using OysterFarming
using Arrow, DataFrames, Dates, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "outputs", "p4_2")
const CALIB = ["mddnr_hog_island", "mddnr_broome_island", "mddnr_swan_point", "mddnr_sandy_hill"]
const HOLDOUT = ["mddnr_hackett_point", "mddnr_oyster_shell_pt", "mddnr_sharkfin_shoal", "mddnr_goose_creek"]
const K_GRID = 0.1:0.05:3.0
const FIRST_YEAR = 1999

base = DermoParams(load_params(joinpath(ROOT, "params", "dermo_ragone_calvo2001.toml")))
forcing = read_contract(joinpath(ROOT, "data", "processed", "forcing_daily.arrow"), "forcing_daily")
disease = read_contract(joinpath(ROOT, "data", "processed", "disease_obs.arrow"), "disease_obs")
dermo = disease[(disease.pathogen .== "P_marinus") .& .!ismissing.(disease.intensity_mean), :]

scaled(p::DermoParams, k) = DermoParams((getfield(p, f) * (f in (:div_intercept, :div_slope) ? k : 1)
                                         for f in fieldnames(DermoParams))...)

"One-year-ahead cases for a bar: (year, I_y, I_next, forcing slice)."
function cases(site)
    f = forcing_for_site(forcing, site)
    obs = sort(dermo[dermo.site_id .== site, :], :date)
    out = []
    for i in 1:nrow(obs)-1
        y0, y1 = obs.date[i], obs.date[i+1]
        y1 == y0 + Year(1) || continue
        year(y0) >= FIRST_YEAR || continue
        sl = f[(f.date .> y0) .& (f.date .<= y1), :]
        (nrow(sl) == Dates.value(y1 - y0) && all(diff(sl.date) .== Day(1))) || continue
        push!(out, (year=year(y1), I0=Float64(obs.intensity_mean[i]), I1=Float64(obs.intensity_mean[i+1]), forcing=sl))
    end
    return out
end

function predict(p, c)
    P0 = c.I0 > 0 ? 10.0^rank_to_log10_cells(c.I0) : 0.0
    P = simulate_dermo(p, c.forcing, P0).P_cells_g[end]
    return cells_to_intensity(P)
end

allcases = Dict(s => cases(s) for s in vcat(CALIB, HOLDOUT))
sse(p, sites) = sum((predict(p, c) - c.I1)^2 for s in sites for c in allcases[s])

errs = [(k, sse(scaled(base, k), CALIB)) for k in K_GRID]
k_best = errs[argmin(last.(errs))][1]
calibrated = scaled(base, k_best)

rows = DataFrame(site_id=String[], role=String[], year=Int[], I_prev=Float64[], I_obs=Float64[],
                 I_published=Float64[], I_calibrated=Float64[])
for s in vcat(CALIB, HOLDOUT), c in allcases[s]
    push!(rows, (s, s in CALIB ? "calibration" : "held out", c.year, c.I0, c.I1, predict(base, c), predict(calibrated, c)))
end

rmse(a, b) = sqrt(mean((a .- b) .^ 2))
ranks(x) = invperm(sortperm(x))
spearman(a, b) = cor(ranks(a), ranks(b))
log10cells(I) = I > 0 ? rank_to_log10_cells(I) : 0.0
at_edge = k_best in (first(K_GRID), last(K_GRID))
mkpath(OUT)
Arrow.write(joinpath(OUT, "predictions.arrow"), rows)
fmt(x) = string(round(x; digits=3))
open(joinpath(OUT, "validation.md"), "w") do io
    println(io, "# P4.2 Dermo model one-year-ahead test (scripts/p4_2_dermo_validation.jl)\n")
    println(io, "Division-rate multiplier fitted on calibration bars: k = $(k_best) (grid $(first(K_GRID))-$(last(K_GRID))).\n")
    at_edge && println(io, "**The fit is at the edge of the grid, so k is not an estimate and no parameter file is written.**
")
    println(io, "Survey years $(FIRST_YEAR) onward. Intensity is the MD DNR 0-7 mean over all oysters.
")
    println(io, "| bars | n | RMSE persistence | RMSE published | RMSE calibrated | bias published | bias calibrated | Spearman published | Spearman calibrated | RMSE log10 cells published |")
    println(io, "|---|---|---|---|---|---|---|---|---|---|")
    for (label, sel) in (("calibration", rows.role .== "calibration"), ("held out", rows.role .== "held out"))
        r = rows[sel, :]
        println(io, "| $label | $(nrow(r)) | $(fmt(rmse(r.I_prev, r.I_obs))) | $(fmt(rmse(r.I_published, r.I_obs))) | ",
                "$(fmt(rmse(r.I_calibrated, r.I_obs))) | $(fmt(mean(r.I_published .- r.I_obs))) | ",
                "$(fmt(mean(r.I_calibrated .- r.I_obs))) | $(fmt(spearman(r.I_published, r.I_obs))) | ",
                "$(fmt(spearman(r.I_calibrated, r.I_obs))) | $(fmt(rmse(log10cells.(r.I_published), log10cells.(r.I_obs)))) |")
    end
    missing_bars = [s for s in vcat(CALIB, HOLDOUT) if isempty(allcases[s])]
    isempty(missing_bars) || println(io, "
Bars with no usable year pair (forcing gaps): ", join(missing_bars, ", "), ".")
    println(io, "\nPer bar (calibrated model):\n")
    println(io, "| bar | role | n | RMSE persistence | RMSE calibrated | corr |")
    println(io, "|---|---|---|---|---|---|")
    for g in groupby(rows, :site_id)
        println(io, "| $(g.site_id[1]) | $(g.role[1]) | $(nrow(g)) | $(fmt(rmse(g.I_prev, g.I_obs))) | ",
                "$(fmt(rmse(g.I_calibrated, g.I_obs))) | $(nrow(g) > 2 ? fmt(cor(g.I_calibrated, g.I_obs)) : "") |")
    end
end

at_edge || open(joinpath(ROOT, "params", "dermo_division_multiplier.toml"), "w") do io
    println(io, "title = \"Multiplier on the Ragone Calvo et al. (2001) parasite division rate, calibrated for Maryland (P4.2)\"\n")
    println(io, "[division_multiplier]")
    println(io, "value = $(k_best)")
    println(io, "units = \"1\"")
    println(io, "source = \"Calibrated in this project: scripts/p4_2_dermo_validation.jl, grid search minimising one-year-ahead squared error in mean Dermo intensity on MD DNR disease bars $(join(CALIB, ", ")) (Fall Survey Tables 3), forcing from CO-OPS and CBP stations; held-out results in outputs/p4_2/validation.md\"")
end
println(read(joinpath(OUT, "validation.md"), String))
