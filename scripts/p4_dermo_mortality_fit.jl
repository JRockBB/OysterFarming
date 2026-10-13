# P4.1 (oyster mortality from Dermo): fit the excess annual mortality rate per unit of Dermo intensity
# from Maryland DNR Fall Survey data, and write params/dermo_mortality_mddnr.toml.
#
# Data: data/processed/disease_obs.arrow (Table 3 Dermo mean intensity, 0-7 scale; Table 4 MSX prevalence)
#       data/processed/mortality_obs.arrow (Table 5 Observed Mortality Index, annual)
#       43 bars, 1990-2024, same fall survey for intensity and mortality.
# Model: annual instantaneous mortality Z = -ln(1 - m) = z0 + z_dermo * I + z_msx * MSX, ordinary least squares.
# MSX is included because it also kills oysters; leaving it out would bias z_dermo upward.
# Standard errors are ordinary OLS errors and ignore clustering by bar and year (they understate uncertainty);
# a bar-level block bootstrap interval is reported as well.
#
# Usage: julia --project scripts/p4_dermo_mortality_fit.jl
# Writes outputs/p4_dermo_mortality/fit.md and params/dermo_mortality_mddnr.toml.

using OysterFarming
using DataFrames, Dates, LinearAlgebra, Random, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "outputs", "p4_dermo_mortality")
const N_BOOT = 2000

d = read_contract(joinpath(ROOT, "data", "processed", "disease_obs.arrow"), "disease_obs")
m = read_contract(joinpath(ROOT, "data", "processed", "mortality_obs.arrow"), "mortality_obs")
dermo = rename(d[d.pathogen .== "P_marinus", [:site_id, :date, :intensity_mean]], :intensity_mean => :I)
msx = rename(d[d.pathogen .== "H_nelsoni", [:site_id, :date, :prevalence]], :prevalence => :msx)
j = innerjoin(innerjoin(m[:, [:site_id, :date, :mortality_fraction]], dermo, on=[:site_id, :date]),
              msx, on=[:site_id, :date])
j = dropmissing(j)
n_total = nrow(j)
j = j[j.mortality_fraction .< 1, :]            # Z undefined at 100 % mortality
j.Z = -log.(1 .- j.mortality_fraction)

design(df) = hcat(ones(nrow(df)), Float64.(df.I), Float64.(df.msx))
b = design(j) \ j.Z
r = j.Z - design(j) * b
s2 = sum(abs2, r) / (nrow(j) - 3)
se = sqrt.(diag(s2 * inv(design(j)' * design(j))))
r2 = 1 - sum(abs2, r) / sum(abs2, j.Z .- mean(j.Z))

rng = Xoshiro(20261013)
sites = unique(j.site_id)
boot = zeros(N_BOOT, 3)
for k in 1:N_BOOT
    pick = rand(rng, sites, length(sites))
    bs = vcat([j[j.site_id .== s, :] for s in pick]...)
    boot[k, :] = design(bs) \ bs.Z
end
ci = [quantile(boot[:, c], q) for c in 1:3, q in (0.025, 0.975)]

mkpath(OUT)
fmt(x) = string(round(x; digits=4))
open(joinpath(OUT, "fit.md"), "w") do io
    println(io, "# Dermo intensity vs observed oyster mortality, Maryland disease bars (scripts/p4_dermo_mortality_fit.jl)\n")
    println(io, "n = $(nrow(j)) bar-years ($(n_total - nrow(j)) with 100 % mortality excluded), ",
            "$(minimum(year.(j.date)))-$(maximum(year.(j.date))), $(length(sites)) bars. Z = -ln(1 - m), per year.\n")
    println(io, "| term | estimate (yr^-1) | OLS se | 95 % bar-bootstrap interval |")
    println(io, "|---|---|---|---|")
    for (c, name) in enumerate(("intercept", "Dermo intensity (per unit, 0-7 scale)", "MSX prevalence (per unit, 0-1)"))
        println(io, "| $name | $(fmt(b[c])) | $(fmt(se[c])) | $(fmt(ci[c, 1])) to $(fmt(ci[c, 2])) |")
    end
    println(io, "\nR^2 = $(fmt(r2)).\n")
    println(io, "| mean Dermo intensity | n | mean observed annual mortality |")
    println(io, "|---|---|---|")
    for lo in 0:0.5:4.5
        s = j[(j.I .>= lo) .& (j.I .< lo + 0.5), :]
        nrow(s) > 0 && println(io, "| $lo-$(lo + 0.5) | $(nrow(s)) | $(fmt(mean(s.mortality_fraction))) |")
    end
end

open(joinpath(ROOT, "params", "dermo_mortality_mddnr.toml"), "w") do io
    println(io, "title = \"Excess oyster mortality per unit of Dermo intensity (fitted in this project)\"\n")
    println(io, "[dermo_mortality_per_intensity_yr]")
    println(io, "value = $(round(b[2]; sigdigits=4))")
    println(io, "units = \"yr^-1 per unit mean intensity (0-7 scale)\"")
    println(io, "source = \"Fitted in this project: scripts/p4_dermo_mortality_fit.jl, OLS of annual instantaneous observed mortality on mean Dermo intensity and MSX prevalence, Maryland DNR 2024 Fall Survey report Tables 3-5 (43 bars, $(minimum(year.(j.date)))-$(maximum(year.(j.date))), n = $(nrow(j))); results in outputs/p4_dermo_mortality/fit.md\"")
    println(io, "notes = \"95 % bar-bootstrap interval $(fmt(ci[2, 1])) to $(fmt(ci[2, 2])); OLS se $(fmt(se[2])); R^2 $(fmt(r2)) with MSX in the model. Intercept $(fmt(b[1])) yr^-1 (background) and MSX term $(fmt(b[3])) yr^-1 are not used by the farm model. Observational: confounded by other causes that co-vary with Dermo.\"")
end
println(read(joinpath(OUT, "fit.md"), String))
