# P2.2 check: run the DEB model on real forcing_daily data for a Maine site.
#
# METHOD DEMONSTRATION ONLY. It uses the Pacific oyster (C. gigas) parameter set
# of Pouvreau et al. (2006), because the eastern oyster set is not complete yet
# (params/deb_c_virginica.toml). The output says nothing about eastern oyster
# growth in Maine.
#
# Usage (from the repository root):
#   julia --project scripts/p2_2_damariscotta_run.jl [site_id]
#
# Writes outputs/p2_2/<site_id>_gigas_demo.arrow and prints a monthly summary.

using OysterFarming
using Arrow, DataFrames, Dates

const ROOT = normpath(joinpath(@__DIR__, ".."))
site_id = isempty(ARGS) ? "dam_upper" : ARGS[1]

# Scenario inputs, not model parameters: a 20 mm seed oyster with half-full
# reserve and an empty reproduction buffer at the start of the forcing block.
const SEED_SHELL_HEIGHT_MM = 20.0
const SEED_SCALED_RESERVE = 0.5

forcing = read_contract(joinpath(ROOT, "data", "processed", "forcing_daily.arrow"), "forcing_daily")
site = forcing_for_site(forcing, site_id)
block = longest_block(site)
println("site $site_id: $(nrow(site)) forcing days; longest unbroken block ",
        first(block.date), " to ", last(block.date), " (", nrow(block), " days, ",
        round(100 * count(block.qc_flag .> 0) / nrow(block), digits=1), "% gap-filled)")

p = DEBParams(load_params(joinpath(ROOT, "params", "deb_c_gigas_pouvreau2006.toml")))
u0 = initial_state(p; shell_height_mm=SEED_SHELL_HEIGHT_MM, e=SEED_SCALED_RESERVE)
out = simulate_individual(p, block, u0)

monthly = combine(groupby(transform(out, :date => ByRow(d -> Date(year(d), month(d), 1)) => :month), :month),
                  :shell_height_mm => last => :shell_height_mm,
                  :dry_tissue_g => last => :dry_tissue_g,
                  :spawned_J => (x -> count(>(0), x)) => :spawn_events)
show(stdout, MIME("text/plain"), monthly; allrows=true)
println()

outdir = joinpath(ROOT, "outputs", "p2_2")
mkpath(outdir)
dest = joinpath(outdir, "$(site_id)_gigas_demo.arrow")
Arrow.write(dest, out)
println("wrote $dest")
