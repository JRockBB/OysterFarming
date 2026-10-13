# Tests for src/decision.jl (PLAN.md P5.2-P5.3) with SYNTHETIC biology and forcing and the
# project economics, nitrogen and Vibrio parameter files.

using DataFrames, Dates

@testset "decision" begin
    syn = DEBParams(load_params(joinpath(@__DIR__, "fixtures", "synthetic_deb_params.toml")))
    fp = FarmParams(load_params(joinpath(@__DIR__, "fixtures", "synthetic_farm_params.toml")))
    econ = EconParams(load_params(joinpath(@__DIR__, "..", "params", "economics_reference_farm.toml")))
    nit = NitrogenParams(load_params(joinpath(@__DIR__, "..", "params", "nitrogen_chesapeake_bmp.toml")))
    vr = VibrioRule(load_params(joinpath(@__DIR__, "..", "params", "vibrio_maine.toml")))
    n = 4 * 365
    d = Date(2020, 1, 1) .+ Day.(0:n-1)
    t = 0:n-1
    forcing = DataFrame(date=d, temp_c=12 .- 7 .* cos.(2pi .* (t .- 30) ./ 365),
                        salinity_psu=fill(30.0, n), chl_ug_l=fill(6.0, n))           # SYNTHETIC
    scen = [(0.5, with_values(syn; X_K=4.0)), (0.5, with_values(syn; X_K=8.0))]
    setup(; vf=1.0, vc=0.0) = DecisionSetup(forcing, scen, fp, econ, nit, vr, 5000.0, vf, vc, 15.0, 2, 0.1)

    @test weighted_quantile([3.0, 1.0, 2.0], [0.2, 0.5, 0.3], 0.5) == 1.0
    @test weighted_quantile([3.0, 1.0, 2.0], [0.2, 0.5, 0.3], 0.6) == 2.0

    r = evaluate_decision(FarmDecision(5, 50_000, 60.0), setup())
    @test length(r.npvs) == 2
    @test r.expected_npv ≈ 0.5 * r.npvs[1] + 0.5 * r.npvs[2]
    @test r.lower_npv == minimum(r.npvs)                 # 10 % quantile of two equal-weight cases
    @test r.npvs[1] > r.npvs[2]                          # more food-limited case earns less
    @test r.expected_harvested > 0 && r.expected_nitrogen_kg > 0

    # a Vibrio handling cost lowers NPV; it cannot raise it
    rv = evaluate_decision(FarmDecision(5, 50_000, 60.0), setup(vc=0.05))
    @test rv.expected_npv < r.expected_npv

    g = grid_search(setup(); seed_months=[4, 6], seeds_per_cohort=[20_000, 50_000], market_sizes_mm=[60.0, 76.0])
    @test nrow(g) == 8
    @test issorted(g.expected_npv; rev=true)
    @test_throws ArgumentError evaluate_decision(FarmDecision(13, 1.0, 60.0), setup())
end
