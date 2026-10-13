# Tests for src/stochastic.jl and the closure / extra-mortality hooks in simulate_farm
# (PLAN.md P4.6). SYNTHETIC parameters and forcing.

using DataFrames, Dates, Random, Statistics

@testset "stochastic" begin
    deb = DEBParams(load_params(joinpath(@__DIR__, "fixtures", "synthetic_deb_params.toml")))
    fp = FarmParams(load_params(joinpath(@__DIR__, "fixtures", "synthetic_farm_params.toml")))
    n = 2 * 365
    d = Date(2020, 1, 1) .+ Day.(0:n-1)
    t = 0:n-1
    forcing = DataFrame(date=d, temp_c=12 .- 7 .* cos.(2pi .* (t .- 30) ./ 365),
                        salinity_psu=fill(30.0, n), chl_ug_l=fill(6.0, n))     # SYNTHETIC
    cohorts = [Cohort("c", Date(2020, 3, 1), 1e4, 15.0)]
    rule = HarvestRule(market_shell_height_mm=60.0, weekly_capacity=500.0)

    @testset "closures block harvest and are counted" begin
        open_daily, open_h = simulate_farm(deb, fp, forcing, cohorts, rule)
        first_harvest = minimum(open_h.date)
        closed = Set(first_harvest .+ Day.(0:13))                  # two-week closure at harvest start
        daily, h = simulate_farm(deb, fp, forcing, cohorts, rule; closed_dates=closed)
        @test isempty(intersect(Set(h.date), closed))
        out = farm_outcome(daily, h)
        @test out.lost_harvest_days == 14
        @test farm_outcome(open_daily, open_h).lost_harvest_days == 0
        @test out.seeded ≈ 1e4
    end

    @testset "extra mortality reduces survivors" begin
        base, _ = simulate_farm(deb, fp, forcing, cohorts, HarvestRule(market_shell_height_mm=1e6))
        sick, _ = simulate_farm(deb, fp, forcing, cohorts, HarvestRule(market_shell_height_mm=1e6);
                                extra_mortality_d=date -> month(date) in (7, 8) ? 0.01 : 0.0)
        @test sick.n_alive[end] < base.n_alive[end]
        @test sick.n_alive[end] + sick.n_dead[end] ≈ 1e4
        @test_throws DomainError simulate_farm(deb, fp, forcing, cohorts, rule; extra_mortality_d=d -> -1.0)
    end

    @testset "ensemble is reproducible and varies across members" begin
        function member(rng)
            xk = deb.X_K * exp(0.3 * randn(rng))                    # SYNTHETIC parameter uncertainty
            closure_start = Date(2021, 6, 1) + Day(rand(rng, 0:60))
            closed = Set(closure_start .+ Day.(0:9))
            daily, h = simulate_farm(with_values(deb; X_K=xk), fp, forcing, cohorts, rule; closed_dates=closed)
            return merge((X_K=xk,), farm_outcome(daily, h))
        end
        a = run_ensemble(Xoshiro(1), 6, member)
        b = run_ensemble(Xoshiro(1), 6, member)
        c = run_ensemble(Xoshiro(2), 6, member)
        @test a == b                                                 # same seed, identical results
        @test a != c
        @test nrow(a) == 6 && a.member == 1:6
        @test length(unique(a.X_K)) == 6
        @test all(a.seeded .≈ 1e4)
        @test_throws ArgumentError run_ensemble(Xoshiro(1), 0, member)
    end
end
