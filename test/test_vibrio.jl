# Tests for src/vibrio.jl (PLAN.md P4.5), with the Maine rule file and SYNTHETIC farm inputs.

using DataFrames, Dates

@testset "vibrio" begin
    vr = VibrioRule(load_params(joinpath(@__DIR__, "..", "params", "vibrio_maine.toml")))

    @testset "control period boundaries (June 1 - October 15 inclusive)" begin
        @test !in_vibrio_control(Date(2024, 5, 31), vr)
        @test in_vibrio_control(Date(2024, 6, 1), vr)
        @test in_vibrio_control(Date(2024, 10, 15), vr)
        @test !in_vibrio_control(Date(2024, 10, 16), vr)
        @test in_vibrio_control(Date(2017, 8, 1), vr)          # applied to every simulated year
    end

    @testset "capacity and cost apply only inside the period" begin
        cap = vibrio_capacity(vr, 1000, 0.5)
        @test cap(Date(2024, 5, 31)) == 1000 && cap(Date(2024, 6, 1)) == 500
        @test cap(Date(2024, 10, 15)) == 500 && cap(Date(2024, 10, 16)) == 1000
        @test_throws DomainError vibrio_capacity(vr, 1000, 1.5)
        h = DataFrame(date=[Date(2024, 5, 31), Date(2024, 6, 1), Date(2024, 10, 15), Date(2024, 10, 16)],
                      count=[10.0, 20.0, 30.0, 40.0])
        @test vibrio_handling_cost(h, vr; cost_per_oyster=0.1) ≈ 0.1 * (20 + 30)
        @test_throws DomainError vibrio_handling_cost(h, vr; cost_per_oyster=-1)
    end

    @testset "reduced capacity in the farm simulation" begin
        deb = DEBParams(load_params(joinpath(@__DIR__, "fixtures", "synthetic_deb_params.toml")))
        fp = FarmParams(load_params(joinpath(@__DIR__, "fixtures", "synthetic_farm_params.toml")))
        n = 3 * 365
        d = Date(2020, 1, 1) .+ Day.(0:n-1)
        t = 0:n-1
        forcing = DataFrame(date=d, temp_c=12 .- 7 .* cos.(2pi .* (t .- 30) ./ 365),
                            salinity_psu=fill(30.0, n), chl_ug_l=fill(6.0, n))   # SYNTHETIC
        rule = HarvestRule(market_shell_height_mm=60.0, weekly_capacity=400.0)
        cohorts = [Cohort("c", Date(2020, 3, 1), 2e4, 15.0)]
        _, h = simulate_farm(deb, fp, forcing, cohorts, rule; weekly_capacity_fn=vibrio_capacity(vr, 400.0, 0.25))
        weekly = combine(groupby(transform(h, :date => ByRow(firstdayofweek) => :week), :week), :count => sum => :total)
        inside = [in_vibrio_control(w, vr) && in_vibrio_control(w + Day(6), vr) for w in weekly.week]
        @test any(inside)
        @test all(weekly.total[inside] .<= 100.0 + 1e-9)       # 400 x 0.25 inside the period
        @test any(weekly.total[.!inside] .> 100.0)
    end
end
