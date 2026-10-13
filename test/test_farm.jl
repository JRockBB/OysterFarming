# Tests for src/farm.jl (PLAN.md P3.1–P3.3), using SYNTHETIC parameters and forcing.

using DataFrames, Dates

const FARM_DEB = DEBParams(load_params(joinpath(@__DIR__, "fixtures", "synthetic_deb_params.toml")))
const FARM_FP = FarmParams(load_params(joinpath(@__DIR__, "fixtures", "synthetic_farm_params.toml")))

"SYNTHETIC seasonal forcing for `years` years."
function synthetic_farm_forcing(years::Integer; chl=6.0)
    n = 365 * years
    d = Date(2020, 1, 1) .+ Day.(0:n-1)
    t = 0:n-1
    return DataFrame(date=d, temp_c=12 .- 7 .* cos.(2pi .* (t .- 30) ./ 365),
                     salinity_psu=fill(30.0, n), chl_ug_l=fill(chl, n))
end

@testset "farm" begin
    @testset "parameters and rules" begin
        @test FARM_FP.flushing_time_d == 2.0
        # ingestion consistent with the Holling response: {C_R} * X_K * c_T * c_S * V^(2/3)
        @test max_ingestion_ug_d(8.0, 1.0, 1.0, FARM_DEB, FARM_FP) ≈ 20.0 * FARM_DEB.X_K * 4.0
        dam = merge(load_params(joinpath(@__DIR__, "..", "params", "farm_damariscotta.toml")),
                    load_params(joinpath(@__DIR__, "..", "params", "deb_c_virginica_clearance.toml")))
        fp = FarmParams(dam)
        @test fp.flushing_time_d == 35.0 && fp.box_volume_L ≈ 123.4e9 && fp.clearance_max_L_cm2_d == 56.0
        @test fp.mortality_rate_d ≈ 0.22 / 365 rtol = 1e-2
        @test HarvestRule().market_shell_height_mm == 76.0
        @test HarvestRule().weekly_capacity == Inf
        @test_throws ArgumentError HarvestRule(market_shell_height_mm=-1)
    end

    @testset "P3.2 food depletion box model" begin
        @test depleted_food(5.0, 0.0, 5.0, FARM_FP) == 5.0          # (a) no oysters: X = X_in
        @test depleted_food(0.0, 1e6, 5.0, FARM_FP) == 0.0
        X = depleted_food(5.0, 2e7, 5.0, FARM_FP)
        @test 0 < X < 5
        resid = (5.0 - X) / FARM_FP.flushing_time_d - 2e7 * X / (X + 5.0) / FARM_FP.box_volume_L
        @test abs(resid) < 1e-8                                      # steady state of dX/dt
        @test depleted_food(5.0, 4e7, 5.0, FARM_FP) < X              # more demand, less food
    end

    @testset "P3.1 three cohorts for five years conserve counts" begin
        forcing = synthetic_farm_forcing(5)
        cohorts = [Cohort("c1", Date(2020, 5, 1), 1e5, 15.0),
                   Cohort("c2", Date(2021, 5, 1), 1e5, 15.0),
                   Cohort("c3", Date(2022, 6, 1), 5e4, 20.0)]
        daily, harvests = simulate_farm(FARM_DEB, FARM_FP, forcing, cohorts, HarvestRule())
        for c in cohorts
            rows = daily[daily.cohort_id .== c.id, :]
            @test first(rows.date) == c.seed_date
            @test all(isapprox.(rows.n_alive .+ rows.n_dead .+ rows.n_harvested, c.n_seeded; rtol=1e-12))
            @test all(diff(rows.n_dead) .>= 0)
        end
        @test nrow(harvests) > 0
    end

    @testset "background mortality follows exp(-m t)" begin
        forcing = synthetic_farm_forcing(1)
        cohorts = [Cohort("c", Date(2020, 1, 1), 1000.0, 15.0)]
        daily, _ = simulate_farm(FARM_DEB, FARM_FP, forcing, cohorts,
                                 HarvestRule(market_shell_height_mm=1e6))
        @test daily.n_alive[end] ≈ 1000 * exp(-FARM_FP.mortality_rate_d * nrow(forcing))
    end

    @testset "P3.2 (b) higher density lowers per-oyster growth" begin
        forcing = synthetic_farm_forcing(2)
        final_sh = map([1e2, 1e4, 1e5, 1e6, 1e7]) do n
            daily, _ = simulate_farm(FARM_DEB, FARM_FP, forcing, [Cohort("c", Date(2020, 3, 1), n, 15.0)],
                                     HarvestRule(market_shell_height_mm=1e6))
            daily.shell_height_mm[end]
        end
        @test issorted(final_sh; rev=true)
        @test final_sh[end] < final_sh[1]
    end

    @testset "no depletion matches the single-oyster simulation" begin
        forcing = synthetic_farm_forcing(1)
        seed = Date(2020, 3, 1)
        daily, _ = simulate_farm(FARM_DEB, FARM_FP, forcing, [Cohort("c", seed, 1e6, 15.0)],
                                 HarvestRule(market_shell_height_mm=1e6); depletion=false)
        single = simulate_individual(FARM_DEB, forcing[forcing.date .>= seed, :],
                                     initial_state(FARM_DEB; shell_height_mm=15.0, e=1.0))
        @test daily.shell_height_mm ≈ single.shell_height_mm
    end

    @testset "P3.3 harvest at market size within weekly capacity" begin
        forcing = synthetic_farm_forcing(3)
        cohorts = [Cohort("a", Date(2020, 4, 1), 2e4, 15.0), Cohort("b", Date(2020, 4, 1), 2e4, 15.0)]
        rule = HarvestRule(market_shell_height_mm=60.0, weekly_capacity=3000.0)
        daily, harvests = simulate_farm(FARM_DEB, FARM_FP, forcing, cohorts, rule)
        @test nrow(harvests) > 0
        @test all(harvests.shell_height_mm .>= rule.market_shell_height_mm)
        weekly = combine(groupby(transform(harvests, :date => ByRow(firstdayofweek) => :week), :week),
                         :count => sum => :total)
        @test all(weekly.total .<= rule.weekly_capacity + 1e-9)
        # oldest-first priority with a tie on seed date: cohort "a" (id order) is emptied first
        first_b = findfirst(==("b"), harvests.cohort_id)
        @test first_b === nothing || all(harvests.cohort_id[1:first_b-1] .== "a")
    end

    @testset "input checks" begin
        forcing = synthetic_farm_forcing(1)
        @test_throws ArgumentError simulate_farm(FARM_DEB, FARM_FP, forcing,
                                                 [Cohort("c", Date(2030, 1, 1), 1.0, 15.0)], HarvestRule())
        @test_throws ArgumentError simulate_farm(FARM_DEB, FARM_FP, forcing,
                                                 [Cohort("c", Date(2020, 2, 1), 1.0, 15.0),
                                                  Cohort("c", Date(2020, 3, 1), 1.0, 15.0)], HarvestRule())
    end
end
