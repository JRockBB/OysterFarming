# Tests for src/economics.jl (PLAN.md P5.1): a hand-calculated reference example.

using DataFrames, Dates

@testset "economics" begin
    econ = EconParams(load_params(joinpath(@__DIR__, "..", "params", "economics_reference_farm.toml")))
    @test econ.price_per_oyster ≈ 18_217_846 / 23_401_019 rtol = 1e-3
    @test econ.variable_cost_per_oyster ≈ 31.19 / 100

    # Hand example (SYNTHETIC schedule): 100,000 seed in 2020; 30,000 harvested in 2022, 20,000 in 2023.
    cohorts = [Cohort("c", Date(2020, 5, 1), 100_000, 15.0)]
    h = DataFrame(date=[Date(2022, 8, 1), Date(2023, 7, 1)], count=[30_000.0, 20_000.0])
    cf = farm_cashflows(cohorts, h, econ, 2020, 2023)
    margin = econ.price_per_oyster - econ.variable_cost_per_oyster        # 0.7785 - 0.3119 = 0.4666
    @test cf.cashflow[1] ≈ -(100_000 * 0.009) - 4800 - 155_375           # = -161,075
    @test cf.cashflow[2] ≈ -4800
    @test cf.cashflow[3] ≈ 30_000 * margin - 4800                       # = 9,198
    @test cf.cashflow[4] ≈ 20_000 * margin - 4800                       # = 4,532
    by_hand = -161_075 - 4800 / 1.1 + 9_198 / 1.1^2 + 4_532 / 1.1^3
    @test npv(cf, 0.10) ≈ by_hand rtol = 1e-9
    @test npv(cf, 0.0) ≈ sum(cf.cashflow)

    extra = farm_cashflows(cohorts, h, econ, 2020, 2023; extra_cost_per_year=Dict(2022 => 1000.0))
    @test npv(extra, 0.10) ≈ npv(cf, 0.10) - 1000 / 1.1^2
    @test_throws ArgumentError farm_cashflows(cohorts, h, econ, 2021, 2023)   # seeded before horizon
    @test_throws ArgumentError farm_cashflows(cohorts, h, econ, 2020, 2022)   # harvest after horizon

    @testset "Maine MAA file: size tiers and capital per acre" begin
        me = EconParams(load_params(joinpath(@__DIR__, "..", "params", "economics_maine_maa.toml")))
        @test me.price_tiers == [(63.5, 0.6413), (76.2, 0.74)]
        @test price_at(me, 60.0) == 0.0
        @test price_at(me, 63.5) == 0.6413
        @test price_at(me, 80.0) == 0.74
        @test price_at(econ, 60.0) == econ.price_per_oyster          # no tiers: single price
        # SYNTHETIC schedule: two cohorts in one year, one the next
        cs = [Cohort("a", Date(2020, 5, 1), 100_000, 20.0), Cohort("b", Date(2020, 9, 1), 50_000, 20.0),
              Cohort("c", Date(2021, 5, 1), 120_000, 20.0)]
        @test capital_outlay(me, cs) ≈ 40_395 * 150_000 / 214_985
        @test capital_outlay(econ, cs) == econ.capital_investment
        h = DataFrame(date=[Date(2022, 8, 1), Date(2022, 9, 1)], count=[1000.0, 2000.0],
                      shell_height_mm=[70.0, 80.0])
        cf = farm_cashflows(cs, h, me, 2020, 2022)
        @test cf.cashflow[1] ≈ -150_000 * 0.009 - capital_outlay(me, cs)
        # 3-year horizon, 10-year life: 70 % of the outlay comes back in the last year
        @test cf.cashflow[3] - (1000 * (0.6413 - 0.5311) + 2000 * (0.74 - 0.5311)) ≈ 0.7 * capital_outlay(me, cs)
        @test_throws ArgumentError farm_cashflows(cs, select(h, Not(:shell_height_mm)), me, 2020, 2022)
    end
end
