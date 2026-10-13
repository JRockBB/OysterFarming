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
end
