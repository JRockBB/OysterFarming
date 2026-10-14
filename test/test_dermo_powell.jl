# Tests for src/dermo_powell.jl (PLAN.md P4.1b). States below are SYNTHETIC.

using DataFrames, Dates

@testset "dermo_powell" begin
    p = PowellDermoParams(load_params(joinpath(@__DIR__, "..", "params", "dermo_powell1996.toml")))
    nj = length(p.W)
    @test nj == 11
    @test p.r_m0 ≈ log(2) * 24 / 60
    @test p.lethal_level[7] == 28 && p.lethal_level[1] == 23
    @test PowellDermoParams(load_params(joinpath(@__DIR__, "..", "params", "dermo_powell1996.toml"));
                            halving_time_h=30).r_m0 ≈ 0.5545 atol = 1e-4

    @testset "rates (Hofmann et al. 1995 eqs. 9-13)" begin
        @test powell_division_rate(20, 20, p) ≈ p.r_d0
        @test powell_division_rate(30, 20, p) ≈ 2 * p.r_d0 rtol = 1e-3          # Q10 = 2
        @test powell_division_rate(20, 5, p) ≈ p.r_d1 / 2
        @test powell_mortality_rate(20, 20, p) ≈ p.r_m0                           # base point
        @test powell_mortality_rate(5, 20, p) ≈ powell_mortality_rate(10, 20, p)  # constant below 10 degC
        @test powell_mortality_rate(25, 25, p) < powell_mortality_rate(20, 20, p)
        # density effect: no limit at light infections, both rates limited at heavy ones
        @test powell_net_rate(25, 20, 5, 7, p) ≈ powell_division_rate(25, 20, p) - powell_mortality_rate(25, 20, p)
        heavy = powell_net_rate(25, 20, 27, 7, p)
        @test abs(heavy) < abs(powell_division_rate(25, 20, p) - powell_mortality_rate(25, 20, p))
        @test powell_net_rate(25, 20, 1, 7, p) == 0.0
    end

    @testset "Mackin conversion (eq. 23)" begin
        W = p.W[7]
        k = 22                                         # 2^20 cells
        @test powell_mackin(k, 7, p) ≈ log10(2.0^20 / (p.upsilon * p.tau * W)) / p.phi
        @test powell_mackin(14, 7, p) == 0.0           # 2^12 cells: scored negative
        @test powell_mackin(1, 7, p) == 0.0
        @test issorted([powell_mackin(k, 7, p) for k in 1:POWELL_N_LEVELS])
        @test powell_mackin(POWELL_N_LEVELS, 1, p) == 5.0
    end

    @testset "conservation, transmission and deaths" begin
        N = zeros(nj, POWELL_N_LEVELS)
        N[7, 1] = 20.0
        rt = powell_transmission_rate(N, p)
        @test rt ≈ (p.r_tb + p.r_t0 * (20 / p.oyster_density_sat) / 3) / 365
        # cold water: parasites decline, nobody dies, oysters are conserved, level 2 keeps its oysters
        N = zeros(nj, POWELL_N_LEVELS); N[7, 2] = 5.0; N[7, 10] = 5.0; N[7, 1] = 10.0
        d = powell_step!(N, 5.0, 20.0, p)
        @test d == 0.0
        @test sum(N) ≈ 20.0
        @test N[7, 9] > 0 && N[7, 2] >= 5.0
        # warm water: infections intensify until oysters reach the lethal level
        f = DataFrame(date=Date(2020, 6, 1):Day(1):Date(2020, 9, 30), temp_c=fill(28.0, 122), salinity_psu=fill(20.0, 122))
        N0 = zeros(nj, POWELL_N_LEVELS); N0[7, 1] = 10.0; N0[7, 14] = 10.0
        out, Nend = simulate_powell_dermo(p, f, N0)
        @test sum(out.dermo_deaths) > 0
        @test all(isapprox.(out.density, 20.0; rtol=1e-9))          # deaths replaced by uninfected oysters
        @test out.mean_mackin[30] > powell_summary(N0, p).mean_mackin
        out2, _ = simulate_powell_dermo(p, f, N0; replace_dead=false)
        @test out2.density[end] < 20.0
        @test_throws ArgumentError simulate_powell_dermo(p, f, zeros(3, 3))
    end
end
