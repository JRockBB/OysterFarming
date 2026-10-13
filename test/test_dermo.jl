# Tests for src/dermo.jl (PLAN.md P4.1) with the published Ragone Calvo et al. (2001) parameters
# and SYNTHETIC constant or seasonal forcing.

using DataFrames, Dates

const DP = DermoParams(load_params(joinpath(@__DIR__, "..", "params", "dermo_ragone_calvo2001.toml")))

"SYNTHETIC forcing: `n` days from `start` at constant temperature and salinity."
dermo_forcing(n, temp_c, sal; start=Date(2001, 1, 1)) =
    DataFrame(date=start .+ Day.(0:n-1), temp_c=fill(temp_c, n), salinity_psu=fill(sal, n))

@testset "dermo" begin
    @testset "rate functions match the published equations" begin
        @test dermo_division_rate(20.0, 20.0, DP) == 0.0                 # no division at or below 20 degC
        @test dermo_division_rate(25.0, 20.0, DP) ≈ (-0.009 + 0.0035 * 20) * exp(0.06931 * 5)
        @test dermo_mortality_rate(18.0, 20.0, DP) == 0.0                # only for 3 < T < 18
        @test dermo_mortality_rate(3.0, 20.0, DP) == 0.0
        @test dermo_mortality_rate(10.0, 10.0, DP) ≈ 0.02 * exp(0.06931 * 8)
        @test dermo_mortality_rate(10.0, 13.0, DP) ≈ 0.0175 * exp(0.06931 * 8)
        @test dermo_mortality_rate(10.0, 2.0, DP) ≈ 0.006 * exp(0.06931 * 8)
        @test [dermo_transmission_cells(s, DP) for s in (20, 13, 10, 7, 5)] == [1000, 500, 100, 50, 0]
        @test rank_to_log10_cells(0) ≈ 0.176
        @test rank_to_log10_cells(6) ≈ 0.176 - 0.463 * 6 + 0.205 * 36
    end

    @testset "transmission happens once a year on day 218" begin
        out = simulate_dermo(DP, dermo_forcing(365, 19.0, 20.0), 0.0)   # 19 degC: no division, no mortality
        jump = findall(>(0), diff(vcat(0.0, out.P_cells_g)))
        @test jump == [218]
        @test out.P_cells_g[end] ≈ 1000
    end

    @testset "intensity rises with temperature and salinity (P4.1 done-when)" begin
        # days for an infection of 1000 cells/g to reach 5e4 cells/g (below the 1e5 host-death loss
        # threshold, so the comparison is not confounded by the loss term)
        function days_to(T, S)
            out = simulate_dermo(DP, dermo_forcing(365, T, S; start=Date(2001, 5, 1)), 1000.0)
            i = findfirst(>=(5e4), out.P_cells_g)
            return i === nothing ? typemax(Int) : i
        end
        @test days_to(28.0, 20.0) < days_to(24.0, 20.0) < days_to(21.0, 20.0)
        @test days_to(25.0, 25.0) < days_to(25.0, 15.0) < days_to(25.0, 8.0)
        # winter: parasites decline, faster at 8-12 psu than at lower salinity
        wint(S) = simulate_dermo(DP, dermo_forcing(90, 8.0, S; start=Date(2001, 1, 1)), 1e4).P_cells_g[end]
        @test wint(10.0) < wint(2.0) < 1e4
    end

    @testset "density feedback and host-death loss keep P bounded" begin
        out = simulate_dermo(DP, dermo_forcing(400, 30.0, 30.0; start=Date(2001, 1, 1)), 1e4)
        @test all(>=(0), out.P_cells_g)
        @test maximum(out.P_cells_g) < DP.G_p
        @test any(out.P_cells_g .> DP.loss_P_threshold)
    end

    @testset "seasonal cycle: summer build-up, winter decline" begin
        n = 730
        d = Date(2001, 1, 1) .+ Day.(0:n-1)
        t = dayofyear.(d)
        forcing = DataFrame(date=d, temp_c=15 .- 11 .* cos.(2pi .* (t .- 20) ./ 365), salinity_psu=fill(16.0, n))
        out = simulate_dermo(DP, forcing, 500.0)
        aug = out.P_cells_g[d .== Date(2002, 8, 31)][1]
        mar = out.P_cells_g[d .== Date(2002, 3, 31)][1]
        @test aug > mar
    end

    @testset "input checks" begin
        @test_throws ArgumentError simulate_dermo(DP, dermo_forcing(0, 20.0, 20.0), 1.0)
        @test_throws DomainError simulate_dermo(DP, dermo_forcing(5, 20.0, 20.0), -1.0)
        f = dermo_forcing(5, 20.0, 20.0)
        @test_throws ArgumentError simulate_dermo(DP, f[[1, 3], :], 1.0)
    end
end

@testset "dermo intensity and oyster mortality" begin
    for x in (1.5, 3.0, 5.0, 7.0)
        @test cells_to_intensity(10^rank_to_log10_cells(x)) ≈ x
    end
    @test cells_to_intensity(0.0) == 0.0
    @test cells_to_intensity(0.5) == 0.0
    @test cells_to_intensity(1e5) < cells_to_intensity(1e6)
    @test_throws DomainError cells_to_intensity(-1.0)
    @test dermo_excess_mortality_d(10^rank_to_log10_cells(3.0), 0.09) ≈ 0.09 * 3 / 365
end
