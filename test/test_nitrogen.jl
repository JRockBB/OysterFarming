# Tests for src/nitrogen.jl (PLAN.md P3.4) against hand-calculated values.

using DataFrames

@testset "nitrogen" begin
    np = NitrogenParams(load_params(joinpath(@__DIR__, "..", "params", "nitrogen_chesapeake_bmp.toml")))
    # Hand calculation for a 76 mm diploid oyster:
    #   tissue DW = 0.0004 * 76^1.82 = 1.05960 g; N = 0.082 * 1.05960 = 0.086887 g
    #   shell DW  = 0.00147 * 76^2.3964 = 47.2606 g; N = 0.002 * 47.2606 = 0.094521 g
    @test nitrogen_per_oyster_g(76.0, np) ≈ 0.0868870 rtol = 1e-5
    @test nitrogen_per_oyster_g(76.0, np; include_shell=true) ≈ 0.0868870 + 0.0945212 rtol = 1e-5
    # triploid: DW = 0.00005 * 76^2.39 = 1.56356 g; N = 0.128212 g
    @test nitrogen_per_oyster_g(76.0, np; ploidy=:triploid) ≈ 0.128212 rtol = 1e-5
    @test nitrogen_per_oyster_g(0.0, np) == 0.0
    @test_throws ArgumentError nitrogen_per_oyster_g(76.0, np; ploidy=:tetraploid)
    @test_throws DomainError nitrogen_per_oyster_g(-1.0, np)

    harvests = DataFrame(count=[1000.0, 500.0], shell_height_mm=[76.0, 76.0])
    @test nitrogen_removed_g(harvests, np) ≈ 1500 * 0.0868870 rtol = 1e-5
    @test nitrogen_removed_g(harvests[1:0, :], np) == 0.0
end
