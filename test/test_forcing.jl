# Tests for src/forcing.jl using the SYNTHETIC forcing fixture.

using DataFrames, Dates

@testset "forcing helpers" begin
    df = read_contract(joinpath(@__DIR__, "..", "data", "fixtures", "synthetic_forcing_daily_valid.arrow"),
                       "forcing_daily")
    site = forcing_for_site(df, "synthetic_site")
    @test nrow(site) == nrow(df)
    @test issorted(site.date)
    @test_throws ArgumentError forcing_for_site(df, "nope")

    d = Date(2020, 1, 1) .+ Day.([0, 1, 2, 5, 6, 10])
    @test continuous_blocks(d) == [1:3, 4:5, 6:6]
    @test continuous_blocks(Date[]) == UnitRange{Int}[]
    @test_throws ArgumentError continuous_blocks(reverse(d))
    gappy = site[[1, 2, 3, 10, 11, 12, 13], :]
    @test nrow(longest_block(gappy)) == 4
end
