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

@testset "time to market" begin
    p = DEBParams(load_params(joinpath(@__DIR__, "fixtures", "synthetic_deb_params.toml")))
    n = 800
    d = Date(2020, 1, 1) .+ Day.(0:n-1)
    forcing = DataFrame(date=d, temp_c=fill(18.0, n), salinity_psu=fill(30.0, n), chl_ug_l=fill(10.0, n))  # SYNTHETIC
    u0 = initial_state(p; shell_height_mm=20.0, e=0.5)
    days = time_to_market(p, forcing, d[1], u0)
    @test days isa Int
    sh = simulate_individual(p, forcing, u0).shell_height_mm
    @test sh[days] >= 76.0 && sh[days-1] < 76.0
    @test time_to_market(p, forcing, d[1], u0; market_shell_height_mm=1e6) === nothing
    @test_throws ArgumentError time_to_market(p, forcing, d[1] - Day(1), u0)
    gappy = forcing[[1:30; 200:n], :]
    @test time_to_market(p, gappy, d[1], u0) === nothing   # stops at the gap
end

@testset "cycle_years" begin
    # SYNTHETIC daily series 2016-2017; value = day index so re-dating can be checked
    d = collect(Date(2016, 1, 1):Day(1):Date(2017, 12, 31))
    f = DataFrame(site_id="x", date=d, temp_c=Float64.(1:length(d)))
    c = cycle_years(f, [2016, 2017], 5)
    @test first(c.date) == Date(2016, 1, 1) && last(c.date) == Date(2020, 12, 31)
    @test all(diff(c.date) .== Day(1))
    @test nrow(c) == sum(daysinyear(y) for y in 2016:2020)
    @test c.temp_c[c.date .== Date(2018, 3, 1)] == f.temp_c[f.date .== Date(2016, 3, 1)]   # 2016 placed in 2018, Feb 29 dropped
    @test c.temp_c[c.date .== Date(2020, 2, 29)] == f.temp_c[f.date .== Date(2016, 2, 29)]  # 2016 again in 2020 (leap)
    g = cycle_years(f, [2017], 4; start_year=2019)       # 2017 placed in leap 2020: Feb 28 repeated
    @test g.temp_c[g.date .== Date(2020, 2, 29)] == f.temp_c[f.date .== Date(2017, 2, 28)]
    @test_throws ArgumentError cycle_years(f[1:100, :], [2016], 2)
end
