# Tests for src/contract.jl using the SYNTHETIC fixtures in data/fixtures/
# (written by data-pipeline/scripts/make_fixtures.py). Reading the Python-written
# files here also checks that the two languages agree on the Arrow encoding.

using DataFrames, Dates

const FIXTURES = joinpath(@__DIR__, "..", "data", "fixtures")
fixture(name, case) = joinpath(FIXTURES, "synthetic_$(name)_$(case).arrow")

@testset "contract" begin
    @testset "valid fixture accepted: $name" for name in sort(collect(keys(CONTRACT_SPECS)))
        df = read_contract(fixture(name, "valid"), name)
        @test nrow(df) > 0
        validator = getfield(OysterFarming, Symbol("validate_", name))
        @test validator(df) === nothing
    end

    @testset "broken fixtures rejected: $name" for name in sort(collect(keys(CONTRACT_SPECS)))
        for (case, needle) in [("missing_column", "missing column"),
                               ("wrong_type", "has type"),
                               ("out_of_range", "value")]
            err = try
                read_contract(fixture(name, case), name)
                nothing
            catch e
                e
            end
            @test err isa ContractError
            @test err !== nothing && any(p -> occursin(needle, p), err.problems)
        end
    end

    forcing() = read_contract(fixture("forcing_daily", "valid"), "forcing_daily")

    @testset "Arrow types round-trip from Python" begin
        df = forcing()
        @test eltype(df.date) == Date
        @test eltype(df.qc_flag) == Int8
        @test nonmissingtype(eltype(df.chl_ug_l)) == Float64
        @test count(ismissing, df.chl_ug_l) == 1
    end

    @testset "temperature bounds are strict" begin
        df = forcing()
        df.temp_c[1] = -2.0
        @test_throws ContractError validate_forcing_daily(df)
    end

    @testset "NaN rejected" begin
        df = forcing()
        df.salinity_psu[3] = NaN
        @test_throws ContractError validate_forcing_daily(df)
    end

    @testset "missing in required column rejected" begin
        df = forcing()
        df.salinity_psu = allowmissing(df.salinity_psu)
        df.salinity_psu[3] = missing
        @test_throws ContractError validate_forcing_daily(df)
    end

    @testset "duplicate site-day rejected" begin
        df = forcing()
        df.date[2] = df.date[1]
        @test_throws ContractError validate_forcing_daily(df)
    end

    @testset "unexpected column rejected" begin
        df = forcing()
        df.extra = ones(nrow(df))
        @test_throws ContractError validate_forcing_daily(df)
    end

    @testset "closure end before start rejected" begin
        df = read_contract(fixture("closures", "valid"), "closures")
        df.end_date[1] = Date(2020, 3, 1)
        @test_throws ContractError validate_closures(df)
    end

    @testset "unknown table name" begin
        @test_throws ArgumentError read_contract(fixture("sites", "valid"), "nope")
    end

    @testset "error message lists problems" begin
        msg = sprint(showerror, ContractError("t", ["a", "b"]))
        @test occursin("t violates", msg) && occursin("- b", msg)
    end
end
