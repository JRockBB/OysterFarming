# Tests for src/params.jl using SYNTHETIC parameter fixtures in test/fixtures/.

pfix(name) = joinpath(@__DIR__, "fixtures", "synthetic_params_$(name).toml")

@testset "params" begin
    @testset "valid file loads" begin
        p = load_params(pfix("valid"))
        @test Set(keys(p)) == Set(["a", "b"])
        @test p["a"].value === 2.0
        @test p["b"].units == "d^-1"
        @test p["b"].notes == "has notes"
        @test p["a"].notes == ""
        @test param_value(p, "b") == 0.5
        @test_throws KeyError param_value(p, "missing")
    end

    @testset "refusals" begin
        for (case, needle) in [("missing_source", "missing source"),
                               ("todo_source", "TODO(source)"),
                               ("missing_units", "missing units")]
            err = try
                load_params(pfix(case))
                nothing
            catch e
                e
            end
            @test err isa ParamError
            @test err !== nothing && any(p -> occursin(needle, p), err.problems)
        end
    end

    @testset "unknown key and missing value" begin
        path = tempname() * ".toml"
        write(path, "[a]\nunits = \"1\"\nsource = \"synthetic\"\nvalu = 1.0\n")
        err = try
            load_params(path)
            nothing
        catch e
            e
        end
        @test err isa ParamError
        @test any(p -> occursin("unknown key", p), err.problems)
        @test any(p -> occursin("missing value", p), err.problems)
        @test occursin("invalid parameter file", sprint(showerror, err))
    end
end

@testset "incomplete project parameter files are refused" begin
    for f in ("deb_c_virginica.toml", "farm_damariscotta.toml")
        @test_throws ParamError load_params(joinpath(@__DIR__, "..", "params", f))
    end
end
