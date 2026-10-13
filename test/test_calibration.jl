# Tests for src/calibration.jl (PLAN.md P2.5). All data here are SYNTHETIC:
# shell heights simulated from the synthetic parameter set plus seeded noise.

using DataFrames, Dates, Random

"SYNTHETIC calibration setup: seasonal forcing, observations every 10 days."
function synthetic_calibration_problem(; n_days=240, noise_sd_mm=0.5, seed=42)
    p = DEBParams(load_params(joinpath(@__DIR__, "fixtures", "synthetic_deb_params.toml")))
    d = Date(2020, 4, 1) .+ Day.(0:n_days-1)
    t = 0:n_days-1
    forcing = DataFrame(date=d, temp_c=14 .+ 6 .* sin.(2pi .* t ./ 365), salinity_psu=fill(30.0, n_days),
                        chl_ug_l=4 .+ 3 .* sin.(2pi .* t ./ 60))
    u0 = initial_state(p; shell_height_mm=20.0, e=0.5)
    obs_i = 10:10:n_days
    clean = simulate_individual(p, forcing, u0).shell_height_mm[obs_i]
    y = clean .+ noise_sd_mm .* randn(Xoshiro(seed), length(obs_i))
    cp = CalibrationProblem(p, forcing, u0, d[obs_i], y; names=[:X_K, :p_M],
                            lower=[1.0, 10.0], upper=[20.0, 50.0])
    return cp, [p.X_K, p.p_M]
end

@testset "calibration" begin
    cp, truth = synthetic_calibration_problem()

    @testset "with_values and generic types" begin
        q = with_values(cp.base; X_K=7.0)
        @test q.X_K == 7.0 && q.p_M == cp.base.p_M
        @test eltype(with_values(cp.base; X_K=big(7.0))) == BigFloat
        @test_throws ArgumentError with_values(cp.base; nope=1.0)
        @test predict_shell_height(cp, truth) ≈
              simulate_individual(cp.base, cp.forcing, cp.u0).shell_height_mm[cp.obs_idx]
    end

    @testset "problem validation" begin
        d = cp.forcing.date
        @test_throws ArgumentError CalibrationProblem(cp.base, cp.forcing, cp.u0, [d[1] - Day(1)], [1.0];
                                                      names=[:X_K], lower=[1.0], upper=[2.0])
        @test_throws ArgumentError CalibrationProblem(cp.base, cp.forcing, cp.u0, d[1:2], [1.0, 2.0];
                                                      names=[:X_K], lower=[3.0], upper=[2.0])
        @test_throws ArgumentError CalibrationProblem(cp.base, cp.forcing, cp.u0, d[1:2], [1.0, 2.0];
                                                      names=[:nope], lower=[1.0], upper=[2.0])
        @test_throws ArgumentError fit_least_squares(cp, [0.0, 30.0])
    end

    @testset "least squares recovers known parameters" begin
        ls = fit_least_squares(cp, [10.0, 35.0])
        @test ls.converged
        @test all(ls.lower_ci .<= truth .<= ls.upper_ci)
        @test isapprox(ls.estimate, truth; rtol=0.15)
        @test 0.3 < ls.sigma_mm < 0.8    # synthetic noise sd is 0.5 mm
    end

    # NUTS sampling takes minutes (much longer under Pkg.test's forced bounds
    # checking), so it runs only when OYSTER_SLOW_TESTS=1. CI sets this in a
    # separate job. Verified locally 2026-10-11: posterior means (5.22, 23.5)
    # and 95 % intervals contain the true values (5, 25).
    get(ENV, "OYSTER_SLOW_TESTS", "0") == "1" && @testset "Bayesian NUTS recovers known parameters" begin
        b = fit_bayesian(Xoshiro(7), cp; sigma_prior_scale_mm=2.0, n_samples=200, n_adapt=200)
        @test all(b.lower_ci .<= truth .<= b.upper_ci)
        @test isapprox(b.mean, truth; rtol=0.15)
        @test 0.3 < b.sigma_mean_mm < 0.8
    end
end

@testset "initial state as a function of parameters" begin
    cp0, truth = synthetic_calibration_problem(n_days=120)
    u0f = p -> initial_state(p; shell_height_mm=20.0, e=functional_response(4.0, p.X_K))
    cpf = CalibrationProblem(cp0.base, cp0.forcing, u0f, cp0.forcing.date[cp0.obs_idx], cp0.obs_shell_height_mm;
                             names=[:X_K], lower=[1.0], upper=[20.0])
    a = predict_shell_height(cpf, [5.0])
    b = predict_shell_height(cpf, [10.0])
    @test a != b
    @test length(a) == length(cp0.obs_idx)
end
