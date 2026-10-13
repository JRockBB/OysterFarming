# Tests for src/deb.jl. Model-property tests use the SYNTHETIC parameter set in
# test/fixtures/synthetic_deb_params.toml and synthetic constant forcing.

using DataFrames, Dates

const SYN = DEBParams(load_params(joinpath(@__DIR__, "fixtures", "synthetic_deb_params.toml")))
const PARAMS_DIR = joinpath(@__DIR__, "..", "params")

"Synthetic constant forcing for `n` days (test data, not observations)."
function synthetic_forcing(n; temp_c=20.0, salinity_psu=30.0, chl_ug_l=20.0, start=Date(2020, 1, 1))
    return DataFrame(date=start .+ Day.(0:n-1), temp_c=fill(temp_c, n),
                     salinity_psu=fill(salinity_psu, n), chl_ug_l=collect(float.(chl_ug_l isa Number ? fill(chl_ug_l, n) : chl_ug_l)))
end

@testset "deb" begin
    @testset "(a) temperature correction is 1 at the reference temperature" begin
        T_ref_c = SYN.T_ref - 273.15
        @test arrhenius_factor(T_ref_c, SYN) ≈ 1 atol = 1e-12
        @test arrhenius_factor(T_ref_c, SYN; normalize=false) < 1
        @test arrhenius_factor(25.0, SYN) > 1
        @test arrhenius_factor(15.0, SYN) < 1
        @test arrhenius_factor(40.0, SYN) < arrhenius_factor(30.0, SYN)   # above T_H rates drop
        gig = DEBParams(load_params(joinpath(PARAMS_DIR, "deb_c_gigas_pouvreau2006.toml")))
        @test arrhenius_factor(20.0, gig) ≈ 1 atol = 1e-12
    end

    @testset "functional response and salinity factor" begin
        @test functional_response(SYN.X_K, SYN.X_K) == 0.5
        @test functional_response(0.0, SYN.X_K) == 0.0
        @test_throws DomainError functional_response(-1.0, SYN.X_K)
        @test salinity_factor(30.0, SYN) == 1.0
        @test salinity_factor(10.0, SYN) == 1.0
        @test salinity_factor(6.5, SYN) ≈ 0.5
        @test salinity_factor(2.0, SYN) == 0.0
    end

    @testset "(b) constant food and temperature give von Bertalanffy growth" begin
        food = 20.0
        f = functional_response(food, SYN.X_K)
        L_m = SYN.kappa * SYN.p_Am / SYN.p_M          # maximum structural length, cm
        g = SYN.E_G / (SYN.kappa * SYN.E_m)           # energy investment ratio
        k_M = SYN.p_M / SYN.E_G                       # somatic maintenance rate coefficient, d^-1
        r_B = k_M * g / (3 * (f + g))                 # von Bertalanffy growth rate, d^-1
        L_inf = f * L_m
        V0 = 1.0
        u0 = (E=f * SYN.E_m * V0, V=V0, E_R=0.0)      # reserve at equilibrium e = f
        n = 1500
        out = simulate_individual(SYN, synthetic_forcing(n; chl_ug_l=food), u0)
        t = 1:n
        L_sim = cbrt.(out.V_cm3)
        L_vb = L_inf .- (L_inf - cbrt(V0)) .* exp.(-r_B .* t)
        @test maximum(abs.(L_sim .- L_vb) ./ L_vb) < 1e-5
        @test out.E_J ./ out.V_cm3 ≈ fill(f * SYN.E_m, n) rtol = 1e-6   # e stays at f
    end

    @testset "(c) no state goes negative, including long starvation" begin
        food = vcat(fill(20.0, 300), fill(0.0, 900), fill(20.0, 300))
        n = length(food)
        out = simulate_individual(SYN, synthetic_forcing(n; temp_c=22.0, chl_ug_l=food),
                                  initial_state(SYN; shell_height_mm=40.0))
        for c in (:E_J, :V_cm3, :E_R_J, :shell_height_mm, :dry_tissue_g, :gsi, :spawned_J)
            @test all(>=(0), out[!, c])
        end
        starve = 301:1200
        @test out.E_R_J[starve[end]] == 0                          # buffer used up first
        @test out.V_cm3[starve[end]] < out.V_cm3[starve[1]]        # then structure shrinks
    end

    # Standard DEB lets a well-fed animal keep building structure from its reserve
    # for a while after food stops (growth continues while kappa*p_C > p_M). So
    # PLAN.md P2.3(d) is tested as: reserve always shrinks; growth stops once the
    # reserve is depleted; an animal already below the growth threshold never grows.
    @testset "(d) zero food shrinks reserve, not growth" begin
        u0 = initial_state(SYN; shell_height_mm=40.0, e=1.0)
        out = simulate_individual(SYN, synthetic_forcing(400; chl_ug_l=0.0), u0)
        E = vcat(u0.E, out.E_J)
        V = vcat(u0.V, out.V_cm3)
        @test all(diff(E) .< 0)
        peak = argmax(V)
        @test peak < length(V)
        @test all(diff(V[peak:end]) .<= 1e-12)

        lean = initial_state(SYN; shell_height_mm=40.0, e=0.05)
        out = simulate_individual(SYN, synthetic_forcing(60; chl_ug_l=0.0), lean)
        V = vcat(lean.V, out.V_cm3)
        @test all(diff(V) .<= 1e-12)
        @test all(diff(vcat(lean.E, out.E_J)) .< 0)
    end

    @testset "spawning empties the buffer when warm and ripe" begin
        out = simulate_individual(SYN, synthetic_forcing(730; temp_c=22.0, chl_ug_l=20.0),
                                  initial_state(SYN; shell_height_mm=40.0))
        spawn_days = findall(>(0), out.spawned_J)
        @test !isempty(spawn_days)
        @test all(out.E_R_J[spawn_days] .== 0)
        cold = simulate_individual(SYN, synthetic_forcing(730; temp_c=SYN.T_spawn_c, chl_ug_l=20.0),
                                   initial_state(SYN; shell_height_mm=40.0))
        @test all(==(0), cold.spawned_J)   # threshold must be exceeded, not met
    end

    @testset "low salinity stops feeding" begin
        u0 = initial_state(SYN; shell_height_mm=40.0)
        fresh = simulate_individual(SYN, synthetic_forcing(30; salinity_psu=1.0), u0)
        @test fresh.E_J[end] < u0.E
    end

    @testset "shell height and initial state" begin
        u0 = initial_state(SYN; shell_height_mm=55.0, e=0.5)
        @test shell_height_mm(u0.V, SYN.delta_M) ≈ 55.0
        @test u0.E ≈ 0.5 * SYN.E_m * u0.V
        @test dry_tissue_g(u0.E, u0.V, 0.0, SYN) ≈ SYN.d_V * u0.V + u0.E / SYN.mu_E
        @test gonado_somatic_index(0.0, 0.0, 0.0, SYN) == 0.0
        @test_throws DomainError initial_state(SYN; shell_height_mm=10.0, e=1.5)
    end

    @testset "runs on a forcing_daily contract file" begin
        df = read_contract(joinpath(@__DIR__, "..", "data", "fixtures", "synthetic_forcing_daily_valid.arrow"),
                           "forcing_daily")
        @test_throws ArgumentError simulate_individual(SYN, df, initial_state(SYN; shell_height_mm=20.0))
        df.chl_ug_l = coalesce.(df.chl_ug_l, 2.0)   # synthetic gap fill for the test
        out = simulate_individual(SYN, df, initial_state(SYN; shell_height_mm=20.0))
        @test names(out) == ["date", "E_J", "V_cm3", "E_R_J", "shell_height_mm", "dry_tissue_g", "gsi", "spawned_J"]
        @test nrow(out) == nrow(df)
        @test out.date == df.date
        @test_throws ArgumentError simulate_individual(SYN, df[[1, 3], :], initial_state(SYN; shell_height_mm=20.0))
    end

    @testset "parameter files" begin
        gig = DEBParams(load_params(joinpath(PARAMS_DIR, "deb_c_gigas_pouvreau2006.toml")))
        @test gig.p_Am ≈ 0.75 * 560
        @test isnan(gig.S_L)
        @test salinity_factor(1.0, gig) == 1.0
        # Eastern oyster set (Lavaud et al. 2017): complete except X_K, which comes from calibration
        vir = load_params(joinpath(PARAMS_DIR, "deb_c_virginica.toml"))
        @test !haskey(vir, "X_K")
        @test_throws KeyError DEBParams(vir)
        pv = DEBParams(vir; X_K=3.0)
        @test pv.p_Am ≈ 0.75 * 249.5                 # = 187.125, as in Lavaud et al. 2021 Table A2
        @test pv.p_Am / pv.E_m ≈ 0.03453 rtol = 1e-3 # energy conductance v of Lavaud et al. 2021
        @test pv.V_p ≈ (0.2 * 3.25)^3
        @test salinity_factor(6.5, pv) ≈ 0.5
        @test arrhenius_factor(293.0 - 273.15, pv) ≈ 1 atol = 1e-12
        @test_throws ArgumentError DEBParams(vir; X_K=-1.0)
        d2 = copy(vir); d2["X_K"] = Param(3.0, "ug chl-a L^-1", "synthetic: test", "")
        @test_throws ArgumentError DEBParams(d2; X_K=3.0)
        # E_m from energy conductance
        d = load_params(joinpath(@__DIR__, "fixtures", "synthetic_deb_params.toml"))
        delete!(d, "E_m")
        d["v"] = Param(0.1, "cm d^-1", "synthetic: test fixture", "")
        @test DEBParams(d).E_m ≈ SYN.p_Am / 0.1
        delete!(d, "S_H")
        @test_throws ArgumentError DEBParams(d)
    end
end
