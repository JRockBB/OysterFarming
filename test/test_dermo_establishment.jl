# Tests for src/dermo_establishment.jl (PLAN.md P4.2e). Forcing and parameter values below are SYNTHETIC
# except the lab anchors, which are read from params/dermo_establishment_lab.toml.

@testset "dermo_establishment" begin
    lab = load_params(joinpath(@__DIR__, "..", "params", "dermo_establishment_lab.toml"))
    pw = load_params(joinpath(@__DIR__, "..", "params", "dermo_powell1996.toml"))
    mk(; beta=1e-9, mu=0.1, epsilon=1e-4, host=1.0) = EstablishmentParams(lab, pw; beta=beta, mu=mu, epsilon=epsilon, host=host)
    p = mk()

    @testset "lab-anchored rates" begin
        for (t, r) in zip(lab["rate_anchor_temp_c"].value, lab["rate_anchor_net_d"].value)
            @test establishment_rate(t, 25.0, p) ≈ r
        end
        @test establishment_rate(2.0, 25.0, p) ≈ establishment_rate(6.0, 25.0, p)     # constant below 6 degC
        @test establishment_rate(32.0, 25.0, p) ≈ establishment_rate(28.0, 25.0, p)   # constant above 28 degC
        @test establishment_rate(28.0, 7.0, p) == 0.0                                  # salinity factor 0 at 7 ppt
        @test establishment_rate(28.0, 15.0, p) ≈ establishment_rate(28.0, 30.0, p)    # factor 1 at >= 15 ppt
        @test establishment_rate(28.0, 25.0, mk(host=0.5)) ≈ 0.5 * establishment_rate(28.0, 25.0, p)
        @test establishment_rate(4.0, 3.5, p) ≈ establishment_rate(4.0, 20.0, p) - lab["cold_low_salinity_loss_d"].value
        @test establishment_rate(4.0, 7.0, p) ≈ establishment_rate(4.0, 20.0, p)       # no extra loss at 7 ppt
    end

    @testset "Mackin conversion" begin
        @test est_mackin(1, p) == 0.0
        @test est_mackin(14, p) == 0.0                 # 2^12 cells: at the detection limit, scored negative
        @test est_mackin(EST_N_LEVELS, p) > 0
        @test issorted([est_mackin(k, p) for k in 1:EST_N_LEVELS])
    end

    nd = 365 * 3
    warm = fill(26.0, nd); cold = fill(4.0, nd); sal = fill(25.0, nd)
    N0 = zeros(EST_N_LEVELS); N0[1] = 20.0

    @testset "transmission switches" begin
        r = simulate_establishment(mk(beta=0.0, epsilon=0.0), warm, sal, N0, 0.0; record=[nd])
        @test r.true_prev[1] == 0.0 && sum(r.deaths) == 0.0
        late = simulate_establishment(mk(beta=0.0, epsilon=1e-3), warm, sal, N0, 0.0; intro_day=nd, record=[nd - 1, nd])
        @test late.true_prev[1] == 0.0 && late.true_prev[2] > 0.0                      # import starts on intro_day
    end

    @testset "conservation and warm-water deaths" begin
        r = simulate_establishment(p, warm, sal, N0, 0.0; record=[100, nd])
        @test sum(r.N) ≈ 20.0 rtol = 1e-10                                             # deaths replaced
        @test sum(r.deaths) > 0
        @test r.prev[2] > 0
        c = simulate_establishment(p, cold, sal, N0, 0.0; record=[nd])
        @test sum(c.deaths) == 0.0                                                     # no growth, no deaths in cold water
        @test c.prev[1] == 0.0                                                         # infections stay below detection
        @test_throws ArgumentError simulate_establishment(p, warm, sal[1:10], N0, 0.0)
        @test_throws ArgumentError simulate_establishment(p, warm, sal, zeros(5), 0.0)
    end
end
