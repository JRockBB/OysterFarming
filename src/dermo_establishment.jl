# Dermo establishment model (PLAN.md P4.2e): infection classes of market-size oysters, within-host parasite
# growth fixed from laboratory rates, and transmission through a pool of waterborne parasites shed by
# infected and dying oysters. Built because P4.2d showed that in cool, salty water (Great Bay, NH) prevalence
# follows establishment (invasion) dynamics that climate-only models cannot represent.
#
# State per site: N[k], oysters m^-2 at infection level k (k = 1 uninfected; level k holds 2^(k-2) cells),
# market size only (Powell et al. 1996 size class 7), and W, waterborne parasites (cells m^-2).
# Daily step:
#   within-host: net rate r(T,S) (lab anchors, params/dermo_establishment_lab.toml), positive rates scaled by
#     the fitted host factor h and limited at high density as in Hofmann et al. (1995, eq. 11). Oysters move up
#     one level at rate r and down at rate 2|r| (cell-conserving; Powell et al. 1996 p. 144). Level 2 never
#     moves down (infections are not lost). Entering the lethal level kills the oyster.
#   shedding: live infected oysters shed 10^(a) C^b cells d^-1 (Bushek et al. 2002); dead oysters release all
#     their cells. W decays at the fitted rate mu.
#   infection: uninfected oysters become infected (level 2) at rate beta * W + epsilon (fitted transmission
#     coefficient and background import).
#   host turnover: non-Dermo mortality m_b; every death is replaced by an uninfected oyster (constant density).

"""
    EstablishmentParams{T}

Lab-anchored constants (from `params/dermo_establishment_lab.toml`) and fitted parameters `beta`
(m² cell⁻¹ d⁻¹), `mu` (d⁻¹), `epsilon` (d⁻¹) and `host` (0-1].
"""
struct EstablishmentParams{T<:Real}
    anchor_T::Vector{Float64}
    anchor_r::Vector{Float64}
    anchor_S::Vector{Float64}
    anchor_fS::Vector{Float64}
    cold_lowS_loss::Float64
    shed_a::Float64
    shed_b::Float64
    detection_cells::Float64
    lethal_level::Int
    W_afdw::Float64
    dens_beta::Float64          # Hofmann et al. 1995 eq. 11 beta
    dens_gamma::Float64         # Hofmann et al. 1995 eq. 11 gamma
    upsilon::Float64            # Mackin conversion (Hofmann et al. 1995 eq. 23)
    tau::Float64
    phi::Float64
    beta::T
    mu::T
    epsilon::T
    host::T
end

const EST_N_LEVELS = 28

function EstablishmentParams(lab::AbstractDict{String,<:Param}, powell::AbstractDict{String,<:Param};
                             beta, mu, epsilon, host)
    v(n) = Float64.(lab[n].value)
    T = promote_type(typeof(beta), typeof(mu), typeof(epsilon), typeof(host))
    return EstablishmentParams{T}(v("rate_anchor_temp_c"), v("rate_anchor_net_d"), v("salinity_anchor_ppt"),
                                  v("salinity_anchor_factor"), param_value(lab, "cold_low_salinity_loss_d"),
                                  param_value(lab, "shed_log10_intercept"), param_value(lab, "shed_log10_slope"),
                                  param_value(lab, "detection_cells"), Int(lab["lethal_level"].value),
                                  param_value(lab, "host_afdw_g"), param_value(powell, "beta"),
                                  param_value(powell, "gamma"), param_value(powell, "upsilon"),
                                  param_value(powell, "tau"), param_value(powell, "phi"),
                                  T(beta), T(mu), T(epsilon), T(host))
end

"Piecewise-linear interpolation, constant beyond the end anchors."
function _interp(x, xs, ys)
    x <= xs[1] && return ys[1]
    x >= xs[end] && return ys[end]
    i = searchsortedlast(xs, x)
    return ys[i] + (ys[i+1] - ys[i]) * (x - xs[i]) / (xs[i+1] - xs[i])
end

"""
    establishment_rate(T, S, p) -> Real

Net in-host parasite rate (d⁻¹, natural log) before density limitation: lab rate at temperature T; if positive,
times the salinity factor and the host factor; if negative, an extra loss below 7 ppt (La Peyre et al. 2010).
"""
function establishment_rate(T::Real, S::Real, p::EstablishmentParams)
    r = _interp(T, p.anchor_T, p.anchor_r)
    if r > 0
        return r * _interp(S, p.anchor_S, p.anchor_fS) * p.host
    else
        return S < 7 ? r - p.cold_lowS_loss * min(1.0, (7 - S) / 3.5) : r
    end
end

est_cells(k::Integer) = k == 1 ? 0.0 : 2.0^(k - 2)

"Mackin rank 0-5 of level k (Hofmann et al. 1995 eq. 23; 0 at or below the detection limit)."
function est_mackin(k::Integer, p::EstablishmentParams)
    C = est_cells(k)
    C <= p.detection_cells && return 0.0
    return clamp(log10(C / (p.upsilon * p.tau * p.W_afdw)) / p.phi, 0.0, 5.0)
end

"""
    simulate_establishment(p, T, S, N0, W0; intro_day=1, background_mortality_d=0.0302/365, record=nothing,
                           transmission_mod=(T, S) -> 1.0, smooth=false) -> NamedTuple

Run daily over temperature and salinity vectors `T`, `S`. Background import `epsilon` acts from day
`intro_day` on. `transmission_mod(T, S)` (>= 0) multiplies the waterborne infection rate beta * W
(PLAN.md P4.2f: a learned or symbolic term); its output type joins the state's element type, so
ForwardDiff dual numbers in it propagate. With `smooth = true` the daily infection probability is
1 - exp(-rate) and the waterborne pool decays by exp(-mu) per day (exact daily conversions, smooth in the
parameters, for gradient-based fitting; PLAN.md P4.2g); the default keeps the earlier capped forms. `record` is a vector of day indices at which to return detectable prevalence (`prev`),
true prevalence (`true_prev`) and mean Mackin rank including zeros (`mackin`); `deaths` is the Dermo death
count per day (oysters m⁻²) and `density` the constant density.
"""
function simulate_establishment(p::EstablishmentParams{P}, T::AbstractVector, S::AbstractVector,
                                N0::AbstractVector, W0::Real; intro_day::Integer=1,
                                background_mortality_d::Real=0.0302 / 365, record=nothing,
                                transmission_mod=(T, S) -> 1.0, smooth::Bool=false) where {P}
    nd = length(T)
    nd > 0 || throw(ArgumentError("forcing is empty"))
    R = promote_type(P, typeof(transmission_mod(first(T), first(S))))
    length(S) == nd || throw(ArgumentError("T and S differ in length"))
    length(N0) == EST_N_LEVELS || throw(ArgumentError("N0 must have $(EST_N_LEVELS) levels"))
    N = R.(collect(N0)); W = R(W0)
    dens = sum(N)
    cells = [est_cells(k) for k in 1:EST_N_LEVELS]
    shed = [k == 1 ? 0.0 : 10.0^p.shed_a * cells[k]^p.shed_b for k in 1:EST_N_LEVELS]
    dlim = [k == 1 ? 1.0 : p.dens_beta * (cells[k] / p.W_afdw)^p.dens_gamma for k in 1:EST_N_LEVELS]
    detect = [cells[k] > p.detection_cells for k in 1:EST_N_LEVELS]
    mk = [est_mackin(k, p) for k in 1:EST_N_LEVELS]
    rec = record === nothing ? Int[] : collect(record)
    prev = zeros(R, length(rec)); tprev = zeros(R, length(rec)); mack = zeros(R, length(rec))
    deaths = zeros(R, nd)
    up = zeros(R, EST_N_LEVELS); down = zeros(R, EST_N_LEVELS)
    L = p.lethal_level
    ri = 1
    for d in 1:nd
        r0 = establishment_rate(T[d], S[d], p)
        for k in 2:EST_N_LEVELS
            g = r0 * min(one(R), dlim[k])
            up[k] = g > 0 ? min(one(R), g) * N[k] : zero(R)
            down[k] = (g < 0 && k > 2) ? min(one(R), 2 * abs(g)) * N[k] : zero(R)
        end
        lam = p.beta * transmission_mod(T[d], S[d]) * W + (d >= intro_day ? p.epsilon : zero(R))
        newinf = (smooth ? 1 - exp(-lam) : min(one(R), lam)) * N[1]
        dead = zero(R); released = zero(R)
        N[1] -= newinf; N[2] += newinf
        for k in 2:EST_N_LEVELS
            N[k] -= up[k] + down[k]
            down[k] > 0 && (N[k-1] += down[k])
            if k + 1 >= L || k == EST_N_LEVELS
                dead += up[k]; released += up[k] * cells[min(k + 1, EST_N_LEVELS)]
            else
                N[k+1] += up[k]
            end
        end
        shedding = zero(R)
        for k in 2:EST_N_LEVELS
            shedding += N[k] * shed[k]
        end
        W = W * (smooth ? exp(-p.mu) : 1 - min(one(R), p.mu)) + shedding + released
        # background mortality, then replace all deaths with uninfected oysters
        bg = zero(R)
        for k in 1:EST_N_LEVELS
            x = N[k] * background_mortality_d
            N[k] -= x; bg += x
        end
        N[1] += dead + bg
        deaths[d] = dead
        while ri <= length(rec) && rec[ri] == d
            tot = sum(N)
            prev[ri] = sum(N[k] for k in 1:EST_N_LEVELS if detect[k]) / tot
            tprev[ri] = (tot - N[1]) / tot
            mack[ri] = sum(N[k] * mk[k] for k in 1:EST_N_LEVELS) / tot
            ri += 1
        end
    end
    return (prev=prev, true_prev=tprev, mackin=mack, deaths=deaths, density=dens, N=N, W=W)
end
