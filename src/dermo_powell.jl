# Perkinsus marinus (Dermo) infection-class model of Hofmann et al. (1995) and Powell et al. (1996)
# (PLAN.md P4.1b). A second Dermo model; the Ragone Calvo et al. (2001) model in dermo.jl is kept.
#
# Hofmann, E.E., Powell, E.N., Klinck, J.M., Saunders, G. (1995) J. Shellfish Res. 14(1):121-151 (Part I).
# Powell, E.N., Klinck, J.M., Hofmann, E.E. (1996) J. Shellfish Res. 15(1):141-165 (Part II).
#
# Reduced form (a stated simplification): the parasite and transmission parts are implemented as
# published; the host energy budget is not. Host size class weights are fixed, so the parasite's energy
# drain and filtration reduction do not feed back on host growth. The caller supplies host growth
# between size classes, recruitment and non-Dermo mortality. State: numbers of oysters per m^2 by size
# class j (11 classes, Part II Table 1) and infection level k (k = 1 uninfected; level k holds
# 2^(k-2) cells per oyster; 28 levels).

"""
    PowellDermoParams

Parameters of the Hofmann et al. (1995) / Powell et al. (1996) parasite model, read with
`PowellDermoParams(load_params("params/dermo_powell1996.toml"))`. `halving_time_h` sets the parasite
mortality rate at 20 °C and 20 ppt: r_m0 = ln 2 · 24 / halving_time_h.
"""
struct PowellDermoParams
    r_d0::Float64
    r_d1::Float64
    alpha::Float64
    T0::Float64
    S0::Float64
    beta::Float64
    gamma::Float64
    r_m0::Float64
    delta::Float64
    upsilon::Float64
    tau::Float64
    phi::Float64
    r_tb::Float64
    r_t0::Float64
    infective_load::Float64
    oyster_density_sat::Float64
    detection_cells::Float64
    W::Vector{Float64}               # class weight (band midpoint), g AFDW
    lethal_level::Vector{Int}        # infection level index at which oysters die (doublings + 2)
end

const POWELL_N_LEVELS = 28

function PowellDermoParams(params::AbstractDict{String,<:Param}; halving_time_h=nothing)
    pv(n) = param_value(params, n)
    bounds = Float64.(params["size_class_afdw_g"].value)
    lethal = Int.(params["lethal_doublings"].value)
    length(bounds) == length(lethal) + 1 || throw(ArgumentError("size class bounds and lethal levels differ"))
    W = [(bounds[i] + bounds[i+1]) / 2 for i in 1:length(lethal)]
    h = halving_time_h === nothing ? pv("halving_time_h") : Float64(halving_time_h)
    h > 0 || throw(ArgumentError("halving time must be positive"))
    return PowellDermoParams(pv("r_d0"), pv("r_d1"), pv("alpha"), pv("T0"), pv("S0"), pv("beta"), pv("gamma"),
                             log(2) * 24 / h, pv("delta"), pv("upsilon"), pv("tau"), pv("phi"), pv("r_tb"),
                             pv("r_t0"), pv("infective_load"), pv("oyster_density_sat"), pv("detection_cells"),
                             W, lethal .+ 2)
end

"Cells per oyster at infection level k (level 1 uninfected): 2^(k-2)."
powell_cells(k::Integer) = k == 1 ? 0.0 : 2.0^(k - 2)

"""
    powell_division_rate(T, S, p) -> Float64

Specific parasite division rate r_d(T, S), d⁻¹ (Hofmann et al. 1995, eqs. 9-10):
r_d0·exp(α(T − T0)) for S ≥ 10 ppt, r_d1·(S/10)·exp(α(T − T0)) below.
"""
function powell_division_rate(T::Real, S::Real, p::PowellDermoParams)
    f = exp(p.alpha * (T - p.T0))
    return S >= 10 ? p.r_d0 * f : p.r_d1 * max(S, 0.0) / 10 * f
end

"""
    powell_mortality_rate(T, S, p) -> Float64

Specific parasite mortality rate r_m(T, S), d⁻¹ (Hofmann et al. 1995, eqs. 12-13, in the reading
that reproduces their Fig. 2b): r_m0·exp(10δ)·exp(−δ(S/S0)(max(T, 10) − 10)).
"""
powell_mortality_rate(T::Real, S::Real, p::PowellDermoParams) =
    p.r_m0 * exp(10 * p.delta) * exp(-p.delta * (max(S, 0.0) / p.S0) * (max(T, 10.0) - 10))

"""
    powell_net_rate(T, S, k, j, p) -> Float64

Net specific rate of change of the parasite population in an oyster of size class j at infection
level k: min(r_d, β r_d (C/W)^γ) − min(r_m, β r_m (C/W)^γ) (Hofmann et al. 1995, eqs. 11, 13, 15).
"""
function powell_net_rate(T::Real, S::Real, k::Integer, j::Integer, p::PowellDermoParams)
    k == 1 && return 0.0
    dens = p.beta * (powell_cells(k) / p.W[j])^p.gamma
    rd = powell_division_rate(T, S, p)
    rm = powell_mortality_rate(T, S, p)
    return min(rd, dens * rd) - min(rm, dens * rm)
end

"""
    powell_mackin(k, j, p) -> Float64

Mackin rank (0-5) of an oyster in size class j at level k (Choi et al. 1989 as used by Hofmann
et al. 1995, eq. 23): log10(C / (υ τ W_j)) / φ, clipped to [0, 5]; 0 when C is at or below the
assay detection limit (`detection_cells`).
"""
function powell_mackin(k::Integer, j::Integer, p::PowellDermoParams)
    C = powell_cells(k)
    C <= p.detection_cells && return 0.0
    return clamp(log10(C / (p.upsilon * p.tau * p.W[j])) / p.phi, 0.0, 5.0)
end

"""
    powell_transmission_rate(N, p) -> Float64

Rate (d⁻¹) at which uninfected oysters become infected (Hofmann et al. 1995, eqs. 24-28):
(r_tb + r_t0 (P1 + P2 + P3)/3) / 365, with P1 = prevalence, P2 = min(1, total cells / total wet
weight / IL), P3 = min(1, total density / OD).
"""
function powell_transmission_rate(N::AbstractMatrix, p::PowellDermoParams)
    total = sum(N)
    total > 0 || return p.r_tb / 365
    P1 = sum(@view N[:, 2:end]) / total
    cells = sum(N[j, k] * powell_cells(k) for j in axes(N, 1), k in 2:size(N, 2))
    wet = sum(N[j, k] * p.upsilon * p.W[j] for j in axes(N, 1), k in axes(N, 2))
    P2 = min(1.0, cells / wet / p.infective_load)
    P3 = min(1.0, total / p.oyster_density_sat)
    return (p.r_tb + p.r_t0 * (P1 + P2 + P3) / 3) / 365
end

"""
    powell_step!(N, T, S, p; dt=1.0, nsub=4) -> Float64

Advance the infection-class state `N` (oysters m⁻², size classes × 28 levels) by `dt` days with
`nsub` explicit sub-steps. Per sub-step: uninfected oysters become infected at the transmission
rate; at each level k ≥ 2 the net parasite rate g moves oysters up one level at rate g (g > 0) or down
at rate 2|g| (g < 0), the cell-conserving scaling of Powell et al. 1996 (p. 144); level 2 never moves
down (infections are never lost); oysters entering their size class's lethal level die. Returns the
number of Dermo deaths (oysters m⁻²).
"""
function powell_step!(N::AbstractMatrix, T::Real, S::Real, p::PowellDermoParams; dt::Real=1.0, nsub::Integer=4)
    nj, nk = size(N)
    nk == POWELL_N_LEVELS || throw(ArgumentError("state must have $(POWELL_N_LEVELS) infection levels"))
    h = dt / nsub
    deaths = 0.0
    up = zeros(nj, nk); down = zeros(nj, nk)
    for _ in 1:nsub
        rt = powell_transmission_rate(N, p)
        fill!(up, 0.0); fill!(down, 0.0)
        for j in 1:nj, k in 2:nk
            g = powell_net_rate(T, S, k, j, p)
            if g > 0
                up[j, k] = min(1.0, g * h) * N[j, k]
            elseif g < 0 && k > 2
                down[j, k] = min(1.0, 2 * abs(g) * h) * N[j, k]
            end
        end
        for j in 1:nj
            newinf = min(1.0, rt * h) * N[j, 1]
            N[j, 1] -= newinf
            N[j, 2] += newinf
            L = p.lethal_level[j]
            for k in 2:nk
                N[j, k] -= up[j, k] + down[j, k]
                down[j, k] > 0 && (N[j, k-1] += down[j, k])
                if k + 1 >= L || k == nk
                    deaths += up[j, k]
                else
                    N[j, k+1] += up[j, k]
                end
            end
            # oysters already at or beyond the lethal level (e.g. from initial conditions) die
            for k in L:nk
                deaths += N[j, k]
                N[j, k] = 0.0
            end
        end
    end
    return deaths
end

"""
    powell_summary(N, p; classes=axes(N, 1)) -> NamedTuple

Field-comparable summaries over the given size classes: `prevalence` (fraction with more than
`detection_cells` cells, as scored by the thioglycollate assay), `true_prevalence` (any cells),
`mean_mackin` (mean Mackin rank 0-5 including zeros) and `density` (oysters m⁻²).
"""
function powell_summary(N::AbstractMatrix, p::PowellDermoParams; classes=axes(N, 1))
    tot = sum(N[j, k] for j in classes, k in axes(N, 2))
    tot > 0 || return (prevalence=NaN, true_prevalence=NaN, mean_mackin=NaN, density=0.0)
    det = sum(N[j, k] for j in classes, k in axes(N, 2) if powell_cells(k) > p.detection_cells; init=0.0)
    inf = sum(N[j, k] for j in classes, k in 2:size(N, 2))
    mk = sum(N[j, k] * powell_mackin(k, j, p) for j in classes, k in axes(N, 2))
    return (prevalence=det / tot, true_prevalence=inf / tot, mean_mackin=mk / tot, density=tot)
end

"""
    simulate_powell_dermo(p, forcing, N0; background_mortality_d=0.0, replace_dead=true, nsub=4) -> (DataFrame, N)

Run the infection-class model over consecutive forcing days (`date`, `temp_c`, `salinity_psu`) from
state `N0` (copied). Each day: [`powell_step!`](@ref), then non-Dermo mortality at
`background_mortality_d` across all classes and levels. With `replace_dead = true` every death
(Dermo and background) is replaced by an uninfected oyster in the same size class, holding density
constant (a stand-in for recruitment and growth into the class; a stated assumption, not from the
papers). Returns daily summaries over all classes (`powell_summary`) plus `dermo_deaths`, and the
final state.
"""
function simulate_powell_dermo(p::PowellDermoParams, forcing, N0::AbstractMatrix;
                               background_mortality_d::Real=0.0, replace_dead::Bool=true, nsub::Integer=4)
    df = DataFrame(forcing; copycols=false)
    nrow(df) > 0 || throw(ArgumentError("forcing is empty"))
    all(diff(df.date) .== Day(1)) || throw(ArgumentError("forcing dates must be consecutive days"))
    size(N0) == (length(p.W), POWELL_N_LEVELS) ||
        throw(ArgumentError("N0 must be $(length(p.W)) × $(POWELL_N_LEVELS)"))
    N = Matrix{Float64}(N0)
    out = DataFrame(date=Date[], prevalence=Float64[], true_prevalence=Float64[], mean_mackin=Float64[],
                    density=Float64[], dermo_deaths=Float64[])
    for i in 1:nrow(df)
        before = vec(sum(N; dims=2))
        dd = powell_step!(N, df.temp_c[i], df.salinity_psu[i], p; nsub=nsub)
        if background_mortality_d > 0
            N .*= (1 - background_mortality_d)
        end
        if replace_dead
            N[:, 1] .+= before .- vec(sum(N; dims=2))
        end
        s = powell_summary(N, p)
        push!(out, (df.date[i], s.prevalence, s.true_prevalence, s.mean_mackin, s.density, dd))
    end
    return out, N
end
