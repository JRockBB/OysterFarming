# Standard Dynamic Energy Budget (DEB) model for one oyster.
#
# Structure follows Pouvreau et al. (2006) "Application of a dynamic energy
# budget model to the Pacific oyster, Crassostrea gigas, reared under various
# environmental conditions", J. Sea Res. 56:156-167 (eqs. 1-9), which applies
# Kooijman (2000) "Dynamic Energy and Mass Budgets in Biological Systems".
# The salinity correction on feeding follows Lavaud et al. (2017), Ecol. Model.
# 363:221-233, as reproduced in Lavaud et al. (2021), Estuar. Coast. Shelf Sci.,
# Table A1.
#
# State variables (per individual):
#   E    reserve energy, J
#   V    structural volume, cm^3
#   E_R  reproduction buffer energy, J
#
# Food X is chlorophyll-a in µg/L (the forcing_daily column chl_ug_l).

"""
    DEBParams

Parameter set for the standard DEB model. Build it from a parameter file with
`DEBParams(load_params(path))`; never construct it from guessed values. The
element type `T` is `Float64` normally; calibration replaces selected values
with other number types (e.g. ForwardDiff duals) via [`with_values`](@ref).

Fields (units):
- `p_Am`   maximum surface-area-specific assimilation rate, J cm⁻² d⁻¹
- `E_m`    maximum reserve density, J cm⁻³
- `p_M`    volume-specific somatic maintenance rate, J cm⁻³ d⁻¹
- `E_G`    volume-specific cost of structure, J cm⁻³
- `kappa`  fraction of mobilised reserve to maintenance + growth, –
- `kappa_R` fraction of reproduction buffer energy fixed in gametes, –
- `V_p`    structural volume at sexual maturity, cm³
- `delta_M` shape coefficient (V = (δ_M L)³ with L in cm), –
- `X_K`    half-saturation food density, µg chl-a L⁻¹
- `T_A`, `T_ref`, `T_L`, `T_H`, `T_AL`, `T_AH`  Arrhenius parameters, K
  (`T_L`…`T_AH` are `NaN` when the parameter file has no tolerance bounds)
- `mu_E`   energy content of reserve and buffer per dry mass, J g⁻¹
- `d_V`    dry mass of structure per structural volume, g cm⁻³
- `GSI_spawn` gonado-somatic index (fraction) that must be exceeded to spawn, –
- `T_spawn_c` water temperature that must be exceeded to spawn, °C
- `S_L`, `S_H` salinity thresholds of the feeding correction (`NaN` = no salinity effect)
"""
struct DEBParams{T<:Real}
    p_Am::T
    E_m::T
    p_M::T
    E_G::T
    kappa::T
    kappa_R::T
    V_p::T
    delta_M::T
    X_K::T
    T_A::T
    T_ref::T
    T_L::T
    T_H::T
    T_AL::T
    T_AH::T
    mu_E::T
    d_V::T
    GSI_spawn::T
    T_spawn_c::T
    S_L::T
    S_H::T
end

Base.eltype(::DEBParams{T}) where {T} = T

_optional(params, name) = haskey(params, name) ? param_value(params, name) : NaN

"""
    DEBParams(params::AbstractDict{String,<:Param})

Build a [`DEBParams`](@ref) from a loaded parameter file.

- `p_Am` is read directly, or computed as `ae * p_Xm` (Pouvreau et al. 2006,
  eq. 2: {p_Am} = ae·µ_X·{J_Xm}, with {p_Xm} already in energy units).
- `E_m` is read directly, or computed as `p_Am / v` from the energy
  conductance `v` (cm d⁻¹), since v = {p_Am}/[E_m] (Kooijman 2010, Table 2.4).
- Temperature-tolerance bounds (`T_L`, `T_H`, `T_AL`, `T_AH`) must be all
  present or all absent. Salinity thresholds (`S_L`, `S_H`) likewise.
- `X_K` may be omitted from the file and given as a keyword instead, e.g. as a
  starting value for calibration (PLAN.md P2.6). Giving it in both places is an error.
"""
function DEBParams(params::AbstractDict{String,<:Param}; X_K::Union{Nothing,Real}=nothing)
    pv(n) = param_value(params, n)
    p_Am = haskey(params, "p_Am") ? pv("p_Am") : pv("ae") * pv("p_Xm")
    E_m = haskey(params, "E_m") ? pv("E_m") : p_Am / pv("v")
    bounds = [_optional(params, n) for n in ("T_L", "T_H", "T_AL", "T_AH")]
    all(isnan, bounds) || !any(isnan, bounds) ||
        throw(ArgumentError("T_L, T_H, T_AL and T_AH must be given together"))
    sal = [_optional(params, n) for n in ("S_L", "S_H")]
    all(isnan, sal) || !any(isnan, sal) ||
        throw(ArgumentError("S_L and S_H must be given together"))
    if X_K === nothing
        X_K = pv("X_K")
    elseif haskey(params, "X_K")
        throw(ArgumentError("X_K is in the parameter file and also given as a keyword"))
    end
    X_K > 0 || throw(ArgumentError("X_K must be positive"))
    p = DEBParams{Float64}(p_Am, E_m, pv("p_M"), pv("E_G"), pv("kappa"), pv("kappa_R"), pv("V_p"),
                  pv("delta_M"), Float64(X_K), pv("T_A"), pv("T_ref"), bounds...,
                  pv("mu_E"), pv("d_V"), pv("GSI_spawn"), pv("T_spawn_c"), sal...)
    0 < p.kappa < 1 || throw(ArgumentError("kappa must be in (0, 1)"))
    0 < p.kappa_R <= 1 || throw(ArgumentError("kappa_R must be in (0, 1]"))
    return p
end

"""
    with_values(p::DEBParams; kwargs...) -> DEBParams

Copy of `p` with the named fields replaced, e.g. `with_values(p; X_K=4.0)`.
The element type is promoted to fit the new values, so automatic
differentiation through the model works. Used by calibration (P2.5).
"""
function with_values(p::DEBParams; kwargs...)
    names = fieldnames(DEBParams)
    for k in keys(kwargs)
        k in names || throw(ArgumentError("DEBParams has no field $k"))
    end
    vals = [haskey(kwargs, n) ? kwargs[n] : getfield(p, n) for n in names]
    T = promote_type(map(typeof, vals)...)
    return DEBParams{T}(T.(vals)...)
end

"""
    arrhenius_factor(temp_c, p::DEBParams; normalize=true) -> Float64

Temperature correction factor c_T applied to physiological rates. Temperature
is converted from °C to K here and only here (T = temp_c + 273.15).

Without tolerance bounds (Kooijman 2000; Pouvreau et al. 2006, p. 3):

    c_T = exp(T_A/T_ref − T_A/T)

With tolerance bounds T_L, T_H (Pouvreau et al. 2006, eq. 9):

    c_T = exp(T_A/T_ref − T_A/T) · s(T),   s(T) = [1 + exp(T_AL/T − T_AL/T_L) + exp(T_AH/T_H − T_AH/T)]⁻¹

With `normalize=true` (default) the bounded form is divided by s(T_ref) so
that c_T(T_ref) = 1 exactly, as in Kooijman (2010, eq. 1.3) and Lavaud et al.
(2017; Table A1 of Lavaud et al. 2021). `normalize=false` reproduces Pouvreau et
al. (2006) eq. 9 literally, which differs by the constant factor s(T_ref).
"""
function arrhenius_factor(temp_c::Real, p::DEBParams; normalize::Bool=true)
    T = temp_c + 273.15
    base = exp(p.T_A / p.T_ref - p.T_A / T)
    isnan(p.T_L) && return base
    s(x) = 1 / (1 + exp(p.T_AL / x - p.T_AL / p.T_L) + exp(p.T_AH / p.T_H - p.T_AH / x))
    return normalize ? base * s(T) / s(p.T_ref) : base * s(T)
end

"""
    functional_response(food, X_K) -> Float64

Holling type II scaled functional response (Pouvreau et al. 2006, eq. 1):

    f = X / (X + X_K)

`food` and `X_K` in the same units (here µg chl-a L⁻¹). Negative food is an error.
"""
function functional_response(food::Real, X_K::Real)
    food >= 0 || throw(DomainError(food, "food density must be non-negative"))
    return food / (food + X_K)
end

"""
    salinity_factor(salinity_psu, p::DEBParams) -> Float64

Salinity correction of feeding (Lavaud et al. 2017; reproduced in Lavaud et al.
2021, Table A1):

    c_S = 1                        if S ≥ S_H
    c_S = (S − S_L) / (S_H − S_L)  if S_L ≤ S < S_H
    c_S = 0                        if S < S_L

Returns 1 when the parameter set has no salinity thresholds (e.g. Pouvreau et
al. 2006, which has no salinity effect).
"""
function salinity_factor(salinity_psu::Real, p::DEBParams)
    isnan(p.S_L) && return one(p.S_L)
    salinity_psu >= p.S_H && return one(p.S_L)
    salinity_psu < p.S_L && return zero(p.S_L)
    return (salinity_psu - p.S_L) / (p.S_H - p.S_L)
end

"""
    deb_fluxes(E, V, E_R, f, c_T, c_S, p) -> NamedTuple

Energy fluxes (J d⁻¹) and state derivatives for one oyster.

Pouvreau et al. (2006) eqs. 2–7:

    p_A = c_S · c_T · {p_Am} · f · V^(2/3)
    p_M = c_T · [p_M] · V
    p_C = [E]/([E_G] + κ[E]) · ([E_G] · c_T{p_Am}/[E_m] · V^(2/3) + p_M),   [E] = E/V
    p_J = (1−κ)/κ · c_T · [p_M] · min(V, V_p)
    dE/dt   = p_A − p_C
    dV/dt   = (κ p_C − p_M) / [E_G]
    dE_R/dt = (1−κ) p_C − p_J                       (only once V ≥ V_p)

Before maturity (V < V_p) the flux (1−κ)p_C − p_J goes to maturation and is
dissipated, so dE_R/dt = 0.

Starvation (Pouvreau et al. 2006, section 2.1): when κ p_C < p_M, growth stops
and the somatic maintenance deficit is paid from E_R; once E_R is empty it is
paid by shrinking structure. Two implementation choices are not specified by
the paper and are made here: (1) structure is lysed at [E_G] J per cm³, and
(2) with an empty buffer, a shortfall in maturity maintenance is left unpaid.
"""
function deb_fluxes(E, V, E_R, f, c_T, c_S, p::DEBParams)
    V23 = V^(2 / 3)
    p_A = c_S * c_T * p.p_Am * f * V23
    p_Mflux = c_T * p.p_M * V
    Ed = V > 0 ? E / V : zero(V)
    p_C = Ed / (p.E_G + p.kappa * Ed) * (p.E_G * c_T * p.p_Am / p.E_m * V23 + p_Mflux)
    p_J = (1 - p.kappa) / p.kappa * c_T * p.p_M * min(V, p.V_p)
    somatic = p.kappa * p_C - p_Mflux
    z = zero(p_C)
    to_buffer = V >= p.V_p ? (1 - p.kappa) * p_C - p_J : z
    deficit = max(z, -somatic)
    dE_R = to_buffer - deficit
    dV = max(somatic, z) / p.E_G
    if E_R <= 0 && dE_R < 0
        dE_R = z
        dV -= deficit / p.E_G
    end
    return (p_A=p_A, p_C=p_C, p_M=p_Mflux, p_J=p_J, dE=p_A - p_C, dV=dV, dE_R=dE_R)
end

function _deb_rhs!(du, u, env, t)
    z = zero(eltype(u))
    E, V, E_R = max(u[1], z), max(u[2], z), max(u[3], z)
    fl = deb_fluxes(E, V, E_R, env.f, env.c_T, env.c_S, env.p)
    du[1] = fl.dE
    du[2] = fl.dV
    du[3] = fl.dE_R
    return nothing
end

"""
    shell_height_mm(V, delta_M) -> Float64

Shell length from structural volume (Pouvreau et al. 2006, section 2.3;
Bourlès et al. 2009, eq. 8): L = V^(1/3) / δ_M with L in cm, returned in mm.
"Shell length" in the French C. gigas literature is the longest shell
dimension, which is "shell height" in North American usage.
"""
shell_height_mm(V::Real, delta_M::Real) = 10 * cbrt(V) / delta_M

"""
    dry_tissue_g(E, V, E_R, p::DEBParams) -> Float64

Total soft-tissue dry mass, g:

    DW = d_V · V + (E + E_R) / µ_E

Structure converted with d_V (g cm⁻³) and reserve plus buffer with µ_E (J g⁻¹),
following the conversions in Pouvreau et al. (2006), section 2.3.
"""
dry_tissue_g(E::Real, V::Real, E_R::Real, p::DEBParams) = p.d_V * V + (E + E_R) / p.mu_E

"""
    gonado_somatic_index(E, V, E_R, p::DEBParams) -> Float64

Gonad mass over total flesh mass (Bourlès et al. 2009, section 2.3):
GSI = (E_R/µ_E) / DW. Returns 0 for an empty animal.
"""
function gonado_somatic_index(E::Real, V::Real, E_R::Real, p::DEBParams)
    dw = dry_tissue_g(E, V, E_R, p)
    return dw > 0 ? (E_R / p.mu_E) / dw : zero(dw)
end

"""
    spawns(temp_c, gsi, p::DEBParams) -> Bool

Spawning rule (Pouvreau et al. 2006, section 2.1): spawning occurs when the
gonado-somatic index exceeds `GSI_spawn` and water temperature exceeds `T_spawn_c`.
The buffer is then emptied and κ_R·E_R is released as gametes.
"""
spawns(temp_c::Real, gsi::Real, p::DEBParams) = temp_c > p.T_spawn_c && gsi > p.GSI_spawn

"""
    initial_state(p::DEBParams; shell_height_mm, e=1.0, E_R=0.0) -> NamedTuple

State `(E, V, E_R)` for an oyster of the given shell height with scaled reserve
density `e = [E]/[E_m]` (0–1). V = (δ_M L)³ with L in cm (Pouvreau et al. 2006,
section 2.3). The choice of `e` and `E_R` is the caller's; there is no default
initial condition in the literature.
"""
function initial_state(p::DEBParams; shell_height_mm::Real, e::Real=1.0, E_R::Real=0.0)
    0 <= e <= 1 || throw(DomainError(e, "scaled reserve density e must be in [0, 1]"))
    V = (p.delta_M * shell_height_mm / 10)^3
    return (E=e * p.E_m * V, V=V, E_R=float(E_R))
end

"""
    simulate_individual(p::DEBParams, forcing, u0; food_col=:chl_ug_l,
                        normalize_temperature=true, reltol=1e-8, abstol=1e-8) -> DataFrame

Run the DEB model for one oyster over a daily forcing table (one site, columns
`date`, `temp_c`, `salinity_psu` and `food_col`, as in `forcing_daily`). Forcing
is held constant within each day; each day is integrated with Tsit5, then the
spawning rule is applied. Dates must be consecutive and food must not be
missing (gap-fill first, PLAN.md P1.4).

`u0` is a NamedTuple `(E, V, E_R)` such as from [`initial_state`](@ref).

Returns one row per day with the state at the end of that day: `date`, `E_J`,
`V_cm3`, `E_R_J`, `shell_height_mm`, `dry_tissue_g`, `gsi`, `spawned_J`.
"""
function simulate_individual(p::DEBParams, forcing, u0; food_col::Symbol=:chl_ug_l,
                             normalize_temperature::Bool=true, reltol=1e-8, abstol=1e-8)
    df = DataFrame(forcing; copycols=false)
    n = nrow(df)
    n > 0 || throw(ArgumentError("forcing is empty"))
    for c in (:date, :temp_c, :salinity_psu, food_col)
        hasproperty(df, c) || throw(ArgumentError("forcing has no column $c"))
    end
    length(unique(df.date)) == n && all(diff(df.date) .== Day(1)) ||
        throw(ArgumentError("forcing dates must be consecutive days for a single site"))
    any(ismissing, df[!, food_col]) &&
        throw(ArgumentError("food column $food_col has missing values; gap-fill the forcing first"))
    all(>=(0), (u0.E, u0.V, u0.E_R)) || throw(DomainError(u0, "initial state must be non-negative"))

    temps = Vector{Float64}(df.temp_c)
    sal = Vector{Float64}(df.salinity_psu)
    food = Vector{Float64}(df[!, food_col])
    cols = _simulate_days(p, temps, sal, food, u0, normalize_temperature, reltol, abstol)
    return DataFrame(date=collect(df.date), E_J=cols.E, V_cm3=cols.V, E_R_J=cols.E_R,
                     shell_height_mm=cols.shell_height_mm, dry_tissue_g=cols.dry_tissue_g,
                     gsi=cols.gsi, spawned_J=cols.spawned)
end

"""
    IndividualStepper(p, u0; normalize_temperature=true, reltol=1e-8, abstol=1e-8)

Advances one oyster's DEB state a day at a time with a single reused Tsit5
integrator. Forcing for each day is passed to [`advance_day!`](@ref). Used by
[`simulate_individual`](@ref) and by the farm layer, so both share the same
daily step and spawning rule.
"""
mutable struct IndividualStepper{P<:DEBParams,T<:Real,I}
    p::P
    u::Vector{T}
    integ::I
    normalize_temperature::Bool
end

function IndividualStepper(p::DEBParams, u0; normalize_temperature::Bool=true, reltol=1e-8, abstol=1e-8)
    all(>=(0), (u0.E, u0.V, u0.E_R)) || throw(DomainError(u0, "initial state must be non-negative"))
    T = promote_type(eltype(p), typeof(u0.E), typeof(u0.V), typeof(u0.E_R))
    u = T[u0.E, u0.V, u0.E_R]
    env = (p=p, f=zero(T), c_T=one(T), c_S=one(T))
    integ = init(ODEProblem(_deb_rhs!, copy(u), (0.0, Inf), env), Tsit5();
                 reltol=reltol, abstol=abstol, save_everystep=false, save_start=false, save_end=false)
    return IndividualStepper(p, u, integ, normalize_temperature)
end

"""
    advance_day!(s::IndividualStepper, temp_c, salinity_psu, f) -> NamedTuple

Integrate one day with constant temperature, salinity and scaled functional
response `f` (0–1), then apply the spawning rule ([`spawns`](@ref)). Returns
`(E, V, E_R, gsi, spawned)` at the end of the day. Forcing changes at the day
boundary by swapping the parameter tuple, which makes Tsit5 re-evaluate the
derivative.
"""
function advance_day!(s::IndividualStepper{P,T}, temp_c::Real, salinity_psu::Real, f::Real) where {P,T}
    p = s.p
    s.integ.p = (p=p, f=T(f), c_T=T(arrhenius_factor(temp_c, p; normalize=s.normalize_temperature)),
                 c_S=T(salinity_factor(salinity_psu, p)))
    derivative_discontinuity!(s.integ, true)
    step!(s.integ, 1.0, true)
    u = s.u
    u .= max.(s.integ.u, zero(T))
    gsi = gonado_somatic_index(u[1], u[2], u[3], p)
    spawned = zero(T)
    if spawns(temp_c, gsi, p)
        spawned = p.kappa_R * u[3]
        u[3] = zero(T)
        gsi = zero(T)
    end
    set_u!(s.integ, u)
    return (E=u[1], V=u[2], E_R=u[3], gsi=gsi, spawned=spawned)
end

# Type-stable inner loop (function barrier) for simulate_individual.
function _simulate_days(p::DEBParams, temps::Vector{Float64}, sal::Vector{Float64}, food::Vector{Float64},
                        u0, normalize_temperature::Bool, reltol, abstol)
    n = length(temps)
    s = IndividualStepper(p, u0; normalize_temperature=normalize_temperature, reltol=reltol, abstol=abstol)
    T = eltype(s.u)
    E_out, V_out, ER_out = zeros(T, n), zeros(T, n), zeros(T, n)
    sh_out, dw_out, gsi_out, sp_out = zeros(T, n), zeros(T, n), zeros(T, n), zeros(T, n)
    for i in 1:n
        r = advance_day!(s, temps[i], sal[i], functional_response(food[i], p.X_K))
        E_out[i], V_out[i], ER_out[i] = r.E, r.V, r.E_R
        sh_out[i] = shell_height_mm(r.V, p.delta_M)
        dw_out[i] = dry_tissue_g(r.E, r.V, r.E_R, p)
        gsi_out[i] = r.gsi
        sp_out[i] = r.spawned
    end
    return (E=E_out, V=V_out, E_R=ER_out, shell_height_mm=sh_out, dry_tissue_g=dw_out, gsi=gsi_out,
            spawned=sp_out)
end
