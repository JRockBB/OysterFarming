# Perkinsus marinus (Dermo) within-host model (PLAN.md P4.1).
#
# Ragone Calvo, L.M., Wetzel, R.L., Burreson, E.M., 2001. Development and verification of a model
# for the population dynamics of the protistan parasite, Perkinsus marinus, within its host, the
# eastern oyster, Crassostrea virginica, in Chesapeake Bay. J. Shellfish Res. 20(1):231-241.
#
# State: average within-host parasite density P, cells per g wet tissue, daily time step.

"""
    DermoParams

Parameters of the Ragone Calvo et al. (2001) model, read with
`DermoParams(load_params("params/dermo_ragone_calvo2001.toml"))`.
"""
struct DermoParams
    A_p::Float64
    G_p::Float64
    div_intercept::Float64
    div_slope::Float64
    q10_a::Float64
    div_T_min_c::Float64
    mort_T_ref_c::Float64
    mort_T_min_c::Float64
    mort_rate_S_gt_12::Float64
    mort_rate_S_8_12::Float64
    mort_rate_S_3_8::Float64
    mort_rate_S_lt_3::Float64
    transmission_doy::Int
    transmission_S_gt_15::Float64
    transmission_S_12_15::Float64
    transmission_S_9_12::Float64
    transmission_S_6_9::Float64
    loss_rate::Float64
    loss_P_threshold::Float64
end

function DermoParams(params::AbstractDict{String,<:Param})
    pv(n) = param_value(params, n)
    names = fieldnames(DermoParams)
    vals = [n === :transmission_doy ? Int(pv(String(n))) : pv(String(n)) for n in names]
    p = DermoParams(vals...)
    p.G_p > p.A_p > 0 || throw(ArgumentError("need G_p > A_p > 0"))
    return p
end

"""
    dermo_division_rate(temp_c, salinity_psu, p) -> Float64

Specific parasite division rate, d⁻¹ (Ragone Calvo et al. 2001, eq. 3):

    μd = (−0.009 + 0.0035·S) · exp(a·(T − 20))      for T > 20 °C, else 0

with a = 0.06931 (Q10 = 2). Used as published: below S ≈ 2.6 the rate is
negative (the regression was fitted at sites with mean salinity 4.6-11.7).
"""
function dermo_division_rate(temp_c::Real, salinity_psu::Real, p::DermoParams)
    temp_c > p.div_T_min_c || return 0.0
    return (p.div_intercept + p.div_slope * salinity_psu) * exp(p.q10_a * (temp_c - p.div_T_min_c))
end

"""
    dermo_mortality_rate(temp_c, salinity_psu, p) -> Float64

Specific parasite mortality rate, d⁻¹ (Ragone Calvo et al. 2001, eq. 5 and Table 1):

    μm = μm_sal(S) · exp(a·(18 − T))      for 3 < T < 18 °C, else 0

μm_sal = 0.0175 (S > 12), 0.02 (8 ≤ S ≤ 12), 0.0175 (3 ≤ S < 8), 0.006 (S < 3).
"""
function dermo_mortality_rate(temp_c::Real, salinity_psu::Real, p::DermoParams)
    (p.mort_T_min_c < temp_c < p.mort_T_ref_c) || return 0.0
    base = salinity_psu > 12 ? p.mort_rate_S_gt_12 :
           salinity_psu >= 8 ? p.mort_rate_S_8_12 :
           salinity_psu >= 3 ? p.mort_rate_S_3_8 : p.mort_rate_S_lt_3
    return base * exp(p.q10_a * (p.mort_T_ref_c - temp_c))
end

"""
    dermo_transmission_cells(salinity_psu, p) -> Float64

Cells g⁻¹ added by the annual transmission event (Ragone Calvo et al. 2001, eq. 6):
1000 (S > 15), 500 (12–15), 100 (9–12), 50 (6–9), 0 (S < 6). Values on an edge go
to the higher bin (the paper leaves edges undefined).
"""
function dermo_transmission_cells(salinity_psu::Real, p::DermoParams)
    salinity_psu >= 15 && return p.transmission_S_gt_15
    salinity_psu >= 12 && return p.transmission_S_12_15
    salinity_psu >= 9 && return p.transmission_S_9_12
    salinity_psu >= 6 && return p.transmission_S_6_9
    return 0.0
end

"""
    dermo_step(P, temp_c, salinity_psu, doy, p) -> Float64

Advance parasite density one day (Ragone Calvo et al. 2001, eqs. 1-7). The model
was built in STELLA with a daily step; this uses the same forward step:

    P(t+1) = P + P_G + P_T − P_M − P_L
    P_G = P · μd · (1 − (P − A_p)/(G_p − A_p))        (eq. 2)
    P_T = transmission cells on day `transmission_doy`  (eq. 6)
    P_M = P · μm                                         (eq. 4)
    P_L = P · 0.5 when P > 100,000, else 0               (eq. 7)

Eq. 1 is printed as dP/dt = P(P_G + P_T − P_M − P_L), but eqs. 2, 4, 6 and 7
already include P, so the flows are summed (the dimensionally consistent reading).
Eq. 2 is used exactly as published, so for P < A_p the feedback term exceeds 1.
"""
function dermo_step(P::Real, temp_c::Real, salinity_psu::Real, doy::Integer, p::DermoParams)
    P_G = P * dermo_division_rate(temp_c, salinity_psu, p) * (1 - (P - p.A_p) / (p.G_p - p.A_p))
    P_T = doy == p.transmission_doy ? dermo_transmission_cells(salinity_psu, p) : 0.0
    P_M = P * dermo_mortality_rate(temp_c, salinity_psu, p)
    P_L = P > p.loss_P_threshold ? P * p.loss_rate : 0.0
    return max(0.0, P + P_G + P_T - P_M - P_L)
end

"""
    simulate_dermo(p, forcing, P0) -> DataFrame

Run the within-host model over consecutive forcing days (columns `date`,
`temp_c`, `salinity_psu`). Returns `date`, `P_cells_g` (density at the end of the
day) and `log10_P`.
"""
function simulate_dermo(p::DermoParams, forcing, P0::Real)
    df = DataFrame(forcing; copycols=false)
    nrow(df) > 0 || throw(ArgumentError("forcing is empty"))
    all(diff(df.date) .== Day(1)) || throw(ArgumentError("forcing dates must be consecutive days"))
    P0 >= 0 || throw(DomainError(P0, "initial density must be non-negative"))
    P = Float64(P0)
    out = zeros(nrow(df))
    for i in 1:nrow(df)
        P = dermo_step(P, df.temp_c[i], df.salinity_psu[i], dayofyear(df.date[i]), p)
        out[i] = P
    end
    return DataFrame(date=collect(df.date), P_cells_g=out, log10_P=log10.(max.(out, 1.0)))
end

"""
    rank_to_log10_cells(rank) -> Float64

Convert a Ray/Mackin infection rank (0-7 scale) to log10 parasite cells per g wet
tissue with the regression of Bushek et al. (1994), as used by Ragone Calvo et
al. (2001, p. 234):  y = 0.176 − 0.463·x + 0.205·x².
"""
rank_to_log10_cells(rank::Real) = 0.176 - 0.463 * rank + 0.205 * rank^2

"""
    cells_to_intensity(P) -> Float64

Mean infection intensity rank (0-7 scale) for a parasite density P (cells per g
wet tissue), the inverse of [`rank_to_log10_cells`](@ref) (Bushek et al. 1994,
as used by Ragone Calvo et al. 2001): solves 0.205·x² − 0.463·x + 0.176 = log10 P
for the root x ≥ 1.13 (where the quadratic is increasing). Densities at or below
the quadratic's minimum (≈ 0.92 cells/g) map to 0.

Assumes the Maryland DNR 0-7 intensity scale and the ranks used by Bushek et al.
are the same Ray-assay scale (recorded in PLAN.md P4.1).
"""
function cells_to_intensity(P::Real)
    P >= 0 || throw(DomainError(P, "density must be non-negative"))
    a, b0, c = 0.205, -0.463, 0.176
    xmin = -b0 / (2a)
    ymin = c - b0^2 / (4a)
    y = P > 0 ? log10(P) : -Inf
    y <= ymin && return 0.0
    return (-b0 + sqrt(b0^2 - 4a * (c - y))) / (2a)
end

"""
    dermo_excess_mortality_d(P, per_intensity_yr) -> Float64

Daily excess oyster mortality rate from Dermo for parasite density P:
`per_intensity_yr × cells_to_intensity(P) / 365`, with `per_intensity_yr` the
fitted annual rate per unit intensity (params/dermo_mortality_mddnr.toml).
Use as `extra_mortality_d` in [`simulate_farm`](@ref).
"""
dermo_excess_mortality_d(P::Real, per_intensity_yr::Real) = per_intensity_yr * cells_to_intensity(P) / 365
