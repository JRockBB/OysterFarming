# Nitrogen removed with harvested oysters (PLAN.md P3.4).

"""
    NitrogenParams

Shell-height-to-dry-weight regressions and nitrogen fractions, read with
`NitrogenParams(load_params("params/nitrogen_chesapeake_bmp.toml"))`.
"""
struct NitrogenParams
    tissue_a_diploid::Float64
    tissue_b_diploid::Float64
    tissue_a_triploid::Float64
    tissue_b_triploid::Float64
    tissue_n_fraction::Float64
    shell_a::Float64
    shell_b::Float64
    shell_n_fraction::Float64
end

function NitrogenParams(params::AbstractDict{String,<:Param})
    pv(n) = param_value(params, n)
    return NitrogenParams(pv("tissue_dw_a_diploid"), pv("tissue_dw_b_diploid"), pv("tissue_dw_a_triploid"),
                          pv("tissue_dw_b_triploid"), pv("tissue_n_fraction"), pv("shell_dw_a"),
                          pv("shell_dw_b"), pv("shell_n_fraction"))
end

"""
    nitrogen_per_oyster_g(shell_height_mm, np; ploidy=:diploid, include_shell=false) -> Float64

Nitrogen (g) in one harvested oyster (Oyster BMP Expert Panel 2016, Table 2b;
second report 2023, section 7.2.8 and Table 7.2):

    tissue N = a_t · SH^b_t · n_tissue            (a_t, b_t by ploidy)
    shell N  = a_s · SH^b_s · n_shell             (only if include_shell)

with SH the shell height in mm and dry weights in g.
"""
function nitrogen_per_oyster_g(shell_height_mm::Real, np::NitrogenParams; ploidy::Symbol=:diploid,
                               include_shell::Bool=false)
    shell_height_mm >= 0 || throw(DomainError(shell_height_mm, "shell height must be non-negative"))
    a, b = if ploidy === :diploid
        np.tissue_a_diploid, np.tissue_b_diploid
    elseif ploidy === :triploid
        np.tissue_a_triploid, np.tissue_b_triploid
    else
        throw(ArgumentError("ploidy must be :diploid or :triploid"))
    end
    n = a * shell_height_mm^b * np.tissue_n_fraction
    include_shell && (n += np.shell_a * shell_height_mm^np.shell_b * np.shell_n_fraction)
    return n
end

"""
    nitrogen_removed_g(harvests, np; kwargs...) -> Float64

Total nitrogen (g) in a harvest log from [`simulate_farm`](@ref): the sum over
rows of `count × nitrogen_per_oyster_g(shell_height_mm)`. Keyword arguments are
passed to [`nitrogen_per_oyster_g`](@ref).
"""
nitrogen_removed_g(harvests, np::NitrogenParams; kwargs...) =
    sum((r.count * nitrogen_per_oyster_g(r.shell_height_mm, np; kwargs...) for r in eachrow(harvests)); init=0.0)
