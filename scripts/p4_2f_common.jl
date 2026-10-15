# Shared definitions for the P4.2f scripts that read the trained UDE (identical to those in scripts/p4_2f_ude_fit.jl):
# data assembly, the transmission network g_net, parameter mapping and fold settings. Expects ROOT to be defined and
# OysterFarming, DataFrames, Dates and Statistics to be loaded.

const MD_START, MD_STOP = Date(1996, 10, 1), Date(2024, 11, 30)
const GB_START, GB_STOP = Date(1995, 1, 1), Date(2011, 12, 31)
const FIRST_MD_YEAR = 1999
const DENSITY0 = 20.0
const BACKGROUND_D = 0.0302 / 365
const ARRIVAL = Dict("full" => 2003.15, "nh_adams_point" => 2004.1, "nh_nannie_island" => 2003.2,
                     "nh_oyster_river" => 2002.2, "nh_squamscott_river" => 2002.5, "nh_woodman_point" => 2003.2,
                     "md_upper" => 2003.1, "md_middle" => 2003.2, "md_lower" => 2002.9)

lab = load_params(joinpath(ROOT, "params", "dermo_establishment_lab.toml"))
pw = load_params(joinpath(ROOT, "params", "dermo_powell1996.toml"))
forcing = read_contract(joinpath(ROOT, "data", "processed", "forcing_daily.arrow"), "forcing_daily")
disease = read_contract(joinpath(ROOT, "data", "processed", "disease_obs.arrow"), "disease_obs")
sites_tbl = read_contract(joinpath(ROOT, "data", "processed", "sites.arrow"), "sites")
lat = Dict(r.site_id => r.lat for r in eachrow(sites_tbl))
dermo = disease[disease.pathogen .== "P_marinus", :]

function held(site, a, b)
    f = forcing_for_site(forcing, site)
    have = Dict(r.date => (Float64(r.temp_c), Float64(r.salinity_psu)) for r in eachrow(f))
    dates = collect(a:Day(1):b)
    last = have[minimum(keys(have))]
    T = Float64[]; S = Float64[]
    for d in dates
        haskey(have, d) && (last = have[d])
        push!(T, last[1]); push!(S, last[2])
    end
    return dates, T, S
end

struct SiteData
    id::String
    region::String
    T::Vector{Float64}
    S::Vector{Float64}
    rec::Vector{Int}
    obs::Vector{Float64}
    prev_obs::Vector{Union{Missing,Float64}}
end

function build(site, region, a, b, keep)
    _, T, S = held(site, a, b)
    o = sort(dermo[(dermo.site_id .== site) .& keep.(dermo.date), :], :date)
    o = o[(o.date .>= a) .& (o.date .<= b), :]
    nrow(o) == 0 && return nothing
    rec = [Dates.value(d - a) + 1 for d in o.date]
    prevs = Union{Missing,Float64}[i == 1 ? missing : Float64(o.prevalence[i-1]) for i in 1:nrow(o)]
    return SiteData(site, region, T, S, rec, Float64.(o.prevalence), prevs)
end

mdregion(s) = lat[s] >= 38.8 ? "md_upper" : lat[s] >= 38.45 ? "md_middle" : "md_lower"
data = let fs = Set(forcing.site_id), out = SiteData[]
    for s in sort(unique(dermo.site_id))
        s in fs || continue
        x = startswith(s, "mddnr_") ? build(s, mdregion(s), MD_START, MD_STOP, d -> year(d) >= FIRST_MD_YEAR) :
            startswith(s, "nh_") ? build(s, "gb", GB_START, GB_STOP, d -> true) : nothing
        x === nothing || push!(out, x)
    end
    out
end

function initial(site::SiteData)
    N0 = zeros(EST_N_LEVELS)
    site.region == "gb" ? (N0[1] = DENSITY0) : (N0[1] = DENSITY0 / 2; N0[15] = DENSITY0 / 2)
    return N0
end

const NH = 6
const NW = NH * 2 + NH + NH * NH + NH + NH + 1
softplus(x) = log1p(exp(-abs(x))) + max(x, zero(x))
function g_net(w, T, S)
    x1 = S / 30; x2 = (T - 15) / 10
    o = 0
    W1 = reshape(view(w, o+1:o+2NH), NH, 2); o += 2NH
    b1 = view(w, o+1:o+NH); o += NH
    W2 = reshape(view(w, o+1:o+NH*NH), NH, NH); o += NH * NH
    b2 = view(w, o+1:o+NH); o += NH
    W3 = view(w, o+1:o+NH); o += NH
    b3 = w[o+1]
    h1 = tanh.(W1[:, 1] .* x1 .+ W1[:, 2] .* x2 .+ b1)
    h2 = tanh.(W2 * h1 .+ b2)
    return softplus(sum(W3 .* h2) + b3)
end
nn_weights(θ) = view(θ, 6:length(θ))
logistic(x) = 1 / (1 + exp(-x))
function params_from(θ, region)
    eps = region == "gb" ? exp(θ[5]) : exp(θ[4])
    return EstablishmentParams(lab, pw; beta=exp(θ[1]), mu=exp(θ[2]), epsilon=eps, host=logistic(θ[3]))
end
intro_day(arrival, site) = site.region == "gb" ? clamp(round(Int, (arrival - 1995) * 365.25) + 1, 1, length(site.T)) : 1
rmse(a, b) = sqrt(mean((a .- b) .^ 2))
