# Data assembly shared by the P4.2g/P4.2h scripts: sites from every system (Maryland DNR bars, Great Bay NH beds,
# Ford & Smolowitz 2007 New England sites, Maine qPCR sites) as `Site` records with daily forcing, start state and
# observations. Expects ROOT to be defined and OysterFarming, DataFrames, Dates, Statistics to be loaded.
# Design notes are in scripts/p4_2g_nested.jl.

const LIT = joinpath(ROOT, "data", "literature", "dermo")
const DENSITY0 = 20.0
const BACKGROUND_D = 0.0302 / 365
const SYSTEMS = ["md", "gb", "ne", "me"]
lab = load_params(joinpath(ROOT, "params", "dermo_establishment_lab.toml"))
pw = load_params(joinpath(ROOT, "params", "dermo_powell1996.toml"))
forcing = read_contract(joinpath(ROOT, "data", "processed", "forcing_daily.arrow"), "forcing_daily")
disease = read_contract(joinpath(ROOT, "data", "processed", "disease_obs.arrow"), "disease_obs")
sites_tbl = read_contract(joinpath(ROOT, "data", "processed", "sites.arrow"), "sites")
dermo = disease[disease.pathogen .== "P_marinus", :]

struct Site
    id::String
    system::String
    group::String
    start::Date                     # date of day index 1
    T::Vector{Float64}
    S::Vector{Float64}
    N0::Vector{Float64}
    prev_rec::Vector{Int}
    prev_obs::Vector{Float64}
    prev_n::Vector{Float64}        # oysters per survey (binomial likelihoods, P4.2i)
    pcr::Bool                       # true: compare with true prevalence (qPCR); false: RFTM-detectable prevalence
    wp_rec::Vector{Int}
    wp_obs::Vector{Float64}
    prev_persist::Vector{Union{Missing,Float64}}
end

"Daily T, S over [a, b] from forcing_daily, repeating the last value on days without a row."
function held(site, a, b)
    f = forcing_for_site(forcing, site)
    have = Dict(r.date => (Float64(r.temp_c), Float64(r.salinity_psu)) for r in eachrow(f))
    last = have[minimum(keys(have))]
    T = Float64[]; S = Float64[]
    for d in a:Day(1):b
        haskey(have, d) && (last = have[d])
        push!(T, last[1]); push!(S, last[2])
    end
    return T, S
end
dayidx(d, a) = Dates.value(d - a) + 1
persist(v) = Union{Missing,Float64}[i == 1 ? missing : v[i-1] for i in eachindex(v)]
function start_state(p0, level)
    N0 = zeros(EST_N_LEVELS); N0[1] = DENSITY0 * (1 - p0); N0[level] += DENSITY0 * p0
    return N0
end

sites = Site[]
fsites = Set(forcing.site_id)
src_of = Dict(r.site_id => r.source for r in eachrow(forcing) if startswith(r.site_id, "mddnr_"))
notes_of = Dict(r.site_id => coalesce(r.notes, "") for r in eachrow(sites_tbl))

# Maryland
let a = Date(1996, 10, 1), b = Date(2024, 11, 30)
    for s in sort(unique(dermo.site_id))
        (startswith(s, "mddnr_") && s in fsites) || continue
        o = sort(dermo[(dermo.site_id .== s) .& (year.(dermo.date) .>= 1999), :], :date)
        nrow(o) < 3 && continue
        T, S = held(s, a, b)
        m = match(r"cbp:([A-Z0-9.]+)", src_of[s])
        grp = "md_" * (m === nothing ? s : m[1])
        push!(sites, Site(s, "md", grp, a, T, S, start_state(0.5, 15), dayidx.(o.date, a), Float64.(o.prevalence),
                          fill(30.0, nrow(o)), false, Int[], Float64[], persist(Float64.(o.prevalence))))   # 30 oysters per bar: MD DNR 2024 Fall Survey report p. 14
    end
end
# Great Bay
let a = Date(1995, 1, 1), b = Date(2011, 12, 31)
    for s in sort(unique(dermo.site_id))
        (startswith(s, "nh_") && s in fsites) || continue
        o = sort(dermo[(dermo.site_id .== s) .& (dermo.date .<= b), :], :date)
        T, S = held(s, a, b)
        m = match(r"sonde (grb[a-z]+)", notes_of[s])
        push!(sites, Site(s, "gb", "gb_" * (m === nothing ? s : m[1]), a, T, S, start_state(0.0, 15), dayidx.(o.date, a),
                          Float64.(o.prevalence), Float64.(o.n), false, dayidx.(o.date, a), Float64.(o.intensity_mean),
                          persist(Float64.(o.prevalence))))
    end
end
# New England (Ford & Smolowitz 2007)
"Split one CSV line, honouring double-quoted fields that contain commas."
function csvsplit(l)
    out = String[]; buf = IOBuffer(); q = false
    for c in l
        if c == '"'
            q = !q
        elseif c == ',' && !q
            push!(out, String(take!(buf)))
        else
            write(buf, c)
        end
    end
    push!(out, String(take!(buf)))
    return out
end
function read_csv_simple(f)
    lines = [l for l in readlines(joinpath(LIT, f)) if !startswith(l, "#") && !isempty(strip(l))]
    hdr = csvsplit(lines[1]); rows = [csvsplit(l) for l in lines[2:end]]
    return DataFrame([Symbol(h) => [String(r[i]) for r in rows] for (i, h) in enumerate(hdr)])
end
fs = read_csv_simple("ford2007_smolowitz_fig2_fig3_digitized.csv")
tr = read_csv_simple("ford2007_smolowitz_fig2_temperature_trace.csv")
const NE_SAL = Dict("Wellfleet MA" => 31.0, "Cotuit MA" => 26.0, "Oyster Bay NY" => 24.0, "Great Sound NJ" => 29.5,
                    "Delaware Bay NJ" => 20.0, "Martha's Vineyard MA (Butler's Neck)" => 15.0,
                    "Westport CT" => 27.5, "Milford CT" => 27.5)
for site in sort(unique(fs.site))
    haskey(NE_SAL, site) || continue
    t = tr[tr.site .== site, :]
    nrow(t) == 0 && continue
    tdates = Date.(t.date); tvals = parse.(Float64, t.temp_C)
    pv = fs[(fs.site .== site) .& (fs.variable .== "prevalence_pct"), :]
    wv = fs[(fs.site .== site) .& (fs.variable .== "weighted_prevalence_0to5"), :]
    pdates = Date.(pv.date); pvals = parse.(Float64, pv.value) ./ 100
    o = sortperm(pdates); pdates = pdates[o]; pvals = clamp.(pvals[o], 0.0, 1.0)
    a = pdates[1]; b = max(maximum(pdates), maximum(Date.(wv.date)))
    days = collect(a:Day(1):b)
    # daily temperature: mean of trace vertices per date, linear interpolation between dates, held at the ends
    ud = sort(unique(tdates)); uv = [mean(tvals[tdates .== d]) for d in ud]
    T = [d <= ud[1] ? uv[1] : d >= ud[end] ? uv[end] :
         (i = searchsortedlast(ud, d); uv[i] + (uv[i+1] - uv[i]) * Dates.value(d - ud[i]) / Dates.value(ud[i+1] - ud[i])) for d in days]
    S = fill(NE_SAL[site], length(days))
    wd = Date.(wv.date); wvals = parse.(Float64, wv.value); keepw = wd .> a
    grp = site in ("Westport CT", "Milford CT") ? "ne_LIS" : "ne_" * site
    sid = "ne_" * lowercase(replace(site, r"[^A-Za-z]+" => "_"))
    push!(sites, Site(sid, "ne", grp, a, T, S, start_state(pvals[1], 15), dayidx.(pdates[2:end], a), pvals[2:end], fill(25.0, length(pvals) - 1), false,   # n = 25 per sample (Ford & Smolowitz 2007)
                      dayidx.(wd[keepw], a), wvals[keepw], persist(pvals)[2:end]))
end
# Maine (Marquis et al. 2020), if digitized. Spin-up from 2015-01-01 uninfected with the system's import (a start from
# the first monthly sample, 0/24 or 24/24, pinned the model at 0 or 1); every sample is scored, against model true
# prevalence (qPCR detects light infections; a qPCR detection threshold is P4.2h step 3).
const MAINE_FILE = "marquis2020_maine_qpcr_monthly.csv"
const MAINE_START = Date(2015, 1, 1)
if isfile(joinpath(LIT, MAINE_FILE))
    me = read_csv_simple(MAINE_FILE)
    for (site, model_site) in (("Jack's Point", "dam_upper"), ("Prentiss Island", "dam_lower"))
        r = me[occursin.(lowercase(site), lowercase.(me.site)), :]
        nrow(r) == 0 && continue
        d = Date.(r.date); v = parse.(Float64, r.prevalence); o = sortperm(d); d = d[o]; v = v[o]
        v = any(v .> 1) ? v ./ 100 : v
        T, S = held(model_site, MAINE_START, d[end])
        sid = "me_" * lowercase(replace(site, r"[^A-Za-z]+" => "_"))
        push!(sites, Site(sid, "me", "me_" * model_site, MAINE_START, T, S, start_state(0.0, 1), dayidx.(d, MAINE_START), v, fill(24.0, length(v)), true,   # n = 24 per sample (Marquis et al. 2020)
                          Int[], Float64[], persist(v)))
    end
end
