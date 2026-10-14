# Reference-farm economics and net present value (PLAN.md P5.1).

"""
    EconParams

Prices and costs for a farm, read with
`EconParams(load_params("params/economics_maine_maa.toml"))` (or another economics file).

Optional parts, used when the parameter file has them:
- size-dependent prices: `price_tiers` holds (minimum shell height mm, USD per oyster) pairs; a
  harvested oyster gets the price of the highest tier it reaches and 0 below the lowest tier
  (files: `petite_min_shell_height_mm`, `price_petite`, `legal_min_shell_height_mm`, `price_legal`);
- capital scaled by farm size: `capital_per_acre` × acres planted, with acres = largest number of
  seeds planted in one calendar year / `seeds_per_acre` (files: `capital_per_acre`, `seeds_per_acre`);
- salvage value: with `capital_life_years` = L > 0, the capital outlay is credited back in the last
  year at straight-line book value, outlay × max(0, 1 − n_years/L), so a horizon shorter than the
  gear's life is not charged the whole investment.
"""
struct EconParams
    price_per_oyster::Float64
    seed_cost_per_seed::Float64
    variable_cost_per_oyster::Float64
    fixed_cost_per_year::Float64
    capital_investment::Float64
    discount_rate::Float64
    price_tiers::Vector{Tuple{Float64,Float64}}
    capital_per_acre::Float64
    seeds_per_acre::Float64
    capital_life_years::Float64
end

EconParams(price, seed, variable, fixed, capital, r) =
    EconParams(price, seed, variable, fixed, capital, r, Tuple{Float64,Float64}[], 0.0, 0.0, 0.0)

function EconParams(params::AbstractDict{String,<:Param})
    pv(n) = param_value(params, n)
    opt(n) = haskey(params, n) ? Float64(pv(n)) : 0.0
    tiers = Tuple{Float64,Float64}[]
    for (lo, pr) in (("petite_min_shell_height_mm", "price_petite"), ("legal_min_shell_height_mm", "price_legal"))
        haskey(params, lo) && push!(tiers, (Float64(pv(lo)), Float64(pv(pr))))
    end
    sort!(tiers)
    e = EconParams(pv("price_per_oyster"), pv("seed_cost_per_seed"), pv("variable_cost_per_oyster"),
                   pv("fixed_cost_per_year"), pv("capital_investment"), pv("discount_rate"),
                   tiers, opt("capital_per_acre"), opt("seeds_per_acre"), opt("capital_life_years"))
    e.discount_rate > -1 || throw(ArgumentError("discount rate must exceed -1"))
    (e.capital_per_acre > 0) == (e.seeds_per_acre > 0) ||
        throw(ArgumentError("capital_per_acre and seeds_per_acre must be given together"))
    return e
end

"""
    price_at(econ, shell_height_mm) -> Float64

Price per oyster at harvest size: the single `price_per_oyster` if no size tiers are set,
otherwise the price of the highest tier whose minimum is reached (0 below the lowest tier).
"""
function price_at(econ::EconParams, shell_height_mm::Real)
    isempty(econ.price_tiers) && return econ.price_per_oyster
    i = findlast(t -> shell_height_mm >= t[1], econ.price_tiers)
    return i === nothing ? 0.0 : econ.price_tiers[i][2]
end

"""
    capital_outlay(econ, cohorts) -> Float64

`capital_investment` plus, when set, `capital_per_acre` × (largest number of seeds planted in
one calendar year / `seeds_per_acre`).
"""
function capital_outlay(econ::EconParams, cohorts::AbstractVector{Cohort})
    econ.capital_per_acre > 0 || return econ.capital_investment
    per_year = Dict{Int,Float64}()
    for c in cohorts
        per_year[year(c.seed_date)] = get(per_year, year(c.seed_date), 0.0) + c.n_seeded
    end
    acres = isempty(per_year) ? 0.0 : maximum(values(per_year)) / econ.seeds_per_acre
    return econ.capital_investment + econ.capital_per_acre * acres
end

"""
    farm_cashflows(cohorts, harvests, econ, first_year, last_year; extra_cost_per_year=Dict()) -> DataFrame

Annual cash flows (USD) for years `first_year:last_year`, with year index
t = year − first_year:

    CF_t = Σ harvest count × (price_at(shell height) − variable cost)  (harvests in that year)
         − Σ n_seeded × seed cost                                       (cohorts seeded that year)
         − fixed cost per year
         − capital_outlay(econ, cohorts)                                (t = 0 only)
         + capital_outlay × max(0, 1 − n_years / capital_life_years)   (last year, if life set)
         − extra_cost_per_year[year]                                    (e.g. Vibrio handling, P4.5)

Size-dependent prices need a `shell_height_mm` column in `harvests` (as returned by
[`simulate_farm`](@ref)). Sources are in the parameter files.
"""
function farm_cashflows(cohorts::AbstractVector{Cohort}, harvests::AbstractDataFrame, econ::EconParams,
                        first_year::Integer, last_year::Integer; extra_cost_per_year::AbstractDict=Dict{Int,Float64}())
    last_year >= first_year || throw(ArgumentError("last_year must not precede first_year"))
    years = collect(first_year:last_year)
    revenue = zeros(length(years)); seed = zeros(length(years)); other = zeros(length(years))
    sized = "shell_height_mm" in names(harvests)
    isempty(econ.price_tiers) || sized ||
        throw(ArgumentError("size-dependent prices need a shell_height_mm column in harvests"))
    for r in eachrow(harvests)
        i = year(r.date) - first_year + 1
        1 <= i <= length(years) || throw(ArgumentError("harvest on $(r.date) outside the horizon"))
        p = sized ? price_at(econ, r.shell_height_mm) : econ.price_per_oyster
        revenue[i] += r.count * (p - econ.variable_cost_per_oyster)
    end
    for c in cohorts
        i = year(c.seed_date) - first_year + 1
        1 <= i <= length(years) || throw(ArgumentError("cohort $(c.id) seeded outside the horizon"))
        seed[i] += c.n_seeded * econ.seed_cost_per_seed
    end
    capital = capital_outlay(econ, cohorts)
    for (i, y) in enumerate(years)
        other[i] = econ.fixed_cost_per_year + get(extra_cost_per_year, y, 0.0) + (i == 1 ? capital : 0.0)
    end
    if econ.capital_life_years > 0
        other[end] -= capital * max(0.0, 1 - length(years) / econ.capital_life_years)
    end
    return DataFrame(year=years, t=years .- first_year, net_harvest_revenue=revenue, seed_cost=seed,
                     other_costs=other, cashflow=revenue .- seed .- other)
end

"""
    npv(cashflows, r) -> Float64

Net present value Σ_t CF_t / (1 + r)^t, with cash flows at the end of year t
counted from t = 0 (the first year is not discounted).
"""
npv(cashflows::AbstractDataFrame, r::Real) = sum(cashflows.cashflow ./ (1 + r) .^ cashflows.t)
