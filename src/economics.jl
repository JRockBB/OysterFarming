# Reference-farm economics and net present value (PLAN.md P5.1).

"""
    EconParams

Prices and costs for a farm, read with
`EconParams(load_params("params/economics_reference_farm.toml"))`.
"""
struct EconParams
    price_per_oyster::Float64
    seed_cost_per_seed::Float64
    variable_cost_per_oyster::Float64
    fixed_cost_per_year::Float64
    capital_investment::Float64
    discount_rate::Float64
end

function EconParams(params::AbstractDict{String,<:Param})
    pv(n) = param_value(params, n)
    e = EconParams(pv("price_per_oyster"), pv("seed_cost_per_seed"), pv("variable_cost_per_oyster"),
                   pv("fixed_cost_per_year"), pv("capital_investment"), pv("discount_rate"))
    e.discount_rate > -1 || throw(ArgumentError("discount rate must exceed -1"))
    return e
end

"""
    farm_cashflows(cohorts, harvests, econ, first_year, last_year; extra_cost_per_year=Dict()) -> DataFrame

Annual cash flows (USD) for years `first_year:last_year`, with year index
t = year − first_year:

    CF_t = Σ harvest count × (price − variable cost per oyster)        (harvests in that year)
         − Σ n_seeded × seed cost                                       (cohorts seeded that year)
         − fixed cost per year
         − capital investment                                           (t = 0 only)
         − extra_cost_per_year[year]                                    (e.g. Vibrio handling, P4.5)

Prices and costs follow Engle & van Senten (Table 11) categories; see
params/economics_reference_farm.toml for sources.
"""
function farm_cashflows(cohorts::AbstractVector{Cohort}, harvests::AbstractDataFrame, econ::EconParams,
                        first_year::Integer, last_year::Integer; extra_cost_per_year::AbstractDict=Dict{Int,Float64}())
    last_year >= first_year || throw(ArgumentError("last_year must not precede first_year"))
    years = collect(first_year:last_year)
    revenue = zeros(length(years)); seed = zeros(length(years)); other = zeros(length(years))
    for r in eachrow(harvests)
        i = year(r.date) - first_year + 1
        1 <= i <= length(years) || throw(ArgumentError("harvest on $(r.date) outside the horizon"))
        revenue[i] += r.count * (econ.price_per_oyster - econ.variable_cost_per_oyster)
    end
    for c in cohorts
        i = year(c.seed_date) - first_year + 1
        1 <= i <= length(years) || throw(ArgumentError("cohort $(c.id) seeded outside the horizon"))
        seed[i] += c.n_seeded * econ.seed_cost_per_seed
    end
    for (i, y) in enumerate(years)
        other[i] = econ.fixed_cost_per_year + get(extra_cost_per_year, y, 0.0) + (i == 1 ? econ.capital_investment : 0.0)
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
