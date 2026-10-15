# Decision layer (PLAN.md P5.2, P5.3): evaluate farm management choices by NPV under
# parameter uncertainty, and search a grid of choices.

"""
    FarmDecision(seed_month, seeds_per_cohort, market_shell_height_mm)

Management choices: month in which a cohort is seeded (on the 15th), number of
seeds per cohort, and the shell height at which oysters are harvested.
"""
struct FarmDecision
    seed_month::Int
    seeds_per_cohort::Float64
    market_shell_height_mm::Float64
end

"""
    DecisionSetup

Everything except the decision: a single-site forcing block, weighted DEB
parameter scenarios (e.g. X_K posterior quantiles), farm, economic and nitrogen
parameters, harvest capacity, Vibrio handling (capacity factor and cost per
oyster during the control period; scenario inputs), seed size, number of yearly
cohorts, and the risk level for the lower-quantile objective. `extra_mortality_d`
(date -> daily rate, default 0) adds mortality to every cohort through
[`simulate_farm`](@ref), e.g. Dermo deaths from a disease scenario (PLAN.md P4.2f);
the 12-argument constructor leaves it at 0.
"""
struct DecisionSetup
    forcing::DataFrame
    scenarios::Vector{Tuple{Float64,DEBParams{Float64}}}
    farm::FarmParams
    econ::EconParams
    nitrogen::NitrogenParams
    vibrio::VibrioRule
    weekly_capacity::Float64
    vibrio_capacity_factor::Float64
    vibrio_cost_per_oyster::Float64
    seed_shell_height_mm::Float64
    n_cohorts::Int
    risk_quantile::Float64
    extra_mortality_d::Function
end

DecisionSetup(forcing, scenarios, farm, econ, nitrogen, vibrio, weekly_capacity, vibrio_capacity_factor,
              vibrio_cost_per_oyster, seed_shell_height_mm, n_cohorts, risk_quantile) =
    DecisionSetup(forcing, scenarios, farm, econ, nitrogen, vibrio, weekly_capacity, vibrio_capacity_factor,
                  vibrio_cost_per_oyster, seed_shell_height_mm, n_cohorts, risk_quantile, d -> 0.0)

"Weighted empirical quantile: smallest x with cumulative weight ≥ q."
function weighted_quantile(x::AbstractVector, w::AbstractVector, q::Real)
    o = sortperm(x)
    c = cumsum(w[o]) ./ sum(w)
    return x[o][findfirst(>=(q - 1e-12), c)]
end

"""
    evaluate_decision(d::FarmDecision, s::DecisionSetup) -> NamedTuple

Run the farm under each parameter scenario and return `expected_npv`
(probability-weighted mean), `lower_npv` (weighted `risk_quantile` quantile),
`expected_nitrogen_kg`, `expected_harvested`, and per-scenario `npvs`.

Cohorts are seeded on the 15th of `seed_month` in the first `n_cohorts` years in
which that date falls inside the forcing block. Each starts at
`seed_shell_height_mm` with reserve at equilibrium with the mean chlorophyll of its
first 14 days. The NPV horizon is the calendar years of the forcing block; oysters
still unharvested at the end are given no value.
"""
function evaluate_decision(d::FarmDecision, s::DecisionSetup)
    1 <= d.seed_month <= 12 || throw(ArgumentError("seed_month must be 1-12"))
    d.seeds_per_cohort > 0 || throw(ArgumentError("seeds_per_cohort must be positive"))
    f = s.forcing
    first_date, last_date = first(f.date), last(f.date)
    seed_dates = Date[]
    for y in year(first_date):year(last_date)
        sd = Date(y, d.seed_month, 15)
        (first_date <= sd && sd + Day(14) <= last_date) && push!(seed_dates, sd)
        length(seed_dates) == s.n_cohorts && break
    end
    isempty(seed_dates) && throw(ArgumentError("no seeding date fits in the forcing block"))
    rule = HarvestRule(market_shell_height_mm=d.market_shell_height_mm, weekly_capacity=s.weekly_capacity)
    capfn = vibrio_capacity(s.vibrio, s.weekly_capacity, s.vibrio_capacity_factor)
    npvs = Float64[]; nitro = Float64[]; harvested = Float64[]; weights = Float64[]
    for (w, p) in s.scenarios
        cohorts = map(enumerate(seed_dates)) do (i, sd)
            food0 = mean(f.chl_ug_l[(f.date .>= sd) .& (f.date .< sd + Day(14))])
            Cohort("c$i", sd, d.seeds_per_cohort, s.seed_shell_height_mm;
                   seed_scaled_reserve=functional_response(food0, p.X_K))
        end
        _, h = simulate_farm(p, s.farm, f, cohorts, rule; weekly_capacity_fn=capfn,
                             extra_mortality_d=s.extra_mortality_d)
        vib = Dict{Int,Float64}()
        for r in eachrow(h)
            in_vibrio_control(r.date, s.vibrio) &&
                (vib[year(r.date)] = get(vib, year(r.date), 0.0) + r.count * s.vibrio_cost_per_oyster)
        end
        cf = farm_cashflows(cohorts, h, s.econ, year(first_date), year(last_date); extra_cost_per_year=vib)
        push!(npvs, npv(cf, s.econ.discount_rate))
        push!(nitro, nitrogen_removed_g(h, s.nitrogen) / 1000)
        push!(harvested, sum(h.count; init=0.0))
        push!(weights, w)
    end
    wn = weights ./ sum(weights)
    return (expected_npv=sum(wn .* npvs), lower_npv=weighted_quantile(npvs, wn, s.risk_quantile),
            expected_nitrogen_kg=sum(wn .* nitro), expected_harvested=sum(wn .* harvested), npvs=npvs)
end

"""
    grid_search(s::DecisionSetup; seed_months, seeds_per_cohort, market_sizes_mm) -> DataFrame

Evaluate every combination (PLAN.md P5.3) and return one row per decision,
sorted by expected NPV (best first).
"""
function grid_search(s::DecisionSetup; seed_months, seeds_per_cohort, market_sizes_mm)
    rows = DataFrame(seed_month=Int[], seeds_per_cohort=Float64[], market_shell_height_mm=Float64[],
                     expected_npv=Float64[], lower_npv=Float64[], expected_nitrogen_kg=Float64[],
                     expected_harvested=Float64[])
    for m in seed_months, n in seeds_per_cohort, ms in market_sizes_mm
        r = evaluate_decision(FarmDecision(m, n, ms), s)
        push!(rows, (m, n, ms, r.expected_npv, r.lower_npv, r.expected_nitrogen_kg, r.expected_harvested))
    end
    return sort(rows, :expected_npv; rev=true)
end
