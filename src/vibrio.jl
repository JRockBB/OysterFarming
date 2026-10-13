# Maine Vibrio parahaemolyticus control rule as harvest constraints (PLAN.md P4.5).
#
# Maine DMR Rules Chapter 115, Vibrio Control Plan (as amended March 19, 2024):
# statewide for oysters from June 1 through October 15 inclusive. During that period
# harvested oysters must reach a certified dealer and temperature control within
# 2 hours (option A) or be shaded and iced immediately (option B); over-exposed
# product must be resubmerged for 24 h or 7 d, or destroyed. The rule is calendar
# based (no temperature trigger in the current text).
#
# The model represents the rule as (1) a flag on every harvest inside the control
# period, (2) an optional reduction of weekly harvest capacity during the period,
# and (3) an optional handling cost per oyster harvested during the period. The
# capacity factor and cost are scenario inputs: no published estimate of either
# was found (PLAN.md P4.5, P5.1).

"""
    VibrioRule

Control period of the Maine Vibrio rule, read with
`VibrioRule(load_params("params/vibrio_maine.toml"))`.
"""
struct VibrioRule
    start_month::Int
    start_day::Int
    end_month::Int
    end_day::Int
end

function VibrioRule(params::AbstractDict{String,<:Param})
    pv(n) = Int(param_value(params, n))
    r = VibrioRule(pv("control_start_month"), pv("control_start_day"), pv("control_end_month"), pv("control_end_day"))
    Date(2001, r.start_month, r.start_day) <= Date(2001, r.end_month, r.end_day) ||
        throw(ArgumentError("control period must not wrap the new year"))
    return r
end

"""
    in_vibrio_control(date, r::VibrioRule) -> Bool

True when `date` falls in the control period, both ends inclusive
(Ch. 115 §115.02: "June 1st through October 15th, inclusive").
"""
function in_vibrio_control(date::Date, r::VibrioRule)
    y = year(date)
    return Date(y, r.start_month, r.start_day) <= date <= Date(y, r.end_month, r.end_day)
end

"""
    vibrio_capacity(r, base_weekly_capacity, controlled_factor) -> Function

Weekly harvest capacity as a function of date, for `simulate_farm`'s
`weekly_capacity_fn`: `base × controlled_factor` inside the control period and
`base` outside it. `controlled_factor` (0–1) is a scenario input representing the
extra handling time (icing or 2-hour dealer delivery).
"""
function vibrio_capacity(r::VibrioRule, base_weekly_capacity::Real, controlled_factor::Real)
    0 <= controlled_factor <= 1 || throw(DomainError(controlled_factor, "controlled_factor must be in [0, 1]"))
    return date -> in_vibrio_control(date, r) ? base_weekly_capacity * controlled_factor : float(base_weekly_capacity)
end

"""
    vibrio_handling_cost(harvests, r; cost_per_oyster) -> Float64

Total extra handling cost for harvests made inside the control period:
Σ count × `cost_per_oyster` over harvest rows whose date is controlled.
`cost_per_oyster` (currency per oyster, e.g. ice) is a scenario input.
"""
function vibrio_handling_cost(harvests::AbstractDataFrame, r::VibrioRule; cost_per_oyster::Real)
    cost_per_oyster >= 0 || throw(DomainError(cost_per_oyster, "cost must be non-negative"))
    return sum((row.count * cost_per_oyster for row in eachrow(harvests) if in_vibrio_control(row.date, r)); init=0.0)
end
