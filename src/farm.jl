# Farm layer (PLAN.md P3.1–P3.3): cohorts on one lease, background mortality,
# food depletion in a box around the lease, and harvest at market size.
#
# Each cohort is represented by one DEB individual (all oysters in a cohort
# share the same state) and a deterministic count of survivors.

"""
    FarmParams

Lease-scale parameters, read from a parameter file with
`FarmParams(load_params(path))` (every value needs a source):

- `mortality_rate_d`    background mortality, d⁻¹ (daily survival exp(−m))
- `flushing_time_d`     flushing time τ of the water box around the lease, d
- `box_volume_L`        volume of that water box, L
- `clearance_max_L_cm2_d` maximum surface-specific clearance rate {C_R},
                        L d⁻¹ cm⁻² (per V^(2/3)); converts feeding into chlorophyll removal
"""
struct FarmParams
    mortality_rate_d::Float64
    flushing_time_d::Float64
    box_volume_L::Float64
    clearance_max_L_cm2_d::Float64
end

function FarmParams(params::AbstractDict{String,<:Param})
    pv(n) = param_value(params, n)
    fp = FarmParams(pv("mortality_rate_d"), pv("flushing_time_d"), pv("box_volume_L"),
                    pv("clearance_max_L_cm2_d"))
    fp.mortality_rate_d >= 0 || throw(ArgumentError("mortality_rate_d must be non-negative"))
    fp.flushing_time_d > 0 || throw(ArgumentError("flushing_time_d must be positive"))
    fp.box_volume_L > 0 || throw(ArgumentError("box_volume_L must be positive"))
    fp.clearance_max_L_cm2_d > 0 || throw(ArgumentError("clearance_max_L_cm2_d must be positive"))
    return fp
end

"""
    Cohort(id, seed_date, n_seeded, seed_shell_height_mm; seed_scaled_reserve=1.0)

A group of oysters seeded together. `seed_scaled_reserve` is the scaled
reserve density e (0–1) at seeding, a scenario input.
"""
struct Cohort
    id::String
    seed_date::Date
    n_seeded::Float64
    seed_shell_height_mm::Float64
    seed_scaled_reserve::Float64
end

Cohort(id, seed_date, n_seeded, seed_shell_height_mm; seed_scaled_reserve=1.0) =
    Cohort(String(id), seed_date, Float64(n_seeded), Float64(seed_shell_height_mm), Float64(seed_scaled_reserve))

"""
    HarvestRule(; market_shell_height_mm=76.0, weekly_capacity=Inf)

Harvest a cohort once its shell height reaches `market_shell_height_mm`,
taking at most `weekly_capacity` oysters per calendar week (Monday start)
across all cohorts, oldest cohort first (P3.3). Management choices, not
biological parameters.
"""
struct HarvestRule
    market_shell_height_mm::Float64
    weekly_capacity::Float64
end

function HarvestRule(; market_shell_height_mm::Real=76.0, weekly_capacity::Real=Inf)
    market_shell_height_mm > 0 || throw(ArgumentError("market_shell_height_mm must be positive"))
    weekly_capacity >= 0 || throw(ArgumentError("weekly_capacity must be non-negative"))
    return HarvestRule(Float64(market_shell_height_mm), Float64(weekly_capacity))
end

"""
    max_ingestion_ug_d(V, c_T, c_S, p, fp) -> Float64

Maximum chlorophyll ingestion of one oyster (f = 1), µg chl-a d⁻¹. With a
Holling type II response, clearance is CR = {C_R}·(1 − f)·c_T·c_S·V^(2/3) and
ingestion CR·X equals {C_R}·X_K·f·c_T·c_S·V^(2/3) (Kooijman 2010, section 2.3;
{C_R} from Lavaud et al. 2017, Table 2), so

    J_X,max = {C_R} · X_K · c_T · c_S · V^(2/3)

with X_K in µg chl-a L⁻¹ and {C_R} in L d⁻¹ cm⁻².
"""
max_ingestion_ug_d(V, c_T, c_S, p::DEBParams, fp::FarmParams) =
    c_S * c_T * fp.clearance_max_L_cm2_d * p.X_K * V^(2 / 3)

"""
    depleted_food(X_in, total_max_ingestion_ug_d, X_K, fp) -> Float64

Food density X (µg chl-a L⁻¹) around the lease from the box model of PLAN.md P3.2,

    dX/dt = (X_in − X)/τ − (Σ_i n_i J_X,max,i) · X/(X + X_K) / V_box,

taken at its steady state for the day (quasi-steady assumption: the box
flushes much faster than oysters grow). Solved by bisection on [0, X_in],
where the right-hand side changes sign monotonically.
"""
function depleted_food(X_in::Real, total_max_ingestion_ug_d::Real, X_K::Real, fp::FarmParams)
    X_in >= 0 || throw(DomainError(X_in, "inflow food must be non-negative"))
    total_max_ingestion_ug_d >= 0 || throw(DomainError(total_max_ingestion_ug_d, "ingestion must be non-negative"))
    (X_in == 0 || total_max_ingestion_ug_d == 0) && return float(X_in)
    g(X) = (X_in - X) / fp.flushing_time_d - total_max_ingestion_ug_d * X / (X + X_K) / fp.box_volume_L
    lo, hi = 0.0, float(X_in)
    for _ in 1:100
        mid = (lo + hi) / 2
        g(mid) > 0 ? (lo = mid) : (hi = mid)
        hi - lo < 1e-12 * max(1.0, X_in) && break
    end
    return (lo + hi) / 2
end

"""
    simulate_farm(p, fp, forcing, cohorts, rule; depletion=true, closed_dates=Set{Date}(),
                  extra_mortality_d=d -> 0.0) -> (daily, harvests)

Run cohorts on one lease over consecutive forcing days (one site). Each day:
(1) food around the lease from [`depleted_food`](@ref) using all live oysters
(or the inflow food if `depletion=false`); (2) every seeded cohort advances one
DEB day; (3) mortality removes a fraction 1 − exp(−(m + m_x)), where m is the
background rate and m_x = `extra_mortality_d(date)` an optional additional daily
rate (e.g. disease, PLAN.md P4.6); (4) cohorts at market size are harvested
within the weekly capacity, except on `closed_dates` (harvest closures). The
capacity in force on a day is `weekly_capacity_fn(date)` (default: the rule's
fixed capacity), e.g. reduced during the Vibrio control period (P4.5).

Returns `daily` (one row per cohort per day after seeding: `date`, `cohort_id`,
`n_alive`, `n_dead`, `n_harvested`, `shell_height_mm`, `dry_tissue_g`,
`food_in_ug_l`, `food_ug_l`, `harvest_blocked` = market-size oysters present on a
closed day) and `harvests` (`date`, `cohort_id`, `count`,
`shell_height_mm`, `dry_tissue_g`). Counts are expected numbers (real-valued).
"""
function simulate_farm(p::DEBParams, fp::FarmParams, forcing, cohorts::AbstractVector{Cohort},
                       rule::HarvestRule; depletion::Bool=true, food_col::Symbol=:chl_ug_l,
                       closed_dates::AbstractSet{Date}=Set{Date}(), extra_mortality_d::Function=d -> 0.0,
                       weekly_capacity_fn::Function=d -> rule.weekly_capacity)
    df = DataFrame(forcing; copycols=false)
    n = nrow(df)
    n > 0 || throw(ArgumentError("forcing is empty"))
    all(diff(df.date) .== Day(1)) || throw(ArgumentError("forcing dates must be consecutive days"))
    any(ismissing, df[!, food_col]) && throw(ArgumentError("food column has missing values"))
    length(unique(c.id for c in cohorts)) == length(cohorts) || throw(ArgumentError("cohort ids must be unique"))
    for c in cohorts
        df.date[1] <= c.seed_date <= df.date[end] ||
            throw(ArgumentError("cohort $(c.id) seed date $(c.seed_date) is outside the forcing"))
    end
    order = sortperm(cohorts; by=c -> (c.seed_date, c.id))   # oldest first for harvest priority
    cs = cohorts[order]
    k = length(cs)
    steppers = Vector{Union{Nothing,IndividualStepper}}(nothing, k)
    alive = zeros(k); dead = zeros(k); harvested = zeros(k)
    V = zeros(k); E = zeros(k); E_R = zeros(k)

    daily = DataFrame(date=Date[], cohort_id=String[], n_alive=Float64[], n_dead=Float64[],
                      n_harvested=Float64[], shell_height_mm=Float64[], dry_tissue_g=Float64[],
                      food_in_ug_l=Float64[], food_ug_l=Float64[], harvest_blocked=Bool[])
    harvests = DataFrame(date=Date[], cohort_id=String[], count=Float64[], shell_height_mm=Float64[],
                         dry_tissue_g=Float64[])
    week = Date(0)
    used_this_week = 0.0
    for i in 1:n
        date, temp_c, sal = df.date[i], df.temp_c[i], df.salinity_psu[i]
        X_in = Float64(df[i, food_col])
        for j in 1:k
            if steppers[j] === nothing && cs[j].seed_date == date
                u0 = initial_state(p; shell_height_mm=cs[j].seed_shell_height_mm, e=cs[j].seed_scaled_reserve)
                steppers[j] = IndividualStepper(p, u0)
                alive[j] = cs[j].n_seeded
                E[j], V[j], E_R[j] = u0.E, u0.V, u0.E_R
            end
        end
        c_T = arrhenius_factor(temp_c, p)
        c_S = salinity_factor(sal, p)
        X = X_in
        if depletion
            demand = sum((alive[j] * max_ingestion_ug_d(V[j], c_T, c_S, p, fp) for j in 1:k
                          if steppers[j] !== nothing); init=0.0)
            X = depleted_food(X_in, demand, p.X_K, fp)
        end
        f = functional_response(X, p.X_K)
        m_x = extra_mortality_d(date)
        m_x >= 0 || throw(DomainError(m_x, "extra mortality rate must be non-negative"))
        survival = exp(-(fp.mortality_rate_d + m_x))
        closed = date in closed_dates
        if firstdayofweek(date) != week
            week = firstdayofweek(date)
            used_this_week = 0.0
        end
        for j in 1:k
            s = steppers[j]
            s === nothing && continue
            r = advance_day!(s, temp_c, sal, f)
            E[j], V[j], E_R[j] = r.E, r.V, r.E_R
            died = alive[j] * (1 - survival)
            alive[j] -= died
            dead[j] += died
            sh = shell_height_mm(V[j], p.delta_M)
            dw = dry_tissue_g(E[j], V[j], E_R[j], p)
            blocked = false
            if sh >= rule.market_shell_height_mm && alive[j] > 0 && closed
                blocked = true
            elseif sh >= rule.market_shell_height_mm && alive[j] > 0
                take = min(alive[j], weekly_capacity_fn(date) - used_this_week)
                if take > 0
                    alive[j] -= take
                    harvested[j] += take
                    used_this_week += take
                    push!(harvests, (date, cs[j].id, take, sh, dw))
                end
            end
            push!(daily, (date, cs[j].id, alive[j], dead[j], harvested[j], sh, dw, X_in, X, blocked))
        end
    end
    return daily, harvests
end
