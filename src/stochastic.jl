# Monte Carlo driver (PLAN.md P4.6).

"""
    run_ensemble(rng, n, run) -> DataFrame

Run `run(rng_i)` for `n` independent members. Each member gets its own
`Xoshiro` generator seeded from `rng`, so the same `rng` seed reproduces every
member exactly, and members do not share random streams. `run` must return a
NamedTuple; the result is one row per member with a `member` column added.
"""
function run_ensemble(rng::AbstractRNG, n::Integer, run::Function)
    n > 0 || throw(ArgumentError("n must be positive"))
    seeds = rand(rng, UInt64, n)
    rows = [merge((member=i,), run(Xoshiro(seeds[i]))) for i in 1:n]
    return DataFrame(rows)
end

"""
    farm_outcome(daily, harvests) -> NamedTuple

Summary of one farm run from [`simulate_farm`](@ref): oysters `seeded`,
`harvested`, `dead`, `mortality_fraction` (dead / seeded), `lost_harvest_days`
(days on which market-size oysters were present but the area was closed) and
`harvest_days` (days with any harvest).
"""
function farm_outcome(daily::AbstractDataFrame, harvests::AbstractDataFrame)
    last_rows = combine(groupby(daily, :cohort_id), :n_alive => last => :alive, :n_dead => last => :dead,
                        :n_harvested => last => :harvested)
    seeded = sum(last_rows.alive .+ last_rows.dead .+ last_rows.harvested)
    dead = sum(last_rows.dead)
    return (seeded=seeded, harvested=sum(last_rows.harvested), dead=dead,
            mortality_fraction=seeded > 0 ? dead / seeded : 0.0,
            lost_harvest_days=length(unique(daily.date[daily.harvest_blocked])),
            harvest_days=length(unique(harvests.date)))
end
