# Helpers for selecting model forcing from a forcing_daily contract table.

"""
    forcing_for_site(forcing, site_id) -> DataFrame

Rows of a `forcing_daily` table for one site, sorted by date.
"""
function forcing_for_site(forcing, site_id::AbstractString)
    df = DataFrame(forcing; copycols=false)
    sub = df[df.site_id .== site_id, :]
    nrow(sub) > 0 || throw(ArgumentError("no forcing rows for site \"$site_id\""))
    return sort(sub, :date)
end

"""
    continuous_blocks(dates) -> Vector{UnitRange{Int}}

Index ranges of runs of consecutive days in a sorted date vector. Use it to
pick spans that [`simulate_individual`](@ref) can run on, since forcing may
still contain unfilled gaps (for example winters without observations).
"""
function continuous_blocks(dates::AbstractVector{Date})
    isempty(dates) && return UnitRange{Int}[]
    issorted(dates) || throw(ArgumentError("dates must be sorted"))
    blocks = UnitRange{Int}[]
    start = 1
    for i in 2:length(dates)
        if dates[i] - dates[i-1] != Day(1)
            push!(blocks, start:i-1)
            start = i
        end
    end
    push!(blocks, start:length(dates))
    return blocks
end

"""
    longest_block(df) -> DataFrame

The longest run of consecutive days in a single-site forcing table.
"""
function longest_block(df)
    blocks = continuous_blocks(df.date)
    return df[blocks[argmax(length.(blocks))], :]
end

"""
    time_to_market(p, forcing, seed_date, u0; market_shell_height_mm=76.0) -> Union{Int,Nothing}

Days from `seed_date` until simulated shell height first reaches
`market_shell_height_mm` (default 76 mm, i.e. 3 inches, the market size named
in PLAN.md P2.7). The run starts at `seed_date`, which must fall inside the
forcing, and continues through consecutive forcing days. Returns `nothing` if
market size is not reached before the forcing ends or hits a gap.
"""
function time_to_market(p::DEBParams, forcing, seed_date::Date, u0::NamedTuple;
                        market_shell_height_mm::Real=76.0)
    df = sort(DataFrame(forcing; copycols=false), :date)
    start = findfirst(==(seed_date), df.date)
    start === nothing && throw(ArgumentError("seed_date $seed_date is not in the forcing"))
    block = continuous_blocks(df.date[start:end])[1]
    run = df[start .+ block .- 1, :]
    out = simulate_individual(p, run, u0)
    i = findfirst(>=(market_shell_height_mm), out.shell_height_mm)
    return i === nothing ? nothing : i
end
