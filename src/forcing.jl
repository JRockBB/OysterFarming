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
