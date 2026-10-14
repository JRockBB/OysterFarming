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

"""
    cycle_years(df, source_years, n_years; start_year=first(source_years)) -> DataFrame

Build an `n_years`-long daily forcing series by repeating the complete calendar years
`source_years` of a single-site table in order (year 1 = source_years[1], ...), re-dated to
consecutive calendar years from `start_year`. Used to run farms longer than the measured
record (PLAN.md P5.3); it treats the measured years as a representative weather sequence.
Leap days: a 366-day source year placed in a non-leap year drops 29 February; a 365-day source
year placed in a leap year repeats 28 February. Every source year must be complete.
"""
function cycle_years(df, source_years::AbstractVector{<:Integer}, n_years::Integer;
                     start_year::Integer=first(source_years))
    n_years >= 1 || throw(ArgumentError("n_years must be at least 1"))
    src = DataFrame(df; copycols=false)
    parts = DataFrame[]
    for i in 1:n_years
        sy = source_years[mod1(i, length(source_years))]
        ty = start_year + i - 1
        s = src[year.(src.date) .== sy, :]
        nrow(s) == daysinyear(sy) || throw(ArgumentError("source year $sy is not complete"))
        if isleapyear(sy) && !isleapyear(ty)
            s = s[.!((month.(s.date) .== 2) .& (day.(s.date) .== 29)), :]
        elseif !isleapyear(sy) && isleapyear(ty)
            i28 = findfirst(d -> month(d) == 2 && day(d) == 28, s.date)
            s = vcat(s[1:i28, :], s[i28:i28, :], s[i28+1:end, :])
        end
        s = copy(s)
        s.date = collect(Date(ty, 1, 1):Day(1):Date(ty, 12, 31))
        push!(parts, s)
    end
    return vcat(parts...)
end
