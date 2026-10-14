# Data contract readers and validators.
#
# The schemas mirror docs/PLAN.md ("Data contract") and the Python validators in
# data-pipeline/src/oyster_data/contract.py. Change all three together.

"""
    ContractError(table_name, problems)

Thrown when a table does not follow the data contract. `problems` lists every
violation found, not just the first.
"""
struct ContractError <: Exception
    table_name::String
    problems::Vector{String}
end

function Base.showerror(io::IO, e::ContractError)
    print(io, e.table_name, " violates the data contract:")
    for p in e.problems
        print(io, "\n  - ", p)
    end
end

"""
    ColumnSpec

One contract column. `kind` is one of `:string`, `:float64`, `:int8`, `:int`
(Int64) or `:date`. `lo`/`hi` are numeric bounds (strict when `lo_open`/`hi_open`).
`allowed_str`/`allowed_int` restrict values to a set when non-empty.
"""
struct ColumnSpec
    name::String
    kind::Symbol
    nullable::Bool
    lo::Float64
    hi::Float64
    lo_open::Bool
    hi_open::Bool
    allowed_str::Vector{String}
    allowed_int::Vector{Int}
end

function ColumnSpec(name, kind; nullable=false, lo=-Inf, hi=Inf, lo_open=false, hi_open=false,
                    allowed_str=String[], allowed_int=Int[])
    return ColumnSpec(name, kind, nullable, lo, hi, lo_open, hi_open, allowed_str, allowed_int)
end

"""
    TableSpec

A contract table: its columns, an optional unique key, and an optional pair of
date columns `(start, end)` where `end >= start` whenever both are present.
"""
struct TableSpec
    name::String
    columns::Vector{ColumnSpec}
    unique_key::Vector{String}
    ordered_dates::Vector{String}
end

TableSpec(name, columns; unique_key=String[], ordered_dates=String[]) =
    TableSpec(name, columns, unique_key, ordered_dates)

const CLOSURE_TYPES = ["rainfall", "biotoxin", "wastewater", "pollution", "conservation", "other"]

const CONTRACT_SPECS = Dict{String,TableSpec}(
    "sites" => TableSpec("sites", [
        ColumnSpec("site_id", :string),
        ColumnSpec("name", :string),
        ColumnSpec("lat", :float64; lo=-90, hi=90),
        ColumnSpec("lon", :float64; lo=-180, hi=180),
        ColumnSpec("region", :string; allowed_str=["maine", "chesapeake", "great_bay_nh"]),
        ColumnSpec("growing_area", :string; nullable=true),
        ColumnSpec("notes", :string; nullable=true),
    ]; unique_key=["site_id"]),
    "forcing_daily" => TableSpec("forcing_daily", [
        ColumnSpec("site_id", :string),
        ColumnSpec("date", :date),
        ColumnSpec("temp_c", :float64; lo=-2, hi=35, lo_open=true, hi_open=true),
        ColumnSpec("salinity_psu", :float64; lo=0, hi=45),
        ColumnSpec("chl_ug_l", :float64; nullable=true, lo=0),
        ColumnSpec("pom_mg_l", :float64; nullable=true, lo=0),
        ColumnSpec("source", :string),
        ColumnSpec("qc_flag", :int8; allowed_int=[0, 1, 2]),
    ]; unique_key=["site_id", "date"]),
    "closures" => TableSpec("closures", [
        ColumnSpec("growing_area", :string),
        ColumnSpec("area_id", :string; nullable=true),
        ColumnSpec("start_date", :date),
        ColumnSpec("end_date", :date; nullable=true),
        ColumnSpec("closure_type", :string; allowed_str=CLOSURE_TYPES),
        ColumnSpec("source", :string),
    ]; ordered_dates=["start_date", "end_date"]),
    "growth_obs" => TableSpec("growth_obs", [
        ColumnSpec("dataset_id", :string),
        ColumnSpec("species", :string; allowed_str=["C_virginica", "C_gigas"]),
        ColumnSpec("site_id", :string),
        ColumnSpec("cohort_id", :string),
        ColumnSpec("date", :date),
        ColumnSpec("shell_height_mm", :float64; nullable=true, lo=0),
        ColumnSpec("dry_tissue_g", :float64; nullable=true, lo=0),
        ColumnSpec("n", :int; lo=0),
    ]),
    "disease_obs" => TableSpec("disease_obs", [
        ColumnSpec("site_id", :string),
        ColumnSpec("date", :date),
        ColumnSpec("pathogen", :string; allowed_str=["P_marinus", "H_nelsoni"]),
        ColumnSpec("prevalence", :float64; lo=0, hi=1),
        ColumnSpec("intensity_mean", :float64; nullable=true, lo=0),
        ColumnSpec("intensity_scale", :string),
        ColumnSpec("n", :int; nullable=true, lo=0),
        ColumnSpec("source", :string),
    ]),
    "economics" => TableSpec("economics", [
        ColumnSpec("year", :int; lo=1800, hi=2200),
        ColumnSpec("region", :string),
        ColumnSpec("value_usd", :float64; lo=0),
        ColumnSpec("pounds", :float64; nullable=true, lo=0),
        ColumnSpec("count", :float64; nullable=true, lo=0),
        ColumnSpec("source", :string),
    ]),
    "mortality_obs" => TableSpec("mortality_obs", [
        ColumnSpec("site_id", :string),
        ColumnSpec("date", :date),
        ColumnSpec("period_days", :int; nullable=true, lo=1),
        ColumnSpec("mortality_fraction", :float64; lo=0, hi=1),
        ColumnSpec("method", :string),
        ColumnSpec("source", :string),
    ]),
)

function _type_ok(kind::Symbol, T::Type)
    kind === :string && return T <: AbstractString
    kind === :float64 && return T === Float64
    kind === :int8 && return T === Int8
    kind === :int && return T === Int64
    kind === :date && return T === Date
    throw(ArgumentError("unknown column kind $kind"))
end

function _check_column!(problems::Vector{String}, c::ColumnSpec, col::AbstractVector)
    T = nonmissingtype(eltype(col))
    if !_type_ok(c.kind, T)
        push!(problems, "column \"$(c.name)\" has type $T, expected $(c.kind)")
        return problems
    end
    n_missing = count(ismissing, col)
    if !c.nullable && n_missing > 0
        push!(problems, "column \"$(c.name)\" has $n_missing missing value(s) but is not nullable")
    end
    vals = collect(skipmissing(col))
    if c.kind === :float64
        n_nan = count(isnan, vals)
        n_nan > 0 && push!(problems, "column \"$(c.name)\" has $n_nan NaN(s); use missing for missing values")
    end
    if c.kind in (:float64, :int8, :int)
        if isfinite(c.lo)
            n = count(v -> c.lo_open ? v <= c.lo : v < c.lo, vals)
            n > 0 && push!(problems, "column \"$(c.name)\" has $n value(s) not $(c.lo_open ? ">" : ">=") $(c.lo)")
        end
        if isfinite(c.hi)
            n = count(v -> c.hi_open ? v >= c.hi : v > c.hi, vals)
            n > 0 && push!(problems, "column \"$(c.name)\" has $n value(s) not $(c.hi_open ? "<" : "<=") $(c.hi)")
        end
    end
    if !isempty(c.allowed_str)
        bad = sort!(unique(String(v) for v in vals if !(v in c.allowed_str)))
        isempty(bad) || push!(problems, "column \"$(c.name)\" has values outside $(c.allowed_str): $bad")
    end
    if !isempty(c.allowed_int)
        bad = sort!(unique(Int(v) for v in vals if !(v in c.allowed_int)))
        isempty(bad) || push!(problems, "column \"$(c.name)\" has values outside $(c.allowed_int): $bad")
    end
    return problems
end

"""
    validate_table(tbl, spec::TableSpec) -> Nothing

Check any Tables.jl-compatible table against `spec`. Throws a
[`ContractError`](@ref) listing every problem: missing or unexpected columns,
wrong element types, missing values in non-nullable columns, NaNs, values out
of range or outside an allowed set, duplicate keys, and end dates before start
dates.
"""
function validate_table(tbl, spec::TableSpec)
    cols = Tables.columntable(tbl)
    present = String.(collect(keys(cols)))
    expected = [c.name for c in spec.columns]
    problems = String[]
    for name in expected
        name in present || push!(problems, "missing column \"$name\"")
    end
    for name in present
        name in expected || push!(problems, "unexpected column \"$name\"")
    end
    for c in spec.columns
        c.name in present && _check_column!(problems, c, cols[Symbol(c.name)])
    end
    if isempty(problems) && !isempty(spec.unique_key)
        keyrows = collect(zip((cols[Symbol(k)] for k in spec.unique_key)...))
        n_dup = length(keyrows) - length(unique(keyrows))
        n_dup > 0 && push!(problems, "$n_dup duplicate row(s) for key $(spec.unique_key)")
    end
    if isempty(problems) && !isempty(spec.ordered_dates)
        s, e = cols[Symbol(spec.ordered_dates[1])], cols[Symbol(spec.ordered_dates[2])]
        n_bad = count(i -> !ismissing(s[i]) && !ismissing(e[i]) && e[i] < s[i], eachindex(s))
        n_bad > 0 && push!(problems, "$n_bad row(s) with $(spec.ordered_dates[2]) before $(spec.ordered_dates[1])")
    end
    isempty(problems) || throw(ContractError(spec.name, problems))
    return nothing
end

"Validate a `sites` table. See [`validate_table`](@ref)."
validate_sites(tbl) = validate_table(tbl, CONTRACT_SPECS["sites"])
"Validate a `forcing_daily` table. See [`validate_table`](@ref)."
validate_forcing_daily(tbl) = validate_table(tbl, CONTRACT_SPECS["forcing_daily"])
"Validate a `closures` table. See [`validate_table`](@ref)."
validate_closures(tbl) = validate_table(tbl, CONTRACT_SPECS["closures"])
"Validate a `growth_obs` table. See [`validate_table`](@ref)."
validate_growth_obs(tbl) = validate_table(tbl, CONTRACT_SPECS["growth_obs"])
"Validate a `disease_obs` table. See [`validate_table`](@ref)."
validate_disease_obs(tbl) = validate_table(tbl, CONTRACT_SPECS["disease_obs"])
"Validate an `economics` table. See [`validate_table`](@ref)."
validate_economics(tbl) = validate_table(tbl, CONTRACT_SPECS["economics"])
"Validate a `mortality_obs` table. See [`validate_table`](@ref)."
validate_mortality_obs(tbl) = validate_table(tbl, CONTRACT_SPECS["mortality_obs"])

"""
    read_contract(path, table_name) -> DataFrame

Read an Arrow IPC contract file, validate it as `table_name` (e.g.
`"forcing_daily"`), and return it as a `DataFrame` with copied columns. The
file is read fully into memory so it is not left memory-mapped.
"""
function read_contract(path::AbstractString, table_name::AbstractString)
    haskey(CONTRACT_SPECS, table_name) || throw(ArgumentError("unknown contract table \"$table_name\""))
    df = DataFrame(Arrow.Table(read(path)); copycols=true)
    validate_table(df, CONTRACT_SPECS[table_name])
    return df
end
