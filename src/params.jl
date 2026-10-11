# Parameter files: every model parameter carries its units and a citation.
# Schema documented in params/README.md.

"""
    Param{T}

One model parameter: `value`, `units`, `source` (full citation) and free-text
`notes`. Constructed only by [`load_params`](@ref), which refuses unsourced values.
"""
struct Param{T}
    value::T
    units::String
    source::String
    notes::String
end

"""
    ParamError(path, problems)

Thrown by [`load_params`](@ref) when a parameter file is incomplete. Lists every problem.
"""
struct ParamError <: Exception
    path::String
    problems::Vector{String}
end

function Base.showerror(io::IO, e::ParamError)
    print(io, "invalid parameter file ", e.path, ":")
    for p in e.problems
        print(io, "\n  - ", p)
    end
end

const _PARAM_KEYS = ("value", "units", "source", "notes")

_is_unsourced(s::AbstractString) = isempty(strip(s)) || occursin("TODO(source)", s)

"""
    load_params(path) -> Dict{String,Param}

Read a TOML parameter file in which every top-level table is one parameter:

```toml
[p_Am]
value = 1.0
units = "J cm^-2 d^-1"
source = "Author et al. Year, Journal Vol:pages, Table N"
notes = "optional"
```

Refuses (throws [`ParamError`](@ref)) when any parameter lacks `value`, `units`
or `source`, has an empty source, or has a source containing `TODO(source)`.
Top-level keys that are not tables (e.g. `title`, `species`) are file metadata
and are ignored. Numeric values are returned as `Float64`.
"""
function load_params(path::AbstractString)
    raw = TOML.parsefile(path)
    problems = String[]
    out = Dict{String,Param}()
    for (name, entry) in raw
        entry isa AbstractDict || continue
        for k in keys(entry)
            k in _PARAM_KEYS || push!(problems, "$name: unknown key \"$k\"")
        end
        haskey(entry, "value") || push!(problems, "$name: missing value")
        units = get(entry, "units", nothing)
        units isa AbstractString || push!(problems, "$name: missing units")
        source = get(entry, "source", nothing)
        if !(source isa AbstractString)
            push!(problems, "$name: missing source")
        elseif _is_unsourced(source)
            push!(problems, "$name: source is empty or TODO(source)")
        end
        if haskey(entry, "value") && units isa AbstractString && source isa AbstractString
            v = entry["value"]
            v = v isa Real && !(v isa Bool) ? Float64(v) : v
            out[name] = Param(v, String(units), String(source), String(get(entry, "notes", "")))
        end
    end
    isempty(problems) || throw(ParamError(String(path), sort!(problems)))
    return out
end

"""
    param_value(params, name) -> Float64

Numeric value of parameter `name`, with a clear error if it is absent.
"""
function param_value(params::AbstractDict{String,<:Param}, name::AbstractString)
    haskey(params, name) || throw(KeyError("parameter \"$name\" not in parameter set"))
    v = params[name].value
    v isa Real || throw(ArgumentError("parameter \"$name\" is not numeric"))
    return Float64(v)
end
