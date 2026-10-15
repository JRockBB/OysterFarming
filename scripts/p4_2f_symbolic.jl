# P4.2f symbolic step: distil the learned transmission multiplier g_θ(S, T) of the UDE (scripts/p4_2f_ude_fit.jl)
# into a closed-form expression with symbolic regression, then test that expression in the full model.
#
# Design (fixed before running):
# - Targets: g_θ from outputs/p4_2f/ude_weights_full.csv on a grid S = 0:1:35 ppt, T = 0:1:30 degC, weighted by how
#   often each (rounded S, T) pair occurs in the training forcing (Maryland bars 1996-2024, Great Bay beds 1995-2011);
#   pairs never seen get weight 0.
# - SymbolicRegression.jl: operators + - * / and exp, maximum size 15, 60 iterations, fixed seed. Pareto front of
#   complexity vs weighted loss reported; chosen formula = simplest within 5 % of the best loss.
# - Test: the chosen formula replaces the network (transmission_mod). Its own constants are kept from the symbolic
#   fit; the five mechanistic parameters (log beta, log mu, logit host, log epsilon MD, log epsilon GB) are refit by
#   Nelder-Mead per fold, with the same folds, arrival years and pass rule as P4.2e/P4.2f. Limitation: the formula's
#   form comes from the full-data network, so held-out data influenced its shape (not its fitted constants).
#
# Usage: julia --project=scripts/sciml scripts/p4_2f_symbolic.jl
# Writes outputs/p4_2f/symbolic.md.

using OysterFarming
using Arrow, DataFrames, Dates, Random, Statistics
using SymbolicRegression
using Optimization, OptimizationOptimJL

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "outputs", "p4_2f")
include(joinpath(@__DIR__, "p4_2f_common.jl"))     # data assembly, g_net, prediction helpers shared with the UDE fit

θ̂ = [parse(Float64, split(l, ",")[2]) for l in readlines(joinpath(OUT, "ude_weights_full.csv"))[2:end]]
w = nn_weights(θ̂)

# ---- training-forcing frequency weights ----
counts = zeros(36, 31)
for s in data
    for (t, sal) in zip(s.T, s.S)
        i = clamp(round(Int, sal), 0, 35) + 1; j = clamp(round(Int, t), 0, 30) + 1
        counts[i, j] += 1
    end
end
Sg = repeat(0.0:35.0, 1, 31); Tg = repeat((0.0:30.0)', 36, 1)
keep = vec(counts) .> 0
X = permutedims(hcat(vec(Sg)[keep], vec(Tg)[keep]))          # 2 x n, rows S, T
y = [g_net(w, X[2, k], X[1, k]) for k in 1:size(X, 2)]
wts = vec(counts)[keep] ./ sum(counts)

options = SymbolicRegression.Options(binary_operators=[+, -, *, /], unary_operators=[exp], maxsize=15,
                                     seed=20261014, deterministic=true, save_to_file=false)
hof = equation_search(X, y; weights=wts, options=options, niterations=60, variable_names=["S", "T"], parallelism=:serial)
front = calculate_pareto_frontier(hof)
best_loss = minimum(m.loss for m in front)
chosen = first(m for m in front if m.loss <= 1.05 * best_loss)
expr_str(m) = string_tree(m.tree, options; variable_names=["S", "T"])
gsym(t, s) = max(0.0, first(eval_tree_array(chosen.tree, reshape([s, t], 2, 1), options))[1])

# ---- refit mechanistic parameters with the symbolic term, per fold ----
function objective_sym(φ, sites, arrival)
    θ = vcat(φ, w)          # network weights unused below; keep layout for params_from
    md = 0.0; nmd = 0; gb = 0.0; ngb = 0
    for s in sites
        p = params_from(θ, s.region)
        r = simulate_establishment(p, s.T, s.S, initial(s), 0.0; intro_day=intro_day(arrival, s),
                                   background_mortality_d=BACKGROUND_D, record=s.rec, transmission_mod=gsym)
        e = sum(abs2, r.prev .- s.obs)
        s.region == "gb" ? (gb += e; ngb += length(s.obs)) : (md += e; nmd += length(s.obs))
    end
    return (nmd > 0 ? md / nmd : 0.0) + (ngb > 0 ? gb / ngb : 0.0)
end
function fit_sym(sites, arrival)
    prob = OptimizationProblem((φ, _) -> objective_sym(φ, sites, arrival), θ̂[1:5])
    return solve(prob, NelderMead(); maxiters=800).u
end
function score_sym(φ, test, arrival)
    θ = vcat(φ, w)
    pred = Float64[]; obs = Float64[]; pers = Union{Missing,Float64}[]
    for s in test
        p = params_from(θ, s.region)
        r = simulate_establishment(p, s.T, s.S, initial(s), 0.0; intro_day=intro_day(arrival, s),
                                   background_mortality_d=BACKGROUND_D, record=s.rec, transmission_mod=gsym)
        append!(pred, r.prev); append!(obs, s.obs); append!(pers, s.prev_obs)
    end
    ok = .!ismissing.(pers)
    return (n=length(obs), rmse=rmse(pred, obs), cor=length(obs) > 2 ? cor(pred, obs) : NaN,
            rmse_model=rmse(pred[ok], obs[ok]), rmse_pers=rmse(Float64.(pers[ok]), obs[ok]), n_pers=count(ok))
end

rows = DataFrame(heldout=String[], n=Int[], rmse=Float64[], cor=Float64[], rmse_model=Float64[], rmse_pers=Float64[], n_pers=Int[])
for s in [x for x in data if x.region == "gb"]
    φ = fit_sym([x for x in data if x.id != s.id], ARRIVAL[s.id])
    sc = score_sym(φ, [s], ARRIVAL[s.id])
    push!(rows, (s.id, sc.n, sc.rmse, sc.cor, sc.rmse_model, sc.rmse_pers, sc.n_pers))
    println("leave out ", s.id, ": ", round(sc.rmse_model; digits=3), " vs ", round(sc.rmse_pers; digits=3)); flush(stdout)
end
for reg in ("md_upper", "md_middle", "md_lower")
    φ = fit_sym([x for x in data if x.region != reg], ARRIVAL[reg])
    sc = score_sym(φ, [x for x in data if x.region == reg], ARRIVAL[reg])
    push!(rows, (reg, sc.n, sc.rmse, sc.cor, sc.rmse_model, sc.rmse_pers, sc.n_pers))
    println("leave out ", reg, ": ", round(sc.rmse_model; digits=3), " vs ", round(sc.rmse_pers; digits=3)); flush(stdout)
end
gbr = rows[startswith.(rows.heldout, "nh_"), :]
wavg(x, v) = sum(x .* v) / sum(v)
gb_model = wavg(gbr.rmse_model, gbr.n_pers); gb_pers = wavg(gbr.rmse_pers, gbr.n_pers)

fmt(x; d=3) = (x isa Real && isnan(x)) ? "n/a" : string(round(x; digits=d))
open(joinpath(OUT, "symbolic.md"), "w") do io
    println(io, "# P4.2f symbolic regression of the learned transmission multiplier (scripts/p4_2f_symbolic.jl)\n")
    println(io, "$(size(X, 2)) grid points with training-forcing weights. Pareto front (complexity, weighted loss, expression):\n")
    println(io, "| complexity | loss | expression |\n|---|---|---|")
    for m in front
        println(io, "| $(compute_complexity(m.tree, options)) | $(round(m.loss; sigdigits=3)) | `$(expr_str(m))` |")
    end
    println(io, "\nChosen (simplest within 5 % of best loss): `$(expr_str(chosen))`.\n")
    Ts = [5.0, 15.0, 22.0, 28.0]
    println(io, "| S (ppt) | " * join(["network T=$(Int(t))" for t in Ts], " | ") * " | " * join(["formula T=$(Int(t))" for t in Ts], " | ") * " |")
    println(io, "|---|" * repeat("---|", 2 * length(Ts)))
    for S in (2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0)
        println(io, "| $(Int(S)) | " * join([fmt(g_net(w, t, S)) for t in Ts], " | ") * " | " * join([fmt(gsym(t, S)) for t in Ts], " | ") * " |")
    end
    println(io, "\n## Held-out tests with the formula (mechanistic parameters refit per fold)\n")
    println(io, "| held out | n | RMSE | Pearson | RMSE model (rows with previous survey) | RMSE persistence (same rows) |\n|---|---|---|---|---|---|")
    for r in eachrow(rows)
        println(io, "| $(r.heldout) | $(r.n) | $(fmt(r.rmse)) | $(fmt(r.cor)) | $(fmt(r.rmse_model)) | $(fmt(r.rmse_pers)) |")
    end
    println(io, "\nGreat Bay held-out, survey-weighted: formula $(fmt(gb_model)) vs persistence $(fmt(gb_pers)). Pass: **$(gb_model < gb_pers)**.")
    println(io, "Limitation: the formula's form was found from the full-data network, so held-out data influenced its shape.")
end
println(read(joinpath(OUT, "symbolic.md"), String))
