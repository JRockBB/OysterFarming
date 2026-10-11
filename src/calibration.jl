# Calibration of selected DEB parameters to shell-height observations (PLAN.md P2.5).
#
# Two methods share one problem definition:
#   fit_least_squares  - bounded L-BFGS (Optimization.jl + Optim.jl) on the sum of
#                        squared errors, with Gauss-Newton standard errors.
#   fit_bayesian       - NUTS sampling (Turing.jl) with uniform priors on the
#                        parameter bounds and Gaussian observation error.
# Both differentiate through the ODE solver with ForwardDiff.

"""
    CalibrationProblem(base, forcing, u0, obs_dates, obs_shell_height_mm;
                       names, lower, upper)

Fit the `DEBParams` fields listed in `names` (e.g. `[:X_K, :p_M]`) within
`lower`/`upper`, keeping every other field of `base` fixed. `forcing` is a
single-site run of consecutive days as accepted by [`simulate_individual`](@ref);
observations are shell heights (mm) on dates inside the forcing.
"""
struct CalibrationProblem{P<:DEBParams,F<:AbstractDataFrame,U<:NamedTuple}
    base::P
    forcing::F
    u0::U
    obs_idx::Vector{Int}
    obs_shell_height_mm::Vector{Float64}
    names::Vector{Symbol}
    lower::Vector{Float64}
    upper::Vector{Float64}
end

function CalibrationProblem(base::DEBParams, forcing, u0::NamedTuple,
                            obs_dates::AbstractVector{Date}, obs_shell_height_mm::AbstractVector{<:Real};
                            names::AbstractVector{Symbol}, lower::AbstractVector{<:Real},
                            upper::AbstractVector{<:Real})
    df = DataFrame(forcing; copycols=false)
    length(obs_dates) == length(obs_shell_height_mm) ||
        throw(ArgumentError("obs_dates and obs_shell_height_mm differ in length"))
    length(names) == length(lower) == length(upper) ||
        throw(ArgumentError("names, lower and upper must have the same length"))
    all(lower .< upper) || throw(ArgumentError("each lower bound must be below its upper bound"))
    for n in names
        n in fieldnames(DEBParams) || throw(ArgumentError("DEBParams has no field $n"))
    end
    pos = Dict(d => i for (i, d) in enumerate(df.date))
    missing_dates = [d for d in obs_dates if !haskey(pos, d)]
    isempty(missing_dates) || throw(ArgumentError("observation dates outside forcing: $missing_dates"))
    idx = [pos[d] for d in obs_dates]
    return CalibrationProblem(base, df, u0, idx, Float64.(obs_shell_height_mm), collect(names),
                              Float64.(lower), Float64.(upper))
end

"""
    predict_shell_height(cp::CalibrationProblem, θ) -> Vector

Simulated shell height (mm) on the observation dates with the fitted
parameters set to `θ` (same order as `cp.names`).
"""
function predict_shell_height(cp::CalibrationProblem, θ::AbstractVector)
    p = with_values(cp.base; (cp.names[i] => θ[i] for i in eachindex(cp.names))...)
    out = simulate_individual(p, cp.forcing, cp.u0)
    return out.shell_height_mm[cp.obs_idx]
end

"Sum of squared shell-height residuals, mm²."
sum_squared_error(cp::CalibrationProblem, θ::AbstractVector) =
    sum(abs2, predict_shell_height(cp, θ) .- cp.obs_shell_height_mm)

"""
    fit_least_squares(cp, θ0; z=1.959964) -> NamedTuple

Minimise the sum of squared errors within the bounds (L-BFGS in a box,
ForwardDiff gradients). Approximate standard errors come from the
Gauss-Newton covariance σ̂²(JᵀJ)⁻¹, with σ̂² = SSE/(n − k) and J the Jacobian
of predictions at the optimum (e.g. Bates & Watts 1988, Nonlinear Regression
Analysis, ch. 2). Intervals are estimate ± z·SE (default 95 %).

Returns `(names, estimate, se, lower_ci, upper_ci, sigma_mm, sse, converged)`.
"""
function fit_least_squares(cp::CalibrationProblem, θ0::AbstractVector{<:Real}; z::Real=1.959964)
    all(cp.lower .<= θ0 .<= cp.upper) || throw(ArgumentError("θ0 must lie within the bounds"))
    optf = OptimizationFunction((θ, c) -> sum_squared_error(c, θ), AutoForwardDiff())
    prob = OptimizationProblem(optf, Float64.(θ0), cp; lb=cp.lower, ub=cp.upper)
    sol = solve(prob, LBFGS())
    θ̂ = collect(sol.u)
    J = ForwardDiff.jacobian(θ -> predict_shell_height(cp, θ), θ̂)
    n, k = length(cp.obs_shell_height_mm), length(θ̂)
    n > k || throw(ArgumentError("need more observations than fitted parameters"))
    sse = sum_squared_error(cp, θ̂)
    σ² = sse / (n - k)
    cov = σ² .* inv(Symmetric(J' * J))
    se = sqrt.(diag(cov))
    return (names=cp.names, estimate=θ̂, se=se, lower_ci=θ̂ .- z .* se, upper_ci=θ̂ .+ z .* se,
            sigma_mm=sqrt(σ²), sse=sse, converged=SciMLBase.successful_retcode(sol))
end

@model function _deb_turing_model(cp, y, sigma_prior_scale_mm)
    θ ~ arraydist([Uniform(cp.lower[i], cp.upper[i]) for i in eachindex(cp.names)])
    σ ~ truncated(Normal(0, sigma_prior_scale_mm); lower=0)
    pred = predict_shell_height(cp, θ)
    y ~ MvNormal(pred, σ^2 * I)
end

"""
    fit_bayesian(rng, cp; sigma_prior_scale_mm, n_samples=500, n_adapt=500,
                 target_accept=0.8, level=0.95) -> NamedTuple

Bayesian calibration with Turing.jl NUTS. Model:

    θ_i ~ Uniform(lower_i, upper_i)
    σ   ~ HalfNormal(sigma_prior_scale_mm)
    y   ~ MvNormal(shell_height(θ), σ² I)

`sigma_prior_scale_mm` is a modelling choice the caller must state. Returns
`(names, mean, lower_ci, upper_ci, sigma_mean_mm, chain)` where the intervals
are equal-tailed credible intervals at `level`.
"""
function fit_bayesian(rng::AbstractRNG, cp::CalibrationProblem; sigma_prior_scale_mm::Real,
                      n_samples::Integer=500, n_adapt::Integer=500, target_accept::Real=0.8,
                      level::Real=0.95)
    model = _deb_turing_model(cp, cp.obs_shell_height_mm, sigma_prior_scale_mm)
    chain = sample(rng, model, NUTS(n_adapt, target_accept; adtype=AutoForwardDiff()), n_samples;
                   progress=false)
    k = length(cp.names)
    θdraws = vec(collect(chain[:θ]))   # one parameter vector per draw (FlexiChains, Turing 0.49)
    draws = [[d[i] for d in θdraws] for i in 1:k]
    α = (1 - level) / 2
    return (names=cp.names, mean=mean.(draws),
            lower_ci=[quantile(d, α) for d in draws], upper_ci=[quantile(d, 1 - α) for d in draws],
            sigma_mean_mm=mean(vec(collect(chain[:σ]))), chain=chain)
end
