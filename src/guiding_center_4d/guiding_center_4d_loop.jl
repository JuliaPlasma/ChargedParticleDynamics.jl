
using PoincareInvariants

export guiding_center_4d_loop_ode,
       guiding_center_4d_loop_iode,
       guiding_center_4d_loop_lode

export guiding_center_4d_poincare_invariant_1st,
       guiding_center_4d_loop_ensemble


# The base problems whose flow advects the loop, at the loop's magnetic moment. Their initial
# condition is a placeholder: `PIEnsembleProblem` replaces it with the sampled points of `f_loop`,
# one ensemble member per point. Everything else — time span, time step — is taken from the
# problem, so pass those here.
guiding_center_4d_loop_ode(; kwargs...) =
    guiding_center_4d_ode(f_loop(0.0), (μ = μ_loop(),); periodic = false, kwargs...)

guiding_center_4d_loop_iode(; kwargs...) =
    guiding_center_4d_iode(f_loop(0.0), (μ = μ_loop(),); periodic = false, kwargs...)

guiding_center_4d_loop_lode(; kwargs...) =
    guiding_center_4d_lode(f_loop(0.0), (μ = μ_loop(),); periodic = false, kwargs...)


@doc raw"""
    guiding_center_4d_poincare_invariant_1st(N; DT = Float64)

Set up the computation of the first Poincaré integral invariant

```math
I_{1} = \oint_{\gamma} \vartheta_{i} (q) \, dq^{i}
```

of the four-dimensional guiding centre dynamics, sampling the loop `f_loop` of this equilibrium at
`N` points.

The guiding centre one-form ``\vartheta = A + u b`` is state-dependent, so this is a *noncanonical*
`FirstPoincareInvariant` built from `ϑ`, not the canonical variant. Pair it with one of
`guiding_center_4d_loop_ode`, `..._iode` or `..._lode` through
[`guiding_center_4d_loop_ensemble`](@ref):

```julia
pinv = guiding_center_4d_poincare_invariant_1st(1000)
prob = guiding_center_4d_loop_iode(; tspan = (0.0, 1E3), tstep = 1.0)
sol  = integrate(guiding_center_4d_loop_ensemble(prob, pinv), VPRKGauss(2))
I₁   = compute!(pinv, sol, parameters(prob))
```
"""
guiding_center_4d_poincare_invariant_1st(N; DT = Float64) = FirstPI{DT, 4}(ϑ, N)


"""
    guiding_center_4d_loop_ensemble(prob, pinv)

Sample the loop `f_loop` of this equilibrium with `pinv` and turn each point into the initial
condition of one member of a `GeometricEquations.EnsembleProblem` built from `prob`. Integrate the
result with `GeometricIntegrators.integrate` and pass the solution to
`PoincareInvariants.compute!`.
"""
guiding_center_4d_loop_ensemble(prob, pinv) = PIEnsembleProblem(prob, pinv, f_loop)
