
using PoincareInvariants

export loop_odeproblem,
       loop_iodeproblem,
       loop_lodeproblem

export poincare_invariant_1st,
       loop_ensemble


# The base problems whose flow advects the loop, at the loop's magnetic moment. Their initial
# condition is a placeholder: `PIEnsembleProblem` replaces it with the sampled points of `f_loop`,
# one ensemble member per point. Everything else — time span, time step — is taken from the
# problem, so pass those here.
loop_odeproblem(; kwargs...) =
    odeproblem(f_loop(0.0); parameters = (μ = μ_loop(),), periodic = false, kwargs...)

loop_iodeproblem(; kwargs...) =
    iodeproblem(f_loop(0.0); parameters = (μ = μ_loop(),), periodic = false, kwargs...)

loop_lodeproblem(; kwargs...) =
    lodeproblem(f_loop(0.0); parameters = (μ = μ_loop(),), periodic = false, kwargs...)


@doc raw"""
    poincare_invariant_1st(N; DT = Float64)

Set up the computation of the first Poincaré integral invariant

```math
I_{1} = \oint_{\gamma} \vartheta_{i} (q) \, dq^{i}
```

of the four-dimensional guiding centre dynamics, sampling the loop `f_loop` of this equilibrium at
`N` points.

The guiding centre one-form ``\vartheta = A + u b`` is state-dependent, so this is a *noncanonical*
`FirstPoincareInvariant` built from `ϑ`, not the canonical variant. Pair it with one of
`loop_odeproblem`, `loop_iodeproblem` or `loop_lodeproblem` through
[`loop_ensemble`](@ref):

```julia
pinv = poincare_invariant_1st(1000)
prob = loop_iodeproblem(; timespan = (0.0, 1E3), timestep = 1.0)
sol  = integrate(loop_ensemble(prob, pinv), VPRKGauss(2))
I₁   = compute!(pinv, sol, parameters(prob))
```
"""
poincare_invariant_1st(N; DT = Float64) = FirstPI{DT, 4}(ϑ, N)


"""
    loop_ensemble(prob, pinv)

Sample the loop `f_loop` of this equilibrium with `pinv` and turn each point into the initial
condition of one member of a `GeometricEquations.EnsembleProblem` built from `prob`. Integrate the
result with `GeometricIntegrators.integrate` and pass the solution to
`PoincareInvariants.compute!`.
"""
loop_ensemble(prob, pinv) = PIEnsembleProblem(prob, pinv, f_loop)
