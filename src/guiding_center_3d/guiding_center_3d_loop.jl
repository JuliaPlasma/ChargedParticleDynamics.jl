
using GeometricEquations: timespan
using PoincareInvariants

export loop_hodeproblem,
       loop_hodeproblem_canonical,
       loop_hodeproblem_compact

export poincare_invariant_1st,
       loop_ensemble

# The base problems whose flow advects the loop, at the loop's magnetic moment. Their initial
# condition is a placeholder: `PIEnsembleProblem` replaces it with the sampled points of `f_loop`,
# one ensemble member per point. Everything else — time span, time step, constraint pair — is
# taken from the problem, so pass those here.
for problem in (:hodeproblem, :hodeproblem_canonical, :hodeproblem_compact)
    @eval function $(Symbol(:loop_, problem))(; kwargs...)
        $problem(f_loop(0.0); parameters = (field = FIELD, μ = μ_loop()), periodic = false,
            kwargs...)
    end
end

@doc raw"""
    poincare_invariant_1st(N; DT = Float64)

Set up the computation of the first Poincaré integral invariant

```math
I_{1} = \oint_{\gamma} p_{i} \, dq^{i}
```

of the three-dimensional guiding centre dynamics, sampling the loop `f_loop` of this equilibrium
at `N` points.

The 3D model is canonical in `(q, p)`, so this is the canonical `FirstPoincareInvariant` on the
six-dimensional phase space. On the constraint manifold `p = ϑ(q, u)`, where
[`loop_ensemble`](@ref) places every point, it is the first invariant of the 4D guiding centre
model on the same loop. Pair it with one of `loop_hodeproblem`, `loop_hodeproblem_canonical` or
`loop_hodeproblem_compact` through [`loop_ensemble`](@ref):

```julia
pinv = poincare_invariant_1st(200)
prob = loop_hodeproblem(; timespan = (0.0, 1E2), timestep = 0.1)
sol  = integrate(loop_ensemble(prob, pinv), PartitionedGauss(2))
I₁   = compute!(pinv, sol)
```
"""
poincare_invariant_1st(N; DT = Float64) = CanonicalFirstPI{DT, 6}(N)

"""
    loop_ensemble(prob, pinv)

Sample the loop `f_loop` of this equilibrium with `pinv`, lift each point `(x, u)` to the
position and momentum `(q, p)` of the guiding centre one-form, and turn it into the initial
condition of one member of a `GeometricEquations.EnsembleProblem` built from `prob`. Integrate
the result with `GeometricIntegrators.integrate` and pass the solution to
`PoincareInvariants.compute!`.
"""
function loop_ensemble(prob, pinv)
    t₀ = timespan(prob)[begin]
    PIEnsembleProblem(prob, pinv, s -> (
        ics = initial_conditions(t₀, f_loop(s)); [ics.q; ics.p]))
end
