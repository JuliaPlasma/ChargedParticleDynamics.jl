
using GeometricEquations: timespan
using PoincareInvariants

export surface_hodeproblem,
       surface_hodeproblem_canonical,
       surface_hodeproblem_compact

export poincare_invariant_2nd,
       surface_ensemble

# The base problems whose flow advects the surface, at the surface's magnetic moment. As for the
# loop, their initial condition is a placeholder that `PIEnsembleProblem` replaces with the
# sampled points of `f_surface`.
for problem in (:hodeproblem, :hodeproblem_canonical, :hodeproblem_compact)
    @eval function $(Symbol(:surface_, problem))(; kwargs...)
        $problem(f_surface(0.0, 0.0); parameters = (field = FIELD, μ = μ_surface()),
            periodic = false, kwargs...)
    end
end

@doc raw"""
    poincare_invariant_2nd(N; DT = Float64, plan = SecondChebyshevPlan)

Set up the computation of the second Poincaré integral invariant

```math
I_{2} = \int_{S} \omega_{ij} (z) \, dz^{i} \wedge dz^{j}
```

of the three-dimensional guiding centre dynamics, sampling the surface `f_surface` of this
equilibrium with the point specification `N`.

The 3D model is canonical in `z = (q, p)`, so the two-form is constant. It is taken in the
convention of the 4D model's `ω`, ``\omega_{ij} = \partial_{j} \vartheta_{i} - \partial_{i}
\vartheta_{j}`` with ``\vartheta = (p, 0)``, which is the negative of
`PoincareInvariants.CanonicalSymplecticMatrix`. On the constraint manifold, where
[`surface_ensemble`](@ref) places every point, this is then the second invariant of the 4D
guiding centre model on the same surface, sign included. Pass `N` as a tuple `(nx, ny)` together
with `plan = SecondFinDiffPlan` for a grid layout.

```julia
pinv = poincare_invariant_2nd(351)
prob = surface_hodeproblem(; timespan = (0.0, 2E1), timestep = 0.1)
sol  = integrate(surface_ensemble(prob, pinv), PartitionedGauss(2))
I₂   = compute!(pinv, sol)
```
"""
function poincare_invariant_2nd(N; DT = Float64, plan = PoincareInvariants.SecondChebyshevPlan)
    SecondPI{DT, 6}(-CanonicalSymplecticMatrix{DT}(6), N, plan)
end

"""
    surface_ensemble(prob, pinv)

Sample the surface `f_surface` of this equilibrium with `pinv`, lift each point `(x, u)` to the
position and momentum `(q, p)` of the guiding centre one-form, and turn it into the initial
condition of one member of a `GeometricEquations.EnsembleProblem` built from `prob`. Integrate
the result with `GeometricIntegrators.integrate` and pass the solution to
`PoincareInvariants.compute!`.
"""
function surface_ensemble(prob, pinv)
    t₀ = timespan(prob)[begin]
    PIEnsembleProblem(prob, pinv,
        (s, t) -> (ics = initial_conditions(t₀, f_surface(s, t)); [ics.q; ics.p]))
end
