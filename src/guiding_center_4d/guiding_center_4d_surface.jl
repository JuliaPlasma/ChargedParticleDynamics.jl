
using PoincareInvariants

export surface_odeproblem,
       surface_iodeproblem,
       surface_lodeproblem

export poincare_invariant_2nd,
       surface_ensemble


# The base problems whose flow advects the surface, at the surface's magnetic moment. As for the
# loop, their initial condition is a placeholder that `PIEnsembleProblem` replaces with the sampled
# points of `f_surface`.
surface_odeproblem(; kwargs...) =
    odeproblem(f_surface(0.0, 0.0); parameters = (μ = μ_surface(),), periodic = false, kwargs...)

surface_iodeproblem(; kwargs...) =
    iodeproblem(f_surface(0.0, 0.0); parameters = (μ = μ_surface(),), periodic = false, kwargs...)

surface_lodeproblem(; kwargs...) =
    lodeproblem(f_surface(0.0, 0.0); parameters = (μ = μ_surface(),), periodic = false, kwargs...)


@doc raw"""
    poincare_invariant_2nd(N; DT = Float64, plan = SecondChebyshevPlan)

Set up the computation of the second Poincaré integral invariant

```math
I_{2} = \int_{S} \omega_{ij} (q) \, dq^{i} \wedge dq^{j}
```

of the four-dimensional guiding centre dynamics, sampling the surface `f_surface` of this
equilibrium with the point specification `N`.

As for the first invariant, the guiding centre two-form ``\omega = d\vartheta`` is state-dependent,
so this is a *noncanonical* `SecondPoincareInvariant` built from `ω`. Pass `N` as a tuple
`(nx, ny)` together with `plan = SecondFinDiffPlan` for a grid layout, which is what
`PoincareInvariants.plot_surface` needs in order to draw the advected surface.

```julia
pinv = poincare_invariant_2nd(2000)
prob = surface_iodeproblem(; timespan = (0.0, 1E3), timestep = 1.0)
sol  = integrate(surface_ensemble(prob, pinv), VPRKGauss(2))
I₂   = compute!(pinv, sol, parameters(prob))
```
"""
poincare_invariant_2nd(N; DT = Float64, plan = PoincareInvariants.SecondChebyshevPlan) =
    SecondPI{DT, 4}(ω, N, plan)


"""
    surface_ensemble(prob, pinv)

Sample the surface `f_surface` of this equilibrium with `pinv` and turn each point into the initial
condition of one member of a `GeometricEquations.EnsembleProblem` built from `prob`. Integrate the
result with `GeometricIntegrators.integrate` and pass the solution to
`PoincareInvariants.compute!`.
"""
surface_ensemble(prob, pinv) = PIEnsembleProblem(prob, pinv, f_surface)
