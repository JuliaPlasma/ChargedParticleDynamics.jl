
using PoincareInvariants

export guiding_center_4d_surface_ode,
       guiding_center_4d_surface_iode,
       guiding_center_4d_surface_lode

export guiding_center_4d_poincare_invariant_2nd,
       guiding_center_4d_surface_ensemble


# The base problems whose flow advects the surface, at the surface's magnetic moment. As for the
# loop, their initial condition is a placeholder that `PIEnsembleProblem` replaces with the sampled
# points of `f_surface`.
guiding_center_4d_surface_ode(; kwargs...) =
    guiding_center_4d_ode(f_surface(0.0, 0.0), (μ = μ_surface(),); periodic = false, kwargs...)

guiding_center_4d_surface_iode(; kwargs...) =
    guiding_center_4d_iode(f_surface(0.0, 0.0), (μ = μ_surface(),); periodic = false, kwargs...)

guiding_center_4d_surface_lode(; kwargs...) =
    guiding_center_4d_lode(f_surface(0.0, 0.0), (μ = μ_surface(),); periodic = false, kwargs...)


@doc raw"""
    guiding_center_4d_poincare_invariant_2nd(N; DT = Float64, plan = SecondChebyshevPlan)

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
pinv = guiding_center_4d_poincare_invariant_2nd(2000)
prob = guiding_center_4d_surface_iode(; tspan = (0.0, 1E3), tstep = 1.0)
sol  = integrate(guiding_center_4d_surface_ensemble(prob, pinv), VPRKGauss(2))
I₂   = compute!(pinv, sol, parameters(prob))
```
"""
guiding_center_4d_poincare_invariant_2nd(N; DT = Float64, plan = PoincareInvariants.SecondChebyshevPlan) =
    SecondPI{DT, 4}(ω, N, plan)


"""
    guiding_center_4d_surface_ensemble(prob, pinv)

Sample the surface `f_surface` of this equilibrium with `pinv` and turn each point into the initial
condition of one member of a `GeometricEquations.EnsembleProblem` built from `prob`. Integrate the
result with `GeometricIntegrators.integrate` and pass the solution to
`PoincareInvariants.compute!`.
"""
guiding_center_4d_surface_ensemble(prob, pinv) = PIEnsembleProblem(prob, pinv, f_surface)
