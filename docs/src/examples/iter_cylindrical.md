# ITER Equilibrium im Cylindrical Coordinates

## Load Dependencies

```@example 1
using CairoMakie
using GeometricIntegrators
using LaTeXStrings
using LinearAlgebra

using ChargedParticleDynamics
using ChargedParticleDynamics: InitialConditions, charged_particle, guiding_center, pauli_particle
using ChargedParticleDynamics: md

import ChargedParticleDynamics.ChargedParticle3d
import ChargedParticleDynamics.GuidingCenter3d
import ChargedParticleDynamics.GuidingCenter4d
import ChargedParticleDynamics.PauliParticle3d

import .ChargedParticle3d.TokamakIterCylindrical: from_cartesian, g₁₁, g₂₂, g₃₃
import .ChargedParticle3d.TokamakIterCylindrical: R₀, B, b, aₚ, bₚ, cₚ, b⃗, ḡ, DF̄, J
```


## Magnetic Field

```@example 1
fig, ax = plot_fieldlines(ChargedParticle3d.TokamakIterCylindrical; xrange = (3., 11.), yrange = (-4., +4.),  levels = 25)
fig
```


## Initial conditions and solver options

Set and initialize initial conditions:
```@example 1
X₀ = from_cartesian(0, [7.0, 0, 0])
E₀ = 1E6
θ₀ = 0.
α₀ = π*5/16
m₀ = md

ics = InitialConditions(X₀, θ₀, α₀, E₀, m₀, 1, aₚ, bₚ, cₚ, b⃗, B, ḡ, DF̄, J)
```

Obtain charged particle initial conditions:
```@example 1
x₀, v₀ = charged_particle(ics)
p₀ = ChargedParticle3d.TokamakIterCylindrical.charged_particle_3d_pᵢ(0, x₀, v₀)
```

Obtain Pauli particle initial conditions:
```@example 1
X₀, u₀, μ = pauli_particle(ics)
```

Obtain guiding center initial conditions:
```@example 1
q₀, μ = guiding_center(ics)
```

Modify guiding center initial conditions for 3D models to avoid singularities:
```@example 1
q3₀ = copy(q₀)
q3₀[2] = sqrt(eps())
q3₀
```

Parameters:
```@example 1
parameters = (μ=μ,)
```

Time step size and integration intervals:
```@example 1
Δt1 = 0.01
Δt2 = 0.1
Δt3 = 1.0
Δt4 = 10.0
Δt5 = 50.0
tspan = (0., 10_000.)
tspan_long = (0., 1_000_000.);
```

Nonlinear solver options:
```@example 1
options = (f_abstol = 1E-15, max_iterations = 1_000)
```

## Charged Particle

Integrate charged particle dynamics and convert result to cartesian coordinates:
```@example 1
code = ChargedParticle3d.TokamakIterCylindrical.charged_particle_3d_pode(x₀, p₀; tstep=Δt1, tspan=tspan)
csol = integrate(code, PartitionedGauss(1); options...)
ccar = cartesian_solution(csol, ChargedParticle3d.TokamakIterCylindrical);
```

Plot solution and energy error:
```@example 1
fig, ax = plot_trajectory_poloidal(ccar.R, ccar.Z, ChargedParticle3d.TokamakIterCylindrical; linewidth = 1)
fig
```

```@example 1
fig, ax = plot_trajectory_3d(ccar.X, ccar.Y, ccar.Z)
fig
```

```@example 1
cham = compute_invariant(csol.t, csol.q, csol.p, code.parameters, ChargedParticle3d.TokamakIterCylindrical.hamiltonian)
plot_invariant_error(csol.t, cham; label="H", padding_right=30)
```

Integrate charged particle dynamics with large time step and convert result to cartesian coordinates:
```@example 1
code_Δt2 = ChargedParticle3d.TokamakIterCylindrical.charged_particle_3d_pode(x₀, p₀; tstep=Δt2, tspan=tspan)
csol_Δt2 = integrate(code_Δt2, PartitionedGauss(1); options...)
ccar_Δt2 = cartesian_solution(csol_Δt2, ChargedParticle3d.TokamakIterCylindrical);
```

## Pauli Particle with Symplectic Integrator

Integrate charged particle dynamics and convert result to cartesian coordinates:
```@example 1
hode = PauliParticle3d.TokamakIterCylindrical.pauli_particle_3d_hode(X₀, u₀, parameters; tstep=Δt3, tspan=tspan)
hsol = integrate(hode, PartitionedGauss(1); options...)
hcar = cartesian_solution(hsol, PauliParticle3d.TokamakIterCylindrical);
```

Plot solution and energy error:
```@example 1
fig, ax = plot_trajectory_poloidal(ccar.R, ccar.Z, ChargedParticle3d.TokamakIterCylindrical; label="Charged Particle Δt=0.01", linewidth = 1)
plot_trajectory_scatter!(fig, ax, hcar.R, hcar.Z; label="Pauli Particle Δt=1.0")
axislegend(ax)
fig
```

```@example 1
fig, ax = plot_trajectory_3d(ccar.X, ccar.Y, ccar.Z)
plot_trajectory_3d!(fig, ax, hcar.X, hcar.Y, hcar.Z; linewidth = 3, color=:orange)
fig
```

```@example 1
hham = compute_invariant(hsol.t, hsol.q, hsol.p, hode.parameters, PauliParticle3d.TokamakIterCylindrical.hamiltonian)
plot_invariant_error(hsol.t, hham; label="H", padding_right=30)
```

Integrate charged particle dynamics with large time step and convert result to cartesian coordinates:
```@example 1
hode_Δt5 = PauliParticle3d.TokamakIterCylindrical.pauli_particle_3d_hode(X₀, u₀, parameters; tstep=Δt5, tspan=tspan_long)
hsol_Δt5 = integrate(hode_Δt5, PartitionedGauss(1); initialguess=NoInitialGuess(), options...);
hcar_Δt5 = cartesian_solution(hsol_Δt5, PauliParticle3d.TokamakIterCylindrical);
```

## Pauli Particle with Variational Integrator

Integrate charged particle dynamics and convert result to cartesian coordinates:
```@example 1
pode = PauliParticle3d.TokamakIterCylindrical.pauli_particle_3d_iode(X₀, u₀, parameters; tstep=Δt3, tspan=tspan)
psol = integrate(pode, VPRKGauss(1); options...)
pcar = cartesian_solution(psol, PauliParticle3d.TokamakIterCylindrical);
```

Plot solution and energy error:
```@example 1
fig, ax = plot_trajectory_poloidal(ccar.R, ccar.Z, ChargedParticle3d.TokamakIterCylindrical; label="Charged Particle Δt=0.01", linewidth = 1)
plot_trajectory_scatter!(fig, ax, pcar.R, pcar.Z; label="Pauli Particle Δt=1.0")
axislegend(ax)
fig
```

```@example 1
fig, ax = plot_trajectory_3d(ccar.X, ccar.Y, ccar.Z)
plot_trajectory_3d!(fig, ax, pcar.X, pcar.Y, pcar.Z; linewidth = 3, color=:orange)
fig
```

```@example 1
pham = compute_invariant(psol.t, psol.q, psol.p, pode.parameters, PauliParticle3d.TokamakIterCylindrical.hamiltonian)
plot_invariant_error(psol.t, pham; label="H", padding_right=30)
```

Integrate charged particle dynamics with large time step and convert result to cartesian coordinates:
```@example 1
pode_Δt5 = PauliParticle3d.TokamakIterCylindrical.pauli_particle_3d_iode(X₀, u₀, parameters; tstep=Δt5, tspan=tspan_long)
psol_Δt5 = integrate(pode_Δt5, VPRKGauss(1); initialguess=NoInitialGuess(), options...);
pcar_Δt5 = cartesian_solution(psol_Δt5, PauliParticle3d.TokamakIterCylindrical);
```


## Guiding Center 4D with Variational Lobatto-IIIA-IIIB

Integrate charged particle dynamics and convert result to cartesian coordinates:
```@example 1
vode = GuidingCenter4d.TokamakIterCylindrical.guiding_center_4d_iode(q₀, parameters; tstep=10., tspan=(0,1800))
vsol = integrate(vode, VPRKLobattoIIIAIIIB(2); options...)
vcar = cartesian_solution(vsol, GuidingCenter4d.TokamakIterCylindrical);
```

Plot solution and energy error:
```@example 1
fig, ax = plot_trajectory_poloidal(ccar.R, ccar.Z, ChargedParticle3d.TokamakIterCylindrical; label="Charged Particle Δt=0.01", linewidth = 1)
plot_trajectory_scatter!(fig, ax, vcar.R, vcar.Z; label="Guiding Center Δt=1.0", color=:orange, markersize=10)
axislegend(ax)
fig
```

```@example 1
fig, ax = plot_trajectory_3d(vcar.X, vcar.Y, vcar.Z)
fig
```

```@example 1
vham = compute_invariant(vsol.t, vsol.q, vsol.p, vode.parameters, GuidingCenter4d.TokamakIterCylindrical.hamiltonian)
plot_invariant_error(vsol.t, vham; label="H", padding_right=30)
```


## Guiding Center 4D with Projected Variational Midpoint

Integrate charged particle dynamics and convert result to cartesian coordinates:
```@example 1
gode = GuidingCenter4d.TokamakIterCylindrical.guiding_center_4d_iode(q₀, parameters; tstep=Δt3, tspan=tspan)
gsol = integrate(gode, SymmetricProjection(VPRKGauss(1)); options...)
gcar = cartesian_solution(gsol, GuidingCenter4d.TokamakIterCylindrical);
```

Plot solution and energy error:
```@example 1
fig, ax = plot_trajectory_poloidal(ccar.R, ccar.Z, ChargedParticle3d.TokamakIterCylindrical; label="Charged Particle Δt=0.01", linewidth = 1)
plot_trajectory_scatter!(fig, ax, gcar.R, gcar.Z; label="Guiding Center 4D Δt=1.0")
axislegend(ax)
fig
```

```@example 1
fig, ax = plot_trajectory_3d(ccar.X, ccar.Y, ccar.Z)
plot_trajectory_3d!(fig, ax, gcar.X, gcar.Y, gcar.Z; linewidth = 3, color=:orange)
fig
```

```@example 1
gham = compute_invariant(gsol.t, gsol.q, gsol.p, gode.parameters, GuidingCenter4d.TokamakIterCylindrical.hamiltonian)
plot_invariant_error(gsol.t, gham; label="H", padding_right=30)
```

Integrate charged particle dynamics with large time step and convert result to cartesian coordinates:
```@example 1
gode_Δt5 = GuidingCenter4d.TokamakIterCylindrical.guiding_center_4d_iode(q₀, parameters; tstep=Δt5, tspan=tspan_long)
gsol_Δt5 = integrate(gode_Δt5, SymmetricProjection(VPRKGauss(1)); options...);
gcar_Δt5 = cartesian_solution(gsol_Δt5, GuidingCenter4d.TokamakIterCylindrical);
```


## Guiding Centre 3D (Modified Poisson Structure)

Integrate charged particle dynamics and convert result to cartesian coordinates:
```@example 1
g3ode = GuidingCenter3d.TokamakIterCylindrical.hode(q3₀, parameters; tstep=Δt3, tspan=(0., 19_189.))
g3sol = integrate(g3ode, PartitionedGauss(1); initialguess=NoInitialGuess(), options...)
g3car = cartesian_solution(g3sol, GuidingCenter3d.TokamakIterCylindrical);
```

Plot solution and energy error:
```@example 1
fig, ax = plot_trajectory_poloidal(ccar.R, ccar.Z, ChargedParticle3d.TokamakIterCylindrical; label="Charged Particle Δt=0.01", linewidth = 1)
plot_trajectory_scatter!(fig, ax, g3car.R, g3car.Z; label="Guiding Center 3D Δt=1.0", color=:orange)
axislegend(ax)
fig
```

```@example 1
fig, ax = plot_trajectory_3d(g3car.X, g3car.Y, g3car.Z)
fig
```

```@example 1
g3ham = compute_invariant(g3sol.t, g3sol.q, g3sol.p, g3ode.parameters, GuidingCenter3d.TokamakIterCylindrical.hamiltonian)
plot_invariant_error(g3sol.t, g3ham; label="H", padding_right=40)
```

```@example 1
fig = plot_invariant_error(g3sol.t, g3ham; label="H", padding_right=40)
ylims!(-0.1, 0.2)
fig
```

Compute and plot constraints:
```@example 1
g3g1 = compute_invariant(g3sol.t, g3sol.q, g3sol.p, g3ode.parameters, GuidingCenter3d.TokamakIterCylindrical.g₁)
println()
println("g₁(0) = ", g3g1[begin])
println("g₁(T) = ", g3g1[end])
println()
plot_invariant(g3sol.t, g3g1; label="g₁", padding_right=40)
```

```@example 1
g3g2 = compute_invariant(g3sol.t, g3sol.q, g3sol.p, g3ode.parameters, GuidingCenter3d.TokamakIterCylindrical.g₂)
println()
println("g₂(0) = ", g3g2[begin])
println("g₂(T) = ", g3g2[end])
println()
plot_invariant(g3sol.t, g3g2; label="g₂", padding_right=40)
```

Integrate charged particle dynamics with large time step and convert result to cartesian coordinates:
```@example 1
g3ode_Δt4 = GuidingCenter3d.TokamakIterCylindrical.hode(q3₀, parameters; tstep=Δt4, tspan=(0, 100*Δt4))
g3sol_Δt4 = integrate(g3ode_Δt4, PartitionedGauss(1); initialguess=MidpointExtrapolation(5, 1.0), options...)
g3car_Δt4 = cartesian_solution(g3sol_Δt4, GuidingCenter3d.TokamakIterCylindrical);
```

```@example 1
g3ode_Δt5 = GuidingCenter3d.TokamakIterCylindrical.hode(q3₀, parameters; tstep=Δt5, tspan=(0, 100*Δt4))
g3sol_Δt5 = integrate(g3ode_Δt5, PartitionedGauss(1); initialguess=MidpointExtrapolation(5, 1.0), options...)
g3car_Δt5 = cartesian_solution(g3sol_Δt5, GuidingCenter3d.TokamakIterCylindrical);
```

Plot solution, energy error and constraints:
```@example 1
fig, ax = plot_trajectory_poloidal(ccar.R, ccar.Z, ChargedParticle3d.TokamakIterCylindrical; label="Charged Particle Δt=0.01", linewidth = 1)
plot_trajectory_scatter!(fig, ax, g3car.R, g3car.Z; label="Guiding Center 3D Δt=1.0", markersize=12, color=Makie.wong_colors()[6])
plot_trajectory_scatter!(fig, ax, g3car_Δt4.R, g3car_Δt4.Z; label="Guiding Center 3D Δt=10.0", markersize=12, color=Makie.wong_colors()[3])
plot_trajectory_scatter!(fig, ax, g3car_Δt5.R, g3car_Δt5.Z; label="Guiding Center 3D Δt=50.0", markersize=12, color=Makie.wong_colors()[2])
axislegend(ax)
fig
```

```@example 1
g3ham_Δt4 = compute_invariant(g3sol_Δt4.t, g3sol_Δt4.q, g3sol_Δt4.p, g3ode_Δt4.parameters, GuidingCenter3d.TokamakIterCylindrical.hamiltonian)
plot_invariant_error(g3sol_Δt4.t, g3ham_Δt4; label="H", padding_right=40)
```

```@example 1
g3g1_Δt4 = compute_invariant(g3sol_Δt4.t, g3sol_Δt4.q, g3sol_Δt4.p, g3ode_Δt4.parameters, GuidingCenter3d.TokamakIterCylindrical.g₁)
plot_invariant(g3sol_Δt4.t, g3g1_Δt4; label="g₁", padding_right=40)
```

```@example 1
g3g2_Δt4 = compute_invariant(g3sol_Δt4.t, g3sol_Δt4.q, g3sol_Δt4.p, g3ode_Δt4.parameters, GuidingCenter3d.TokamakIterCylindrical.g₂)
plot_invariant(g3sol_Δt4.t, g3g2_Δt4; label="g₂", padding_right=40)
```

```@example 1
g3ham_Δt5 = compute_invariant(g3sol_Δt5.t, g3sol_Δt5.q, g3sol_Δt5.p, g3ode_Δt5.parameters, GuidingCenter3d.TokamakIterCylindrical.hamiltonian)
plot_invariant_error(g3sol_Δt5.t, g3ham_Δt5; label="H", padding_right=40)
```

```@example 1
g3g1_Δt5 = compute_invariant(g3sol_Δt5.t, g3sol_Δt5.q, g3sol_Δt5.p, g3ode_Δt5.parameters, GuidingCenter3d.TokamakIterCylindrical.g₁)
plot_invariant(g3sol_Δt5.t, g3g1_Δt5; label="g₁", padding_right=40)
```

```@example 1
g3g2_Δt5 = compute_invariant(g3sol_Δt5.t, g3sol_Δt5.q, g3sol_Δt5.p, g3ode_Δt5.parameters, GuidingCenter3d.TokamakIterCylindrical.g₂)
plot_invariant(g3sol_Δt5.t, g3g2_Δt5; label="g₂", padding_right=40)
```
