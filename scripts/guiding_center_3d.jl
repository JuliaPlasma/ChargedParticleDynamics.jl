using GeometricIntegrators
using LinearAlgebra
using SimpleSolvers

using ChargedParticleDynamics
using ChargedParticleDynamics: InitialConditions, charged_particle, guiding_center, pauli_particle
using ChargedParticleDynamics: md

import ChargedParticleDynamics.ChargedParticle3d
import ChargedParticleDynamics.GuidingCenter3d
import ChargedParticleDynamics.GuidingCenter4d
import ChargedParticleDynamics.PauliParticle3d

import .ChargedParticle3d.TokamakIterCylindrical: from_cartesian, g₁₁, g₂₂, g₃₃
import .ChargedParticle3d.TokamakIterCylindrical: R₀, B, b, aₚ, bₚ, cₚ, b⃗, ḡ, DF̄, J


X₀ = from_cartesian(0, [7.0, 0, 0])
E₀ = 1E6
θ₀ = 0.
α₀ = π*5/16
m₀ = md

ics = InitialConditions(X₀, θ₀, α₀, E₀, m₀, 1, aₚ, bₚ, cₚ, b⃗, B, ḡ, DF̄, J)

x₀, v₀ = charged_particle(ics)
p₀ = ChargedParticle3d.TokamakIterCylindrical.charged_particle_3d_pᵢ(0, x₀, v₀)

q₀, μ = guiding_center(ics)

q3₀ = copy(q₀)
q3₀[2] = sqrt(eps())
q3₀

parameters = (μ=μ,)

options = (f_abstol = 1E-14, max_iterations = 10_000)#, linesearch = Bisection()

Δt1 = 0.01
Δt2 = 0.1
Δt3 = 1.0
Δt4 = 10.0
Δt5 = 50.0
# tspan = (0., 5_000.)
# tspan = (0., 10_000.)
tspan = (0., 50_000.)
tspan_long = (0., 1_000_000.);


g3ode = GuidingCenter3d.TokamakIterCylindrical.hode(q3₀, parameters; tstep=Δt3, tspan=tspan)#(0., 19_189.)
g3sol = integrate(g3ode, PartitionedGauss(1); initialguess=NoInitialGuess(),options...)
# g3sol = integrate(g3ode, PartitionedGauss(1); options...)
g3car = cartesian_solution(g3sol, GuidingCenter3d.TokamakIterCylindrical);
