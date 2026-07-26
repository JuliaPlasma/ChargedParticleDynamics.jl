using GeometricIntegrators
using SimpleSolvers: Options

using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCartesian
using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCartesian: hamiltonian, hamiltonian_u, g₁, g₂, λₒ, λ₁, λ₂, b₁, b₂, b₃
using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCartesian: dg₁dq₁, dg₁dq₂, dg₁dq₃, dg₁dp₁, dg₁dp₂, dg₁dp₃
using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCartesian: dg₂dp₁, dg₂dp₂, dg₂dp₃, dg₂dq₁, dg₂dq₂, dg₂dq₃

# using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCylindrical
# using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCylindrical: hamiltonian, hamiltonian_u, g₁, g₂, λₒ, λ₁, λ₂, b₁, b₂, b₃
# using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCylindrical: dg₁dq₁, dg₁dq₂, dg₁dq₃, dg₁dp₁, dg₁dp₂, dg₁dp₃
# using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCylindrical: dg₂dp₁, dg₂dp₂, dg₂dp₃, dg₂dq₁, dg₂dq₂, dg₂dq₃

# using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallToroidal
# using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallToroidal: hamiltonian, hamiltonian_u, g₁, g₂, λₒ, λ₁, λ₂, b₁, b₂, b₃
# using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallToroidal: dg₁dq₁, dg₁dq₂, dg₁dq₃, dg₁dp₁, dg₁dp₂, dg₁dp₃
# using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallToroidal: dg₂dp₁, dg₂dp₂, dg₂dp₃, dg₂dq₁, dg₂dq₂, dg₂dq₃

# `f_abstol` is an absolute bound on the residual and has to stay above its round-off floor, which
# for these models is `‖ϑ‖ eps`. The 1E-14 this script used to ask for is below it, so the solver
# had no reachable criterion and the line search ran to its 1000-iteration limit on every step.
# `f_reltol` is left at the `SimpleSolvers` default deliberately; see `test/guiding_center_3d_tests.jl`.
const options = (f_abstol=1E-12, max_iterations=50, warn_iterations=50)

equ = hodeproblem(initial_conditions_barely_passing())
# equ = hodeproblem(initial_conditions_barely_trapped())
# equ = hodeproblem(initial_conditions_deeply_passing())
# equ = hodeproblem(initial_conditions_deeply_trapped())


sol = integrate(equ, PartitionedGauss(1); options...)

h = [hamiltonian(sol.t[i], sol.q[i], sol.p[i], parameters(equ)) for i in eachindex(sol.t)]
hu = [hamiltonian_u(sol.t[i], sol.q[i], sol.p[i], parameters(equ)) for i in eachindex(sol.t)]
λ0 = [λₒ(sol.t[i], sol.q[i], sol.p[i]) for i in eachindex(sol.t)]
λ1 = [λ₁(sol.t[i], sol.q[i], sol.p[i], parameters(equ)) for i in eachindex(sol.t)]
λ2 = [λ₂(sol.t[i], sol.q[i], sol.p[i], parameters(equ)) for i in eachindex(sol.t)]
g1 = [g₁(sol.t[i], sol.q[i], sol.p[i]) for i in eachindex(sol.t)]
g2 = [g₂(sol.t[i], sol.q[i], sol.p[i]) for i in eachindex(sol.t)]
b1 = [b₁(sol.t[i], sol.q[i]) for i in eachindex(sol.t)]
b2 = [b₂(sol.t[i], sol.q[i]) for i in eachindex(sol.t)]
b3 = [b₃(sol.t[i], sol.q[i]) for i in eachindex(sol.t)]

println()
println("q(0) = ", sol.q[begin])
println("q(T) = ", sol.q[end])
println()
println("p(0) = ", sol.p[begin])
println("p(T) = ", sol.p[end])
println()
println("H(0) = ", h[begin], ", ", hu[begin])
println("H(T) = ", h[end], ", ", hu[end])
println()
println("g₁(0) = ", g1[begin])
println("g₁(T) = ", g1[end])
println()
println("g₂(0) = ", g2[begin])
println("g₂(T) = ", g2[end])
println()
println("λ₀(0) = ", λ0[begin])
println("λ₀(T) = ", λ0[end])
println()
println("λ₁(0) * g₁(0) + λ₂(0) * g₂(0) = ", λ1[begin] * g1[begin] + λ2[begin] * g2[begin])
println("λ₁(T) * g₁(T) + λ₂(T) * g₂(T) = ", λ1[end] * g1[end] + λ2[end] * g2[end])
println()

println()
println("b₁(x₀) = ", b1[begin])
println("b₂(x₀) = ", b2[begin])
println("b₃(x₀) = ", b3[begin])
println()

println()
println("dg₁dq₁(0) * dg₂dp₁(0) = ", dg₁dq₁(sol.t[begin], sol.q[begin], sol.p[begin]) * dg₂dp₁(sol.t[begin], sol.q[begin], sol.p[begin]))
println("dg₁dq₂(0) * dg₂dp₂(0) = ", dg₁dq₂(sol.t[begin], sol.q[begin], sol.p[begin]) * dg₂dp₂(sol.t[begin], sol.q[begin], sol.p[begin]))
println("dg₁dq₃(0) * dg₂dp₃(0) = ", dg₁dq₃(sol.t[begin], sol.q[begin], sol.p[begin]) * dg₂dp₃(sol.t[begin], sol.q[begin], sol.p[begin]))
println("dg₁dp₁(0) * dg₂dq₁(0) = ", dg₁dp₁(sol.t[begin], sol.q[begin], sol.p[begin]) * dg₂dq₁(sol.t[begin], sol.q[begin], sol.p[begin]))
println("dg₁dp₂(0) * dg₂dq₂(0) = ", dg₁dp₂(sol.t[begin], sol.q[begin], sol.p[begin]) * dg₂dq₂(sol.t[begin], sol.q[begin], sol.p[begin]))
println("dg₁dp₃(0) * dg₂dq₃(0) = ", dg₁dp₃(sol.t[begin], sol.q[begin], sol.p[begin]) * dg₂dq₃(sol.t[begin], sol.q[begin], sol.p[begin]))
println()

println()
println("dg₁dq₁(0) = ", dg₁dq₁(sol.t[begin], sol.q[begin], sol.p[begin]))
println("dg₁dq₂(0) = ", dg₁dq₂(sol.t[begin], sol.q[begin], sol.p[begin]))
println("dg₁dq₃(0) = ", dg₁dq₃(sol.t[begin], sol.q[begin], sol.p[begin]))
println("dg₁dp₁(0) = ", dg₁dp₁(sol.t[begin], sol.q[begin], sol.p[begin]))
println("dg₁dp₂(0) = ", dg₁dp₂(sol.t[begin], sol.q[begin], sol.p[begin]))
println("dg₁dp₃(0) = ", dg₁dp₃(sol.t[begin], sol.q[begin], sol.p[begin]))
println("dg₂dp₁(0) = ", dg₂dp₁(sol.t[begin], sol.q[begin], sol.p[begin]))
println("dg₂dp₂(0) = ", dg₂dp₂(sol.t[begin], sol.q[begin], sol.p[begin]))
println("dg₂dp₃(0) = ", dg₂dp₃(sol.t[begin], sol.q[begin], sol.p[begin]))
println("dg₂dq₁(0) = ", dg₂dq₁(sol.t[begin], sol.q[begin], sol.p[begin]))
println("dg₂dq₂(0) = ", dg₂dq₂(sol.t[begin], sol.q[begin], sol.p[begin]))
println("dg₂dq₃(0) = ", dg₂dq₃(sol.t[begin], sol.q[begin], sol.p[begin]))
println()
