using GeometricIntegrators

using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallCartesian
using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallCartesian: hamiltonian, hamiltonian_u, g₁, g₂, g₃, λₒ, λ₁, λ₂, b₁, b₂, b₃

# See `guiding_center_3d_constraint.jl` for why `f_abstol` cannot be pushed to the round-off floor:
# `2eps()` leaves the solver no reachable stopping criterion on these models.
const options = (f_abstol=1E-12, max_iterations=50, warn_iterations=50)

equ = hodeproblem(initial_conditions_default())
# equ = hodeproblem(initial_conditions_default(); timestep=200.0, timespan=(0.0, 2E5))
# equ = hodeproblem(initial_conditions_default(); timestep=400.0, timespan=(0.0, 2E5))


# sol = integrate(equ, PartitionedGauss(1); options...)
sol = integrate(equ, PartitionedGauss(1); initialguess=MidpointExtrapolation(2), options...)
# sol = integrate(equ, PartitionedGauss(1); initialguess=NoInitialGuess(), options...)
# sol = integrate(equ, PartitionedGauss(2); initialguess=NoInitialGuess(), options...)
# sol = integrate(equ, PartitionedGauss(2); initialguess=MidpointExtrapolation(2), options...)
# sol = integrate(equ, IRK3(); options...)


h = [hamiltonian(sol.t[i], sol.q[i], sol.p[i], parameters(equ)) for i in eachindex(sol.t)]
h0 = h[begin]
hu = [hamiltonian_u(sol.t[i], sol.q[i], sol.p[i], parameters(equ)) for i in eachindex(sol.t)]
λ0 = [λₒ(sol.t[i], sol.q[i], sol.p[i]) for i in eachindex(sol.t)]
λ1 = [λ₁(sol.t[i], sol.q[i], sol.p[i], parameters(equ)) for i in eachindex(sol.t)]
λ2 = [λ₂(sol.t[i], sol.q[i], sol.p[i], parameters(equ)) for i in eachindex(sol.t)]
g1 = [g₁(sol.t[i], sol.q[i], sol.p[i]) for i in eachindex(sol.t)]
g2 = [g₂(sol.t[i], sol.q[i], sol.p[i]) for i in eachindex(sol.t)]
g3 = [g₃(sol.t[i], sol.q[i], sol.p[i]) for i in eachindex(sol.t)]
b1 = [b₁(sol.t[i], sol.q[i]) for i in eachindex(sol.t)]
b2 = [b₂(sol.t[i], sol.q[i]) for i in eachindex(sol.t)]
b3 = [b₃(sol.t[i], sol.q[i]) for i in eachindex(sol.t)]
R = [sqrt(sol.q[i, 1]^2 + sol.q[i, 2]^2) for i in eachindex(sol.t)]
Z = sol.q[:, 3]

println()
println("q(0) = ", sol.q[begin])
println("q(T) = ", sol.q[end])
println()
println("p(0) = ", sol.p[begin])
println("p(T) = ", sol.p[end])
println()
println("H(0) = ", h[begin])
println("H(T) = ", h[end])
println()
println("Hᵤ(0) = ", hu[begin])
println("Hᵤ(T) = ", hu[end])
println()
println("g₁(0) = ", g1[begin])
println("g₁(T) = ", g1[end])
println()
println("g₂(0) = ", g2[begin])
println("g₂(T) = ", g2[end])

println("g₃(0) = ", g3[begin])
println("g₃(T) = ", g3[end])
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


using CairoMakie

f = Figure(size=(1000, 800))

axsol = Axis(f[1, 1], xlabel="R", ylabel="Z")
axham = Axis(f[1, 2], xlabel="t", ylabel="[H(t) - H(0)] / H(0)")
axg1 = Axis(f[2, 1], xlabel="t", ylabel="g₁")
axg2 = Axis(f[2, 2], xlabel="t", ylabel="g₂")

scatter!(axsol, R, Z)
plot!(axham, sol.t, (h .- h0) ./ h0)
plot!(axg1, sol.t, g1)
plot!(axg2, sol.t, g2)

save("tokamak-3d.png", f)
