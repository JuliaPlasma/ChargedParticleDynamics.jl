using GeometricIntegrators
using SimpleSolvers: Options

using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCartesian
using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCartesian: hamiltonian, hamiltonian_u, g₁, g₂, g₃, λₒ, λ₁, λ₂, b₁, b₂, b₃
using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCartesian: gᵏ, dgᵏdqₗ, dgᵏdpₗ, constraint_pair, default_constraints

# using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCylindrical
# using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCylindrical: hamiltonian, hamiltonian_u, g₁, g₂, g₃, λₒ, λ₁, λ₂, b₁, b₂, b₃
# using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCylindrical: gᵏ, dgᵏdqₗ, dgᵏdpₗ, constraint_pair, default_constraints

# using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallToroidal
# using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallToroidal: hamiltonian, hamiltonian_u, g₁, g₂, g₃, λₒ, λ₁, λ₂, b₁, b₂, b₃
# using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallToroidal: gᵏ, dgᵏdqₗ, dgᵏdpₗ, constraint_pair, default_constraints

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
g3 = [g₃(sol.t[i], sol.q[i], sol.p[i]) for i in eachindex(sol.t)]
b1 = [b₁(sol.t[i], sol.q[i]) for i in eachindex(sol.t)]
b2 = [b₂(sol.t[i], sol.q[i]) for i in eachindex(sol.t)]
b3 = [b₃(sol.t[i], sol.q[i]) for i in eachindex(sol.t)]

# The pair the problem above was built with, in the order the multipliers see it. `λ₁` divides
# `{g₂, H}` and `λ₂` divides `-{g₁, H}` by `λₒ = {g₁, g₂}`, where `g₁` and `g₂` mean these two.
const c = constraint_pair(default_constraints())

println()
println("constraint pair = ", default_constraints(), "  ", c)
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

# All three vanish along the exact flow, whichever two the pair retains, so all three are worth
# printing: a drift that cancels in the retained pair still shows up in the third.
println("g¹(0) = ", g1[begin], "   g¹(T) = ", g1[end])
println("g²(0) = ", g2[begin], "   g²(T) = ", g2[end])
println("g³(0) = ", g3[begin], "   g³(T) = ", g3[end])
println()
println("λ₀(0) = ", λ0[begin])
println("λ₀(T) = ", λ0[end])
println()
println("λ₁ g₁ + λ₂ g₂ (0) = ", λ1[begin] * gᵏ(c[1], sol.t[begin], sol.q[begin], sol.p[begin]) +
                                λ2[begin] * gᵏ(c[2], sol.t[begin], sol.q[begin], sol.p[begin]))
println("λ₁ g₁ + λ₂ g₂ (T) = ", λ1[end] * gᵏ(c[1], sol.t[end], sol.q[end], sol.p[end]) +
                                λ2[end] * gᵏ(c[2], sol.t[end], sol.q[end], sol.p[end]))
println()

println()
println("b₁(x₀) = ", b1[begin])
println("b₂(x₀) = ", b2[begin])
println("b₃(x₀) = ", b3[begin])
println()

# The six products that make up λₒ = Σₗ (∂g₁/∂qₗ ∂g₂/∂pₗ - ∂g₁/∂pₗ ∂g₂/∂qₗ). Which of them carry the
# vanishing component of b is what decides whether the pair is usable at this initial condition.
let t₀ = sol.t[begin], q₀ = sol.q[begin], p₀ = sol.p[begin]
    println()
    for i in 1:3
        l = Val(i)
        println("∂g₁/∂q$i * ∂g₂/∂p$i = ", dgᵏdqₗ(c[1], l, t₀, q₀, p₀) * dgᵏdpₗ(c[2], l, t₀, q₀, p₀))
    end
    for i in 1:3
        l = Val(i)
        println("∂g₁/∂p$i * ∂g₂/∂q$i = ", dgᵏdpₗ(c[1], l, t₀, q₀, p₀) * dgᵏdqₗ(c[2], l, t₀, q₀, p₀))
    end
    println()

    println()
    for (label, k) in (("g₁", c[1]), ("g₂", c[2]))
        for i in 1:3
            println("∂$label/∂q$i = ", dgᵏdqₗ(k, Val(i), t₀, q₀, p₀))
        end
        for i in 1:3
            println("∂$label/∂p$i = ", dgᵏdpₗ(k, Val(i), t₀, q₀, p₀))
        end
    end
    println()
end
