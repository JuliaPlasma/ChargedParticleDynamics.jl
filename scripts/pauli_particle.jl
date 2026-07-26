using GeometricIntegrators
using CairoMakie


# Choose problem
import ChargedParticleDynamics.PauliParticle3d.TokamakIterCylindrical as prob
# import ChargedParticleDynamics.PauliParticle3d.TokamakSmallCartesian as prob
# import ChargedParticleDynamics.PauliParticle3d.TokamakSmallCylindrical as prob
# import ChargedParticleDynamics.PauliParticle3d.TokamakSmallToroidal as prob
# import ChargedParticleDynamics.PauliParticle3d.SolovevIter as prob
# import ChargedParticleDynamics.PauliParticle3d.SolovevIterXpoint as prob

# Choose method (see https://github.com/JuliaGNI/GeometricIntegrators.jl/blob/main/src/integrators/method_list.jl)
method = Gauss(1)
# method = Gauss(2)

# The following only work with pode but not iode
# method = SymplecticEulerA()
# method = SymplecticEulerB()
# method = LobattoIIIAIIIB(2) # leapfrog
# method = LobattoIIIBIIIA(2) # leapfrog

# Solver options, splatted into `integrate` as keywords.
#
# `f_abstol` is an absolute bound on the residual and has to stay above its round-off floor. The
# Pauli residual carries the momentum equation `p = ϑ(q)`, whose scale is set by the one-form, so
# that floor is `‖ϑ‖ eps` — of order 1E-14 for the ITER-scale equilibria. Asking for 1E-14, as this
# script used to, leaves the solver no reachable criterion and it runs to its iteration limit on
# most steps. `f_reltol` is deliberately left at the `SimpleSolvers` default: its relative term is
# what lets a large-magnitude solve converge at all. See `test/guiding_center_3d_tests.jl`.
options = (f_abstol = 1E-12, max_iterations = 50, warn_iterations = 50)

# Create partitioned ODE (Hamilton's equations) and implicit ODE (Euler-Lagrange equations)
pode = prob.podeproblem()
iode = prob.iodeproblem()

# You may want to adapt the keyword arguments timespan = (t₀,t₁) and timestep = Δt, e.g.
# pode = prob.podeproblem(timespan = (0.0, 1E4), timestep = 0.1)
# iode = prob.iodeproblem(timespan = (0.0, 1E4), timestep = 0.1)

# Integrate pode and iode
psol = integrate(pode, method; options...)
isol = integrate(iode, method; options...)

# Plot solution and energy error
function plot_solution(prob, ode, sol, prefix)
    # Compute Cartesian coordinates from solution
    R = prob.R.(sol.t[:], sol.q[:,1], sol.q[:,2], sol.q[:,3])
    X = prob.X.(sol.t[:], sol.q[:,1], sol.q[:,2], sol.q[:,3])
    Y = prob.Y.(sol.t[:], sol.q[:,1], sol.q[:,2], sol.q[:,3])
    Z = prob.Z.(sol.t[:], sol.q[:,1], sol.q[:,2], sol.q[:,3])

    # Compute energy
    hamiltonian = (t,q,p) -> prob.hamiltonian(t, q, p, ode.parameters)
    energy = hamiltonian.(sol.t, sol.q, sol.p)
    energy_error = (energy .- energy[begin]) / energy[begin]

    # Plot trajectory projected on poloidal plane
    fig = Figure(backgroundcolor = :transparent, size = (600, 600), fontsize = 18)
    ax = Axis(
        fig[1, 1],
        aspect = 1,
        xlabel = "R",
        ylabel = "Z",
    )
    lines!(ax, R, Z; linewidth = 3)
    # scatter!(ax, R, Z)
    save(prefix * "_trajectory.pdf", fig)

    # Plot trajectory in 3D
    fig = Figure(backgroundcolor = :transparent, size = (800, 800), fontsize = 18)
    ax = Axis3(
        fig[1, 1],
        xlabel = "X",
        ylabel = "Y",
        zlabel = "Z",
    )
    lines!(ax, X, Y, Z; linewidth = 3)
    # scatter!(ax, X, Y, Z)
    save(prefix * "_trajectory_3d.pdf", fig)

    # Plot energy error
    fig = Figure(backgroundcolor = :transparent, size = (800, 400), fontsize = 18)
    ax = Axis(
        fig[1, 1],
        xlabel = "t",
        ylabel = "(H(t) - H(0)) / H(0)",
    )
    xlims!(ax, (sol.t[begin], sol.t[end]))
    lines!(ax, sol.t[:], energy_error; linewidth = 3)
    save(prefix * "_energy_error.pdf", fig)
end

plot_solution(prob, pode, psol, "pauli_particle_hamiltonian_$(string(nameof(typeof(method))))")
plot_solution(prob, iode, isol, "pauli_particle_lagrangian_$(string(nameof(typeof(method))))")
