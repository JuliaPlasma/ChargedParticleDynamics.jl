using GeometricIntegrators
using CairoMakie
using SimpleSolvers: Options
using Test


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

# Solver options
options = Options(x_reltol = 1E-14, f_abstol = 1E-14, f_reltol = 1E-14)

# Create partitioned ODE (Hamilton's equations) and implicit ODE (Euler-Lagrange equations)
pode = prob.pauli_particle_3d_pode()
iode = prob.pauli_particle_3d_iode()

# You may want to adapt the keyword arguments tspan = (t₀,t₁) and tstep = Δt, e.g.
# pode = prob.pauli_particle_3d_pode(tspan = (0.0, 1E4), tstep = 0.1)
# iode = prob.pauli_particle_3d_iode(tspan = (0.0, 1E4), tstep = 0.1)

# Integrate pode and iode
psol = integrate(pode, method; options = options)
isol = integrate(iode, method; options = options)

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
