#
# Guiding centre dynamics in 4D: integrate one equilibrium and plot the trajectory and the energy
# error. The companion of `pauli_particle.jl`, which does the same for the Pauli particle.
#
# This file used to be a `Weave` driver for `guiding_center_dynamics_{splitting,vprk}.jmd`. Those
# notebooks have been removed, so it is written out here as a plain script instead.
#

using GeometricIntegrators
using CairoMakie


# Choose problem
import ChargedParticleDynamics.GuidingCenter4d.TokamakIterCylindrical as prob
# import ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCartesian as prob
# import ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCylindrical as prob
# import ChargedParticleDynamics.GuidingCenter4d.TokamakSmallToroidal as prob
# import ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCartesian as prob
# import ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCylindrical as prob
# import ChargedParticleDynamics.GuidingCenter4d.SolovevIter as prob
# import ChargedParticleDynamics.GuidingCenter4d.SolovevIterXpoint as prob

# The gyrokinetic model provides the same equilibria and the same constructors. Note that its
# independent variable is the rescaled time, with dt = B*∥ ds, so its default step is the one above
# divided by B*∥ — do not carry a `timestep` across from one model to the other.
# import ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dTokamakIterCylindrical as prob

# Choose initial condition
ics = prob.initial_conditions_deeply_passing()
# ics = prob.initial_conditions_barely_passing()
# ics = prob.initial_conditions_deeply_trapped()
# ics = prob.initial_conditions_barely_trapped()

# Choose method (see https://github.com/JuliaGNI/GeometricIntegrators.jl/blob/main/src/integrators/method_list.jl)
method = Gauss(2)

# Solver options, splatted into `integrate` as keywords. See the comment in `pauli_particle.jl`
# for why `f_abstol` has to stay above the round-off floor of the residual.
options = (f_abstol = 1E-12, max_iterations = 50, warn_iterations = 50)

# The explicit form, ẋ = v(x), and the implicit one built from the one-form ϑ = A + u b. The two
# describe the same system; the second is the one the variational integrators apply to.
ode = prob.odeproblem(ics)
iode = prob.iodeproblem(ics)

# You may want to adapt the keyword arguments timespan = (t₀,t₁) and timestep = Δt, e.g.
# ode  = prob.odeproblem(ics; timespan = (0.0, 1E4), timestep = 1.0)
# iode = prob.iodeproblem(ics; timespan = (0.0, 1E4), timestep = 1.0)

osol = integrate(ode, method; options...)
isol = integrate(iode, VPRKGauss(2); options...)


# Plot solution and energy error
function plot_solution(prob, problem, sol, prefix)
    # Compute cartesian coordinates from the solution. The first three components of the state are
    # the guiding centre position in the equilibrium's own chart; the fourth is the parallel
    # velocity.
    R = prob.R.(sol.t[:], sol.q[:, 1], sol.q[:, 2], sol.q[:, 3])
    X = prob.X.(sol.t[:], sol.q[:, 1], sol.q[:, 2], sol.q[:, 3])
    Y = prob.Y.(sol.t[:], sol.q[:, 1], sol.q[:, 2], sol.q[:, 3])
    Z = prob.Z.(sol.t[:], sol.q[:, 1], sol.q[:, 2], sol.q[:, 3])

    # Compute energy. The guiding centre Hamiltonian is a function of the state alone — unlike the
    # Pauli particle's, which takes (q, p) — so it is evaluated on `sol.q` only.
    energy = [prob.hamiltonian(sol.t[i], sol.q[i], parameters(problem)) for i in eachindex(sol.t)]
    energy_error = (energy .- energy[begin]) ./ energy[begin]

    # Plot trajectory projected on poloidal plane
    fig = Figure(backgroundcolor = :transparent, size = (600, 600), fontsize = 18)
    ax = Axis(
        fig[1, 1],
        aspect = 1,
        xlabel = "R",
        ylabel = "Z",
    )
    lines!(ax, R, Z; linewidth = 3)
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

plot_solution(prob, ode, osol, "guiding_center_ode_$(string(nameof(typeof(method))))")
plot_solution(prob, iode, isol, "guiding_center_iode_VPRKGauss")
