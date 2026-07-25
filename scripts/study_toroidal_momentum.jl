#
# Which expression is the conserved toroidal momentum?
#
# Every guiding centre equilibrium with a symmetry direction used to define
#
#     toroidal_momentum(t, q) = R(t, q) * ϑ₃(t, q)
#
# This script was written to settle whether the factor of R belongs there. It does not: ϑ₃ is
# already the covariant φ-component of the one-form and is the canonical momentum on its own.
#
# In cartesian coordinates the question does not even arise in that form — the third coordinate is
# z, not an angle — and the conserved quantity is the generator of rotation about the z-axis.
#
# See docs/src/findings.md, "The toroidal momentum of the guiding centre models".
#
# Run with:  julia --project=scripts scripts/study_toroidal_momentum.jl
#

using ChargedParticleDynamics
using GeometricIntegrators
using Printf

const GC = ChargedParticleDynamics.GuidingCenter4d

"Relative variation of a data series over the whole solution."
relvar(x) = (maximum(x) - minimum(x)) / max(abs(x[1]), eps())

function candidates(M, sol)
    ϑ₃   = [M.ϑ₃(sol.t[i], sol.q[i]) for i in eachindex(sol.t)]
    Rϑ₃  = [M.R(sol.t[i], sol.q[i]) * M.ϑ₃(sol.t[i], sol.q[i]) for i in eachindex(sol.t)]
    ang  = [sol.q[i][1] * M.ϑ₂(sol.t[i], sol.q[i]) -
            sol.q[i][2] * M.ϑ₁(sol.t[i], sol.q[i]) for i in eachindex(sol.t)]
    (ϑ₃ = relvar(ϑ₃), Rϑ₃ = relvar(Rϑ₃), angular = relvar(ang))
end

function compare()
    println("Relative variation over the orbit of each candidate, Gauss(2)\n")
    @printf("  %-34s %-10s %-12s %-12s %-12s\n", "equilibrium", "coords", "ϑ₃", "R ϑ₃", "x ϑ₂ - y ϑ₁")

    cases = (("TokamakSmallCylindrical",  GC.TokamakSmallCylindrical,  "cylindrical", (tstep = 10.0, tspan = (0.0, 1E3))),
             ("TokamakSmallToroidal",     GC.TokamakSmallToroidal,     "toroidal",    (tstep = 10.0, tspan = (0.0, 1E3))),
             ("SolovevIterXpoint",        GC.SolovevIterXpoint,        "cylindrical", (tstep = 1.0,  tspan = (0.0, 1E2))),
             ("TokamakSmallCartesian",    GC.TokamakSmallCartesian,    "cartesian",   (tstep = 10.0, tspan = (0.0, 1E3))),
             ("TokamakMediumCartesian",   GC.TokamakMediumCartesian,   "cartesian",   (tstep = 1.0,  tspan = (0.0, 1E3))))

    for (name, M, coords, kw) in cases
        prob = M.guiding_center_4d_ode(M.initial_conditions_barely_passing()...; kw...)
        sol  = integrate(prob, Gauss(2))
        c    = candidates(M, sol)
        @printf("  %-34s %-10s %-12.2e %-12.2e %-12.2e\n", name, coords, c.ϑ₃, c.Rϑ₃, c.angular)
    end
end

"""
The residual of the cartesian candidate is discretisation error, not a defect of the quantity:
halving the step reduces it by a factor of sixteen, as befits the fourth-order Gauss method.
"""
function convergence()
    println("\nCartesian candidate x ϑ₂ - y ϑ₁ under step refinement (TokamakMediumCartesian)\n")
    @printf("  %-10s %-12s %s\n", "Δt", "rel. var.", "ratio")

    M = GC.TokamakMediumCartesian
    previous = NaN
    for Δt in (1.0, 0.5, 0.25, 0.125)
        prob = M.guiding_center_4d_ode(M.initial_conditions_barely_passing()...; tstep = Δt, tspan = (0.0, 1E3))
        sol  = integrate(prob, Gauss(2))
        v    = candidates(M, sol).angular
        @printf("  %-10.3g %-12.2e %s\n", Δt, v, isnan(previous) ? "" : @sprintf("%.1f", previous / v))
        previous = v
    end
end

compare()
convergence()
