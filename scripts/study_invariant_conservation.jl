#
# What do the models actually conserve, and how well?
#
# Every model in this package claims a Hamiltonian structure, and most claim one or more further
# invariants. This tabulates them side by side on the same equilibrium family, which is the check
# that the diagnostics and the equations of motion agree with each other.
#
# For the 3D guiding centre model the interesting quantity is not the energy but the two
# constraints: they vanish along the exact flow, so their magnitude is the drift off the manifold
# on which the momentum is the guiding centre one-form, and keeping it bounded is the whole point
# of the approximately symplectic methods of Li, Zhang & Liu.
#
# See docs/src/findings.md, "Conservation across the model families".
#
# Run with:  julia --project=scripts scripts/study_invariant_conservation.jl
#

using ChargedParticleDynamics
using GeometricIntegrators
using Printf

mx(ds) = maximum(abs(ds[i]) for i in eachindex(ds))

const G3 = ChargedParticleDynamics.GuidingCenter3d
const G4 = ChargedParticleDynamics.GuidingCenter4d
const GK = ChargedParticleDynamics.GyroKinetics4d

function guiding_centre_4d()
    println("4D guiding centre — ODEProblem, Gauss(2)\n")
    @printf("  %-30s %-14s %s\n", "equilibrium", "|ΔH/H|", "|Δp_φ/p_φ|")

    for (name, M, kw) in (("TokamakSmallCylindrical", G4.TokamakSmallCylindrical, (timestep = 10.0, timespan = (0.0, 1E3))),
                          ("TokamakSmallToroidal",    G4.TokamakSmallToroidal,    (timestep = 10.0, timespan = (0.0, 1E3))),
                          ("TokamakSmallCartesian",   G4.TokamakSmallCartesian,   (timestep = 10.0, timespan = (0.0, 1E3))),
                          ("SolovevIterXpoint",       G4.SolovevIterXpoint,       (timestep = 1.0,  timespan = (0.0, 1E2))))
        sol = integrate(M.odeproblem(M.initial_conditions_barely_passing(); kw...), Gauss(2))
        _, eerr = M.compute_energy_error(sol)
        _, perr = M.compute_toroidal_momentum_error(sol)
        @printf("  %-30s %-14.2e %.2e\n", name, mx(eerr), mx(perr))
    end
end

function guiding_centre_3d()
    println("\n3D guiding centre — HODEProblem, PartitionedGauss(2)\n")
    @printf("  %-30s %-14s %-14s %-14s %s\n", "equilibrium", "|ΔH/H|", "|g₁|", "|g₂|", "|Δp_φ/p_φ|")

    for (name, M) in (("SolovevIterXpoint",      G3.SolovevIterXpoint),
                      ("TokamakMediumCartesian", G3.TokamakMediumCartesian))
        prob = M.hodeproblem(M.initial_conditions_barely_passing(); timestep = 0.1, timespan = (0.0, 1E2))
        sol  = integrate(prob, PartitionedGauss(2); initialguess = MidpointExtrapolation(5))
        _, eerr = M.compute_energy_error(sol)
        c = M.compute_constraints(sol)
        pstr = if isdefined(M, :toroidal_momentum)
            _, perr = M.compute_toroidal_momentum_error(sol)
            @sprintf("%.2e", mx(perr))
        else
            "—"
        end
        @printf("  %-30s %-14.2e %-14.2e %-14.2e %s\n", name, mx(eerr), mx(c.g₁), mx(c.g₂), pstr)
    end

    println("""
      Only the equilibria that integrate at all are listed. The cylindrical and toroidal ones fail
      on the singular constraint pair; see study_guiding_center_3d_conditioning.jl.""")
end

function gyrokinetics()
    println("\nGyrokinetic guiding centre — SolovevIterXpoint, rescaled time\n")
    @printf("  %-30s %s\n", "method", "|ΔH/H|")

    M = GK.GuidingCenter4dSolovevIterXpoint
    q₀, params = M.initial_conditions_deeply_passing()

    sol = integrate(M.odeproblem(q₀; parameters = params), Gauss(2))
    _, e = M.compute_energy_error(sol)
    @printf("  %-30s %.2e\n", "Gauss(2), full field", mx(e))

    sol = integrate(M.sodeproblem(q₀; parameters = params),
                    Composition(Tuple(Gauss(1) for _ in 1:6), Strang()))
    _, e = M.compute_energy_error(sol)
    @printf("  %-30s %.2e\n", "Strang split, volume pres.", mx(e))

    println("""
      The splitting is volume preserving rather than symplectic, so its energy is not conserved by
      construction — it is bounded here because the composition is symmetric and the orbit is
      short. Energy behaviour over long runs is the open question for this scheme.""")
end

guiding_centre_4d()
guiding_centre_3d()
gyrokinetics()
