#
# What do the models actually conserve, and how well?
#
# Every model in this package claims a Hamiltonian structure, and most claim one or more further
# invariants. This tabulates them side by side on the same equilibrium family, which is the check
# that the diagnostics and the equations of motion agree with each other.
#
# For the 3D guiding centre model the interesting quantity is not the energy but the three
# constraints: they vanish along the exact flow, so their magnitude is the drift off the manifold
# on which the momentum is the guiding centre one-form, and keeping it bounded is the whole point
# of the approximately symplectic methods of Li, Zhang & Liu. The constraint the retained pair omits
# is reported separately from the two it keeps, because it drifts 1/bₘ times as much.
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
    @printf("  %-30s %-8s %-14s %-14s %-14s %s\n",
            "equilibrium", "pair", "|ΔH/H|", "|gᵏ| retained", "|gᵏ| omitted", "|Δp_φ/p_φ|")

    # `SolovevSymmetricField` is the one equilibrium with an initial condition that is left out. Its
    # pair is well conditioned — `(g², g³)` is the only regular one there, and λₒ ≈ -228 — but it is
    # the stiffest equilibrium in the package, `‖ϑ‖ ≈ 256`, and the constraint drift grows fast enough
    # that Newton meets a NaN a few hundred steps into the thousand this section runs. See the note in
    # its test block in `test/guiding_center_3d_tests.jl`.
    for (name, M) in (("SolovevIterXpoint",       G3.SolovevIterXpoint),
                      ("TokamakMediumCartesian",  G3.TokamakMediumCartesian),
                      ("TokamakSmallCylindrical", G3.TokamakSmallCylindrical),
                      ("TokamakMediumCylindrical",G3.TokamakMediumCylindrical),
                      ("TokamakSmallToroidal",    G3.TokamakSmallToroidal))
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
        # All three components of v × b vanish along the exact flow, but not equally well along a
        # numerical one: the identity b₁g² - b₂g¹ + b₃g³ = 0 determines the omitted constraint from
        # the retained pair divided by the component bₘ of b the pair leaves out, so the omitted one
        # drifts 1/bₘ times as much. Reporting the two separately is what makes that visible.
        gs = (c.g₁, c.g₂, c.g₃)
        retained = map(v -> only(typeof(v).parameters), M.constraint_pair(M.default_constraints()))
        omitted = only(setdiff(1:3, retained))

        @printf("  %-30s %-8s %-14.2e %-14.2e %-14.2e %s\n", name, M.default_constraints(),
                mx(eerr), maximum(mx(gs[k]) for k in retained), mx(gs[omitted]), pstr)
    end

    println("""
      Every equilibrium that ships an initial condition integrates now that the constraint pair is
      selectable; the cylindrical and toroidal ones used to fail on the hard-coded pair. See
      study_guiding_center_3d_conditioning.jl for the conditioning of each pair, and note that the
      omitted constraint is the one to watch: it is the retained pair's drift amplified by 1/bₘ, which
      is why the equilibria with a small bₘ show a worse third column than second.""")
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
