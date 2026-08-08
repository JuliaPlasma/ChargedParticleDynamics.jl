#
# Do the three guiding-centre-related model families integrate the same orbit from the same
# (x, u, μ)?
#
# The question has a different answer for each pair, and the difference is the point:
#
#   * `GuidingCenter3d` and `GuidingCenter4d` are the *same model* in different variables. The 4D
#     form carries the state (x, u) and the 3D form the position and its conjugate momentum
#     p = ϑ(x, u) = A + u b, with the parallel velocity eliminated by the constraints v × b = 0. They
#     must therefore agree to the accuracy of the nonlinear solves, and section 1 measures that they
#     do. This is an equivalence check, and a disagreement would be a bug in one of them.
#
#   * `PauliParticle3d` is a *different model*: a regular Lagrangian on the full six-dimensional
#     phasespace, differing from the charged particle only by the term μ|B|. Guiding centre dynamics
#     is its slow manifold — see the "Slow Manifolds" section of `docs/src/pauli_particle_3d.md` — so
#     the two agree only to the order at which the initial condition is placed on that manifold.
#     Section 2 measures how well, and section 3 identifies what limits it.
#
# Section 3 is the substantive finding. Placing the Pauli particle at v = u b⃗ puts it on the
# *lowest-order* slow manifold: purely parallel, with no drift. The guiding centre moves with the
# ∇B and curvature drifts as well, so the two initial velocities differ by the drift velocity, and
# the Pauli particle carries that difference as a gyration for the whole orbit. Supplying the
# guiding centre velocity itself as the Pauli initial velocity removes it, and buys up to three
# orders of magnitude.
#
# See docs/src/findings.md, "The three families on one initial condition".
#
# Run with:  julia --project=scripts scripts/study_model_agreement.jl
#

using ChargedParticleDynamics
using GeometricIntegrators
using LinearAlgebra
using Printf

const G3 = ChargedParticleDynamics.GuidingCenter3d
const G4 = ChargedParticleDynamics.GuidingCenter4d
const P3 = ChargedParticleDynamics.PauliParticle3d

# The tolerance the test suite settled on; see `scripts/study_solver_tolerances.jl`.
const OPTS = (f_abstol = 1E-12, max_iterations = 50)

# A thousand steps at a step every family resolves. The guiding centre modules declare steps as large
# as 500 for the small tokamak and 1.0 for the ITER-size ones, but the Pauli model has to resolve a
# gyration the guiding centre model has averaged away, so the common step is the Pauli one.
const CASES = (("TokamakSmallCartesian",   10.0,  1000),
               ("TokamakSmallCylindrical", 10.0,  1000),
               ("TokamakSmallToroidal",    10.0,  1000),
               ("TokamakIterCylindrical",  0.01,  1000),
               ("SolovevIter",             0.01,  1000),
               ("SolovevIterXpoint",       0.01,  1000))

const ICS = (:initial_conditions_barely_passing, :initial_conditions_barely_trapped,
             :initial_conditions_deeply_passing, :initial_conditions_deeply_trapped)

label(icn) = replace(string(icn), "initial_conditions_" => "")

# `periodic = false` throughout. The guiding centre problems wrap their angular coordinate into the
# chart's range and the Pauli problems do not, so with wrapping left on, a passing orbit's ϕ is
# compared against the same ϕ plus a multiple of 2π — which reads as a complete disagreement and was
# how these two rows first looked.
function g4(M, x, u, μ, span, step)
    sol = integrate(M.odeproblem([x..., u]; parameters = (μ = μ,), timespan = span, timestep = step,
                                 periodic = false), Gauss(2); OPTS...)
    (sol.q[end][1:3], sol.q[end][4])
end

function g3(M, x, u, μ, span, step)
    sol = integrate(M.hodeproblem([x..., u]; parameters = (μ = μ,), timespan = span, timestep = step,
                                  periodic = false), PartitionedGauss(2); OPTS...)
    (collect(sol.q[end]), M.u(span[end], sol.q[end], sol.p[end]))
end

"""
The Pauli orbit from the initial velocity `v₀`, which is what the sections below vary: `u b⃗` for the
lowest-order slow manifold, the guiding centre velocity for the drift-corrected one.
"""
function pauli(M, x, v₀, μ, span, step)
    sol = integrate(M.hodeproblem(x, v₀, μ; timespan = span, timestep = step),
                    PartitionedGauss(2); OPTS...)
    q, p = sol.q[end], sol.p[end]
    (collect(q), M.v(span[end], q, p)' * M.b(span[end], q))
end

reldiff(a, b) = norm(a - b) / max(norm(b), 1E-30)

"The gyroradius over the field's own gradient scale length, ρ/L with ρ = √(2μ/B) and L = B/|∇B|."
function ρ_over_L(M, x, μ)
    B₀ = M.B(0.0, x)
    ∇B = norm([M.dBdx₁(0.0, x), M.dBdx₂(0.0, x), M.dBdx₃(0.0, x)])
    sqrt(2μ / B₀) * ∇B / B₀
end

"Every case as `(name, condition, x, u, μ, span, step)`, taking the state from the 4D module."
function cases()
    out = []
    for (name, step, nsteps) in CASES, icn in ICS
        M4 = getfield(G4, Symbol(name))
        isdefined(M4, icn) || continue
        ic = getfield(M4, icn)()
        push!(out, (name, icn, ic.q[1:3], ic.q[4], ic.params.μ, (0.0, nsteps * step), step))
    end
    out
end


# ---------------------------------------------------------------------------------------------
# 1. The two guiding centre formulations are the same model
# ---------------------------------------------------------------------------------------------

function equivalence()
    println("""
    1. GuidingCenter3d against GuidingCenter4d
    =========================================

    The same model in different variables, so this is an equivalence check rather than a comparison:
    `Δx` and `Δu` are the relative differences of the final position and parallel velocity, with the
    4D form as the reference. `PartitionedGauss(2)` on the 3D `hodeproblem` against `Gauss(2)` on the
    4D `odeproblem`, a thousand steps.
    """)

    @printf("  %-24s %-16s %12s %12s\n", "equilibrium", "condition", "Δx", "Δu")

    for (name, icn, x, u, μ, span, step) in cases()
        r4 = g4(getfield(G4, Symbol(name)), x, u, μ, span, step)
        r3 = g3(getfield(G3, Symbol(name)), x, u, μ, span, step)
        @printf("  %-24s %-16s %12.3e %12.3e\n", name, label(icn),
                reldiff(r3[1], r4[1]), abs(r3[2] - r4[2]) / max(abs(r4[2]), 1E-30))
        flush(stdout)
    end

    println("""

    The curvilinear charts agree at 1E-13 and below, which is the accuracy of the nonlinear solves
    rather than anything about the models. The cartesian chart is three to five orders looser, and
    that is a property of the chart and not of the equivalence: there the toroidal transit appears as
    a rotation of (x, y) that has to be resolved, where in the other charts it is a coordinate that
    simply winds. The same asymmetry sets that chart's time step; see the note in
    `pauli_particle_3d/tokamak_small_cartesian.jl`.""")
end


# ---------------------------------------------------------------------------------------------
# 2. The Pauli particle on the lowest-order slow manifold
# ---------------------------------------------------------------------------------------------

function slowmanifold()
    println("""


    2. PauliParticle3d against GuidingCenter4d, at v = u b⃗
    ======================================================

    A different model, agreeing only to the order at which its initial condition sits on the slow
    manifold. `v = u b⃗` is purely parallel, which is that manifold to lowest order. `ρ/L` is the
    gyroradius over the field's gradient scale length, the parameter the guiding centre expansion is
    an expansion in.
    """)

    @printf("  %-24s %-16s %12s %12s %12s\n", "equilibrium", "condition", "Δx", "Δu", "ρ/L")

    for (name, icn, x, u, μ, span, step) in cases()
        M4 = getfield(G4, Symbol(name))
        r4 = g4(M4, x, u, μ, span, step)
        rp = pauli(getfield(P3, Symbol(name)), x, u * getfield(P3, Symbol(name)).b⃗(0.0, x), μ, span, step)
        @printf("  %-24s %-16s %12.3e %12.3e %12.3e\n", name, label(icn),
                reldiff(rp[1], r4[1]), abs(rp[2] - r4[2]) / max(abs(r4[2]), 1E-30), ρ_over_L(M4, x, μ))
        flush(stdout)
    end

    println("""

    Across equilibria the agreement is ordered by ρ/L, as the expansion says it should be: 1.7E-6 to
    1.6E-4 for the small tokamak at ρ/L ≈ 2E-3, 2E-4 for the ITER-size tokamak at 2E-2, and 2E-3 to
    1.8E-2 for the ITER Solov'ev at 2E-1.

    Within one equilibrium it is *not* a power of either parameter. Holding the position and span
    fixed on `SolovevIterXpoint` and taking μ down by four decades leaves the difference at 1.7E-2,
    and taking u down by four decades leaves it at 1.5E-3; both saturate rather than converge. That
    floor is what section 3 identifies.""")
end


# ---------------------------------------------------------------------------------------------
# 3. What the floor is
# ---------------------------------------------------------------------------------------------

function driftcorrection()
    println("""


    3. The floor is the initial condition, not the models
    ====================================================

    The guiding centre does not move along `b`: it moves along `b` *plus* the ∇B and curvature
    drifts. So `v = u b⃗` and the guiding centre velocity differ by the drift velocity, the Pauli
    particle starts that far off the manifold, and it carries the difference as a gyration of radius
    ~|v_drift|/B for the whole orbit — an offset that does not shrink when μ or u does, because the
    drift shrinks with them.

    Supplying the guiding centre velocity itself as the Pauli initial velocity removes it. Nothing
    else changes: same position, same μ, same step, same span, same method.
    """)

    @printf("  %-24s %-16s %12s %12s %8s\n", "equilibrium", "condition",
            "v = u b⃗", "v = v_gc", "gain")

    for (name, icn, x, u, μ, span, step) in cases()
        M4, Mp = getfield(G4, Symbol(name)), getfield(P3, Symbol(name))
        ref = g4(M4, x, u, μ, span, step)[1]

        vgc = zeros(4)
        M4.guiding_center_4d_v(vgc, 0.0, [x..., u], (μ = μ,))

        da = reldiff(pauli(Mp, x, u * Mp.b⃗(0.0, x), μ, span, step)[1], ref)
        db = reldiff(pauli(Mp, x, vgc[1:3], μ, span, step)[1], ref)

        @printf("  %-24s %-16s %12.3e %12.3e %8.1f\n", name, label(icn), da, db, da / db)
        flush(stdout)
    end

    println("""

    **The spread across equilibria collapses from four decades to two.** In section 2 the difference
    ran from 1.7E-6 on the small tokamak to 1.8E-2 on the ITER Solov'ev, tracking ρ/L; here every
    equilibrium sits between 1.9E-6 and 1.8E-4, and the ITER Solov'ev — three orders worse than the
    small tokamak before — is no longer distinguishable from it. That is the check that the two models
    really do share a slow manifold: most of what separated them was where the initial condition was
    put on it, not the dynamics.

    The gain is ordered by ρ/L, which is the expected shape, but it is **not** a gain everywhere:

      ITER Solov'ev,  ρ/L = 2E-1   650-1080x better
      ITER tokamak,   ρ/L = 2E-2    14-37x better
      small tokamak,  ρ/L = 2E-3    0.1-0.5x, i.e. 2-10x *worse* on three of every four conditions

    The sign change is not noise — it is consistent across charts and conditions — and it says the
    correction is a partial one. Supplying `v_gc` fixes the O(v_drift) mismatch in the initial
    *velocity*, but the Pauli particle's state is the *particle* position, which differs from the
    guiding centre position by the gyroradius vector, and that O(ρ) mismatch is left in place. Where
    the velocity term dominates the correction wins by three orders; where ρ/L is small enough that
    the residual is already down at the position term, adding one without the other moves the total
    the wrong way. Both corrections together would be the first-order slow manifold proper.

    `PauliParticle3d.initial_conditions(x, u, μ)` deliberately does *not* do any of this. It would
    make the module depend on `GuidingCenter4d` for a drift formula, and the lowest-order condition is
    the self-contained and conventional one. What this section is for is to say what the residual in
    section 2 *is*, so that it is not read as a disagreement between the models. See `TODO.md`.""")
end


function main()
    equivalence()
    slowmanifold()
    driftcorrection()
end

main()
