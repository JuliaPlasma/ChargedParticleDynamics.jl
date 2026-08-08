#
# The three guiding-centre-related families on one initial condition.
#
# `GuidingCenter3d` and `GuidingCenter4d` are the same model in different variables, so they must
# agree to the accuracy of the nonlinear solves. `PauliParticle3d` is a different model whose slow
# manifold is the guiding centre one, so it agrees to the order at which its initial condition is
# placed on that manifold. Both statements are asserted below, and the second is the reason the
# families' `initial_conditions_*` are worth keeping aligned at all.
#
# See `scripts/study_model_agreement.jl` and docs/src/findings.md, "The three families on one
# initial condition", for the measurements the thresholds here are set from.
#

using SafeTestsets


module ModelAgreementUtils

using ChargedParticleDynamics
using GeometricIntegrators
using LinearAlgebra

export G3, G4, P3, final_g3, final_g4, final_pauli, reldiff, ICS, sharedmodules

const G3 = ChargedParticleDynamics.GuidingCenter3d
const G4 = ChargedParticleDynamics.GuidingCenter4d
const P3 = ChargedParticleDynamics.PauliParticle3d

# The tolerance the rest of the suite uses; see `scripts/study_solver_tolerances.jl`.
const OPTS = (f_abstol = 1E-12, max_iterations = 50)

const ICS = (:initial_conditions_barely_passing, :initial_conditions_barely_trapped,
             :initial_conditions_deeply_passing, :initial_conditions_deeply_trapped,
             :initial_conditions_trapped, :initial_conditions_pauli)

# `periodic = false`: the guiding centre problems wrap their angular coordinate into the chart's
# range and the Pauli problems do not, so a passing orbit's ϕ would otherwise be compared against
# the same ϕ plus a multiple of 2π.
function final_g4(M, x, u, μ, span, step)
    sol = integrate(M.odeproblem([x..., u]; parameters = (μ = μ,), timespan = span,
                                 timestep = step, periodic = false), Gauss(2); OPTS...)
    (sol.q[end][1:3], sol.q[end][4])
end

function final_g3(M, x, u, μ, span, step)
    sol = integrate(M.hodeproblem([x..., u]; parameters = (μ = μ,), timespan = span,
                                  timestep = step, periodic = false), PartitionedGauss(2); OPTS...)
    (collect(sol.q[end]), M.u(span[end], sol.q[end], sol.p[end]))
end

function final_pauli(M, x, v₀, μ, span, step)
    sol = integrate(M.hodeproblem(x, v₀, μ; timespan = span, timestep = step),
                    PartitionedGauss(2); OPTS...)
    (collect(sol.q[end]), nothing)
end

reldiff(a, b) = norm(a - b) / max(norm(b), 1E-30)

"Every equilibrium name that all three families define, in a fixed order."
sharedmodules() = sort([string(n) for n in names(P3, all = true)
                        if isa(getfield(P3, n), Module) && n !== nameof(P3) &&
                           isdefined(G3, n) && isdefined(G4, n)])

end


@safetestset "Model families: the declared initial conditions are the same across all three          " begin

    using ChargedParticleDynamics
    using ..ModelAgreementUtils

    # This is the guard on the alignment rather than on any dynamics, and it is cheap. Before it, the
    # Pauli family declared `initial_conditions_*` for three of its eight equilibria, was missing
    # `barely_passing` on one of those three, and `GuidingCenter3d` carried a `y = 0.1` offset on
    # three modules that its 4D counterpart did not — so nothing compared the families at all.
    #
    # All three offsets are gone. `TokamakMediumCartesian`'s cost the most to remove and is not shared
    # with the Pauli family: `y = 0` puts it on the `b₁ = b_x = 0` midplane, where `:g31` is singular,
    # so its `default_constraints` moved to `:g23` — the better conditioned pair in any case — and its
    # `hodeproblem_canonical` now needs a tenth of the step. See `guiding_center_3d_tests.jl`.
    #
    # Guarded against becoming vacuous: the loops below skip an equilibrium that no longer declares a
    # condition in all three families, so an empty module list would pass silently.
    @test !isempty(sharedmodules())

    for name in sharedmodules()
        E3, E4, Ep = getfield(G3, Symbol(name)), getfield(G4, Symbol(name)), getfield(P3, Symbol(name))

        for i in ICS
            all(isdefined(E, i) for E in (E3, E4, Ep)) || continue

            ic4 = getfield(E4, i)()
            ic3 = getfield(E3, i)()
            icp = getfield(Ep, i)()

            x, u, μ = ic4.q[1:3], ic4.q[4], ic4.params.μ

            # position
            @test ic3.q ≈ x
            @test icp.q ≈ x
            # parallel velocity, recovered from each family's own state
            @test E3.u(0.0, ic3.q, ic3.p) ≈ u  atol = 1E-12
            @test icp.v' * Ep.b(0.0, icp.q) ≈ u  atol = 1E-12
            # magnetic moment
            @test ic3.params.μ == μ
            @test icp.params.μ == μ
        end
    end

    # And that there is something to compare: every shared equilibrium carries the same set of names.
    for name in sharedmodules()
        Es = (getfield(G3, Symbol(name)), getfield(G4, Symbol(name)), getfield(P3, Symbol(name)))
        for i in ICS
            @test length(unique(isdefined(E, i) for E in Es)) == 1
        end
    end

end


@safetestset "Model families: GuidingCenter3d and GuidingCenter4d are the same model                 " begin

    using ChargedParticleDynamics
    using ..ModelAgreementUtils

    # 500 steps at a step both resolve. The thresholds are two to three orders above the measured
    # values, which are 1E-13 and below in the curvilinear charts; the cartesian chart of the small
    # tokamak is looser because the toroidal transit appears there as a rotation of (x, y) that has
    # to be resolved rather than a coordinate that winds.
    #
    # The set is every equilibrium `scripts/study_model_agreement.jl` measures, at the same steps, so
    # the assertions and the measurements they are set from cannot drift apart.
    CASES = (("TokamakSmallCartesian",   10.0, 1E-6),
             ("TokamakSmallCylindrical", 10.0, 1E-9),
             ("TokamakSmallToroidal",    10.0, 1E-9),
             ("TokamakIterCylindrical",  0.01, 1E-9),
             ("SolovevIter",             0.01, 1E-9),
             ("SolovevIterXpoint",       0.01, 1E-9))

    for (name, step, tol) in CASES
        M3, M4 = getfield(G3, Symbol(name)), getfield(G4, Symbol(name))
        span = (0.0, 500 * step)

        for i in (:initial_conditions_barely_passing, :initial_conditions_deeply_trapped)
            isdefined(M4, i) || continue
            ic = getfield(M4, i)()
            x, u, μ = ic.q[1:3], ic.q[4], ic.params.μ

            r3 = final_g3(M3, x, u, μ, span, step)
            r4 = final_g4(M4, x, u, μ, span, step)

            @test reldiff(r3[1], r4[1]) < tol
            @test abs(r3[2] - r4[2]) / max(abs(r4[2]), 1E-30) < tol
        end
    end

end


@safetestset "Model families: the Pauli particle tracks the guiding centre on its slow manifold      " begin

    using ChargedParticleDynamics
    using ..ModelAgreementUtils

    # Two statements, and the second is the substantive one.
    #
    # `v = u b⃗` is the slow manifold to lowest order and tracks the guiding centre to 2E-2 at worst,
    # the worst case being the ITER Solov'ev where ρ/L ≈ 0.24 and the expansion is barely valid.
    #
    # Supplying the guiding centre velocity itself removes the O(v_drift) mismatch in the initial
    # velocity and brings *every* equilibrium below 1E-3, independently of ρ/L. That is what says the
    # two models share a slow manifold rather than merely resembling each other: the four decades of
    # spread in the first column are the initial condition, not the dynamics.
    #
    # The set is every equilibrium `scripts/study_model_agreement.jl` measures, at the same steps.
    # `TokamakSmallCartesian` is the row where the drift correction is *worse* than the lowest-order
    # condition — 2.3E-5 against 1.9E-5 on `barely_passing`, 1.0E-5 against 6.9E-6 on
    # `deeply_trapped` — because `v_gc` fixes the `O(v_drift)` mismatch in the velocity while leaving
    # the `O(ρ)` mismatch in the position, and this is where the residual is already down at the
    # position term. It is kept for exactly that reason: the assertion below is that the corrected
    # condition lands under 1E-3 everywhere, not that it improves everywhere.
    CASES = (("TokamakSmallCartesian",   10.0),
             ("TokamakSmallCylindrical", 10.0),
             ("TokamakSmallToroidal",    10.0),
             ("TokamakIterCylindrical",  0.01),
             ("SolovevIter",             0.01),
             ("SolovevIterXpoint",       0.01))

    for (name, step) in CASES
        M4, Mp = getfield(G4, Symbol(name)), getfield(P3, Symbol(name))
        span = (0.0, 500 * step)

        for i in (:initial_conditions_barely_passing, :initial_conditions_deeply_trapped)
            isdefined(M4, i) || continue
            ic = getfield(M4, i)()
            x, u, μ = ic.q[1:3], ic.q[4], ic.params.μ

            ref = final_g4(M4, x, u, μ, span, step)[1]

            parallel = final_pauli(Mp, x, u * Mp.b⃗(0.0, x), μ, span, step)[1]
            @test reldiff(parallel, ref) < 3E-2

            vgc = zeros(4)
            M4.guiding_center_4d_v(vgc, 0.0, [x..., u], (μ = μ,))
            drifting = final_pauli(Mp, x, vgc[1:3], μ, span, step)[1]
            @test reldiff(drifting, ref) < 1E-3
        end
    end

end
