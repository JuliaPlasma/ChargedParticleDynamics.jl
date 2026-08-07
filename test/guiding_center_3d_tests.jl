
using SafeTestsets

module GuidingCenter3dTests

using GeometricIntegrators
using Test

const nl = 100
const nx = 10
const ny = 10

# `f_abstol` is an *absolute* bound on the residual, so it has to sit above the round-off floor of
# the problem. The variational and Pauli residuals contain the momentum equation `p = ϑ(q)`, whose
# scale is set by the one-form `ϑ = A + u b`, so that floor is `‖ϑ‖ eps`: the ITER-like Solov'ev
# equilibria have `‖ϑ‖ ≈ 16` (`≈ 54` for the symmetric one), i.e. a floor of `3.5E-15` (`1.2E-14`).
# An `f_abstol` below that floor is unreachable, and the solver then spends its whole budget failing
# to meet it. The library default does not clear it either: `GeometricIntegratorsBase` scales its own
# with the stage system, `f_abstol = max(8, solversize(method, problem)) * eps(datatype(problem))`,
# which is `1.8E-15`–`2.7E-15` for the methods used here.
#
# `max_iterations` bounds the damage should a future equilibrium have a larger `‖ϑ‖`, and
# `warn_iterations` matches it so that reaching the cap is reported — at its default of 1000 it could
# never fire below the cap. Stagnation is not what the cap is for: `SimpleSolvers` stops a solve that
# cannot move its iterate after `max_stalls = 2` steps. It is for solves that are genuinely still
# iterating, which on `GuidingCenter4d.SolovevSymmetricField` neither converge nor stall.
#
# `f_reltol` is deliberately left at the `SimpleSolvers` default of `√eps`. Its relative term
# `f_reltol ‖F(x₀)‖` is what lets a large-magnitude solve converge at all, and pinning it to `1E-15`
# alongside `f_abstol` removed the only criterion these problems could still meet.
const options = (f_abstol=1E-12, max_iterations=50, warn_iterations=50)

export test_guiding_center_3d
export nl, nx, ny

# ---------------------------------------------------------------------------------------------
# Time step and time span limits of the 3D guiding centre equilibria
# ---------------------------------------------------------------------------------------------
#
# Every call in this file runs at its module's own `DEFAULT_TIMESPAN`/`DEFAULT_TIMESTEP` — a
# thousand steps in every model — except the two noted below. Each module's constants carry the
# measurement behind its own step; what follows is what belongs to the file rather than to one
# module, and none of it is recoverable from the code.
#
# **A formulation, not a module, carries an exception.** Two formulations cannot hold their module's
# declared span, and only those two take a shorter one:
#
#   TokamakMediumCartesian.hodeproblem_canonical loses `deeply_passing` between t = 20 and t = 50
#     (max|g¹| 6.8E-8 at t = 10, 5.3E-5 at 20, 4.3E+20 at 50, NaN at the declared 100), while
#     `hodeproblem` and `hodeproblem_compact` hold the full span at 1E-8.
#   SolovevSymmetricField.hodeproblem loses the orbit past t ≈ 10 at *any* step size tried — a finer
#     step does not rescue it — as does `hodeproblem_canonical`. `hodeproblem_compact` is the only
#     formulation that survives the declared example here, which is why the two are not levelled to
#     a common span: the difference is the thing worth testing.
#
# **`Dipole3d` is deliberately not at its most accurate step.** `Δt = 0.1` exhibits the constraint
# drift rather than minimising it: the three formulations separate legibly there (4.8E-6, 7.1E-4,
# 8.1E-6 in max|g¹|) where at Δt = 0.03 they all fall into the 1E-8 range and the difference stops
# being visible.
#
# **The `TokamakSmall*` charts do not all take the Δt = 500 their 4D and Pauli counterparts use.**
# `Cylindrical` does, in all three formulations (4.3E-3 relative energy; constraints 6.0E-9, 1.8E-6,
# 3.8E-9), and declares it. `Toroidal` does in `hodeproblem` and `hodeproblem_compact` (4.8E-3) but
# not in `hodeproblem_canonical`, which throws a `SingularException` there and takes 250, the
# largest step it tolerates.
#
# `Cartesian` cannot, and **the constraint pair is not what stops it**. Its `:g31` pair is singular
# at the initial condition (λₒ = 0), leaving `:g12` and `:g23`, and at Δt = 500 both fail in every
# formulation: `:g12` throws a NaN direction in all three, and `:g23` — the pair it declares — throws
# a `SingularException` in `hodeproblem` and a NaN in the other two. Giving it the `:g12` its
# siblings use would change nothing.
#
# One subtlety when reproducing that: `hodeproblem_compact(ic)` and
# `hodeproblem_compact(ic; constraints = s)` are different problems. Without the keyword the compact
# form is the regularised, pair-independent variant and the more robust of the two — at Δt = 500 on
# `Cartesian` it is the only thing that runs, and then only to 4.2E-1 in relative energy. With an
# explicit pair it throws like the rest. The tests call it without.
#
# Every step and span above was chosen against energy, and against the constraints where the 3D
# model makes them available — not against solver silence, which proves less than it looks; see the
# header of `quiet_solver_warnings.jl`.

# The tests assert that the integration runs to completion and returns a solution. Solver silence is
# asserted separately and for the file as a whole, by `runtests.jl`, which is the stronger statement
# and does not depend on which orbit happens to trip a warning on which platform. See
# `test/quiet_solver_warnings.jl`.
function test_guiding_center_3d(equ::ODEProblem)
    @test integrate(equ, Gauss(2); options...) isa GeometricSolution
end

# These take the default initial guess. `PartitionedGauss` and `VPRKGauss` both declare
# `default_iguess() = HermiteExtrapolation()`, and requesting `MidpointExtrapolation(5)` instead
# costs 70% of the runtime of a 3D guiding centre integration for nothing: Newton converges in one
# iteration under either, and the trajectories agree to round-off (bit-identical for most equilibria,
# 1.5E-14 relative at worst). The Pauli theta pinch is the one place that needs it, for a reason
# recorded in `pauli_particle_3d_tests.jl`.
function test_guiding_center_3d(equ::Union{HODEProblem,PODEProblem})
    @test integrate(equ, PartitionedGauss(2); options...) isa GeometricSolution
end

function test_guiding_center_3d(equ::Union{IODEProblem,LODEProblem})
    @test integrate(equ, VPRKGauss(2); options...) isa GeometricSolution
end

end



@safetestset "Guiding Centre Dynamics in 3D with ITER-like Solov'ev Equilibrium with X-Point                      " begin

    using ChargedParticleDynamics.GuidingCenter3d.SolovevIterXpoint
    using ..GuidingCenter3dTests


    # test_guiding_center_3d(ode(initial_conditions_trapped(); timestep = 1E4, timespan = (0, 1E6)))
    # test_guiding_center_3d(ode(initial_conditions_barely_passing()))
    # test_guiding_center_3d(ode(initial_conditions_barely_trapped()))
    # test_guiding_center_3d(ode(initial_conditions_deeply_passing()))
    # test_guiding_center_3d(ode(initial_conditions_deeply_trapped()))
    # test_guiding_center_3d(guiding_center_3d_loop_ode(nl), Δt=1.)
    # test_guiding_center_3d(guiding_center_3d_surface_ode(nx, ny), Δt=1.)

    test_guiding_center_3d(hodeproblem(initial_conditions_barely_passing()))
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_trapped()))

    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_passing()))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_trapped()))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_deeply_passing()))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_deeply_trapped()))

    test_guiding_center_3d(hodeproblem_compact(initial_conditions_barely_passing()))
    test_guiding_center_3d(hodeproblem_compact(initial_conditions_deeply_trapped()))

    # test_guiding_center_3d(iode(initial_conditions_trapped(); timestep = 1E4, timespan = (0, 1E6)))
    # test_guiding_center_3d(iode(initial_conditions_barely_passing()))
    # test_guiding_center_3d(iode(initial_conditions_barely_trapped()))
    # test_guiding_center_3d(iode(initial_conditions_deeply_passing()))#; timestep = 1E-2, timespan = (0, 1E1)
    # test_guiding_center_3d(iode(initial_conditions_deeply_trapped()))#; timestep = 1E-2, timespan = (0, 1E1)
    # test_guiding_center_3d(guiding_center_3d_loop_iode(nl), Δt=1.)
    # test_guiding_center_3d(guiding_center_3d_surface_iode(nx, ny), Δt=1.)

end


# `Dipole3d` and `QuadraticPotentials3d` are the two equilibria that are *not* blocked by the
# `b₁ = 0` singularity of the `(g³, g¹)` constraint pair — `b₁` is -0.41 and 2E-3 at their
# respective initial conditions — and until now neither had a test block at all. With these two,
# the 3D model's dynamics is exercised on four equilibria of eleven rather than two.
@safetestset "Guiding Centre Dynamics in 3D in a Dipole Field                                                     " begin

    using ChargedParticleDynamics.GuidingCenter3d.Dipole3d
    using ..GuidingCenter3dTests

    # This module's step is deliberately coarse enough to *show* the constraint drift rather than to
    # minimise it: at `Δt = 0.1` the three formulations separate legibly (4.8E-6, 7.1E-4 and 8.1E-6
    # in `max|g¹|`), where at `Δt = 0.03` they all fall to the 1E-8 range and the difference between
    # them stops being visible. What governs whether `hodeproblem_canonical` converges here is the
    # constraint pair rather than the step: `Dipole3d` is the one equilibrium whose pair is chosen on
    # its conditioning along the orbit rather than at the initial condition. See `dipole.jl`.
    test_guiding_center_3d(hodeproblem(initial_conditions_dipole()))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_dipole()))
    test_guiding_center_3d(hodeproblem_compact(initial_conditions_dipole()))

end


@safetestset "Guiding Centre Dynamics in 3D in Quadratic Potentials                                               " begin

    using ChargedParticleDynamics.GuidingCenter3d.QuadraticPotentials3d
    using ..GuidingCenter3dTests

    # The only equilibrium in the package with a non-zero electrostatic potential, so this is also
    # the only integration test that exercises the `φ` term of `hamiltonian` and the `E` term of
    # `dHdqᵢ` on a field where they do not vanish.
    #
    # The declared span was cut from 2.5E4 to 5E2 — fifty thousand steps down to a thousand — because
    # `hodeproblem_canonical` diverges to `Inf` with 26 capped steps over the longer run. All three
    # formulations sit between 3.1E-9 and 7.4E-9 at the corrected default.
    test_guiding_center_3d(hodeproblem(initial_conditions_quadratic()))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_quadratic()))
    test_guiding_center_3d(hodeproblem_compact(initial_conditions_quadratic()))

end


@safetestset "Guiding Centre Dynamics in 3D with medium-size Tokamak Equilibrium in Cartesian Coordinates         " begin

    using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCartesian
    using ..GuidingCenter3dTests

    # test_guiding_center_3d(ode(initial_conditions_barely_passing()))
    # test_guiding_center_3d(ode(initial_conditions_barely_trapped()))
    # test_guiding_center_3d(ode(initial_conditions_deeply_passing()))
    # test_guiding_center_3d(ode(initial_conditions_deeply_trapped()))
    # test_guiding_center_3d(guiding_center_3d_loop_ode(nl))
    # test_guiding_center_3d(guiding_center_3d_surface_ode(nx, ny))

    test_guiding_center_3d(hodeproblem(initial_conditions_barely_passing()))
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_trapped()))

    # Only the canonicalised form needs a shorter span here, so only it gets one: it loses
    # `deeply_passing` between `t = 20` and `t = 50` (5.3E-5 in `max|g¹|` at 20, 4.3E+20 at 50,
    # `NaN` at the default 100), where the other two formulations hold the full span at 1E-8. The
    # exception belongs to the formulation that blocks, not to the module's declared example.
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_passing(); timespan=(0.0, 1E1)))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_trapped(); timespan=(0.0, 1E1)))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_deeply_passing(); timespan=(0.0, 1E1)))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_deeply_trapped(); timespan=(0.0, 1E1)))

    test_guiding_center_3d(hodeproblem_compact(initial_conditions_barely_passing()))
    test_guiding_center_3d(hodeproblem_compact(initial_conditions_deeply_trapped()))

    # The one equilibrium where no component of `b` vanishes at the initial condition, so all three
    # constraint pairs are usable and can be compared against each other; see the "the constraint
    # formulations agree" block at the end of this file.
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_passing(); constraints=:g12))
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_passing(); constraints=:g23))

    # test_guiding_center_3d(iode(initial_conditions_barely_passing()))
    # test_guiding_center_3d(iode(initial_conditions_barely_trapped()))
    # test_guiding_center_3d(iode(initial_conditions_deeply_passing()))
    # test_guiding_center_3d(iode(initial_conditions_deeply_trapped()))
    # test_guiding_center_3d(guiding_center_3d_loop_iode(nl))
    # test_guiding_center_3d(guiding_center_3d_surface_iode(nx, ny))

end


@safetestset "Guiding Centre Dynamics in 3D with medium-size Tokamak Equilibrium in Cylindrical Coordinates       " begin

    using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCylindrical
    using ..GuidingCenter3dTests

    test_guiding_center_3d(hodeproblem(initial_conditions_barely_passing()))
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_trapped()))

    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_passing()))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_trapped()))

    test_guiding_center_3d(hodeproblem_compact(initial_conditions_barely_passing()))
    test_guiding_center_3d(hodeproblem_compact(initial_conditions_deeply_trapped()))

end


@safetestset "Guiding Centre Dynamics in 3D with small-size Tokamak Equilibrium in Cartesian Coordinates          " begin

    using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallCartesian
    using ..GuidingCenter3dTests

    test_guiding_center_3d(hodeproblem(initial_conditions_barely_passing()))
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_trapped()))

    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_passing()))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_trapped()))

    test_guiding_center_3d(hodeproblem_compact(initial_conditions_barely_passing()))
    test_guiding_center_3d(hodeproblem_compact(initial_conditions_deeply_trapped()))

end


@safetestset "Guiding Centre Dynamics in 3D with small-size Tokamak Equilibrium in Cylindrical Coordinates        " begin

    using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallCylindrical
    using ..GuidingCenter3dTests

    test_guiding_center_3d(hodeproblem(initial_conditions_barely_passing()))
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_trapped()))

    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_passing()))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_trapped()))

    test_guiding_center_3d(hodeproblem_compact(initial_conditions_barely_passing()))
    test_guiding_center_3d(hodeproblem_compact(initial_conditions_deeply_trapped()))

end


@safetestset "Guiding Centre Dynamics in 3D with small-size Tokamak Equilibrium in Toroidal Coordinates           " begin

    using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallToroidal
    using ..GuidingCenter3dTests

    # This module runs at `Δt = 500`, matching its 4D and Pauli counterparts. Only
    # `hodeproblem_canonical` cannot: it throws a `SingularException` there on either regular pair.
    # 250 is the largest step it tolerates, so it takes that and the others run at the default.

    test_guiding_center_3d(hodeproblem(initial_conditions_barely_passing()))
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_trapped()))

    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_passing(); timestep=250.0, timespan=(0.0, 2.5E5)))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_trapped(); timestep=250.0, timespan=(0.0, 2.5E5)))

    test_guiding_center_3d(hodeproblem_compact(initial_conditions_barely_passing()))
    test_guiding_center_3d(hodeproblem_compact(initial_conditions_deeply_trapped()))

end


@safetestset "Guiding Centre Dynamics in 3D with symmetric Solov'ev Equilibrium                                   " begin

    using ChargedParticleDynamics.GuidingCenter3d.SolovevSymmetricField
    using ..GuidingCenter3dTests

    # The only equilibrium for which two of the three constraint pairs are singular at once: `b = e₂`
    # at these initial conditions, so both `(g³, g¹)` and `(g¹, g²)` divide by zero and only the
    # `(g², g³)` this module defaults to survives. See `default_constraints` in `solovev_symmetric.jl`.
    #
    # `hodeproblem` carries a shorter span here and `hodeproblem_compact` does not, because it is the
    # former that blocks: `hodeproblem` and `hodeproblem_canonical` lose this orbit somewhere between
    # `t = 10` and `t = 50` at *any* step size tried — 2.4E+6 and 1.5E+228 over the module's declared
    # `t ∈ [0, 1E3]` — and a finer step does not rescue them. `hodeproblem_compact` holds that span at
    # 1.1E-7 and so runs at the default. That only the compact form survives a long run here is the
    # thing worth testing, and it is visible precisely because the two are not levelled to the same
    # span. As in the 4D model, this is the worst-conditioned equilibrium of the set.
    #
    # Also the one block that cannot be run over the module's own `DEFAULT_TIMESPAN`. This is the
    # stiffest equilibrium in the package — `‖ϑ‖ ≈ 256`, an order of magnitude above the ITER-like
    # Solov'ev equilibria — and the constraint drift grows fast enough that Newton meets a NaN a few
    # hundred steps in. A hundred steps, the same workload as the blocks above, is well inside that.
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_passing(); timespan=(0.0, 1E1)))
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_trapped(); timespan=(0.0, 1E1)))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_passing(); timespan=(0.0, 1E1)))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_trapped(); timespan=(0.0, 1E1)))

    test_guiding_center_3d(hodeproblem_compact(initial_conditions_barely_passing()))
    test_guiding_center_3d(hodeproblem_compact(initial_conditions_deeply_trapped()))

    # `hodeproblem_canonical` is *not* exercised here. It starts fine — the pair is well conditioned,
    # `λₒ ≈ -228` and never falls below -173 along the orbit — but the ∂λ/∂q, ∂λ/∂p terms it adds
    # amplify the constraint drift enough that Newton hits a NaN after some thirty steps. That is a
    # property of the canonicalised formulation on this equilibrium, not of the constraint pair: the
    # same orbit runs to completion with `hodeproblem` and `hodeproblem_compact`. Section 2 of
    # `scripts/study_guiding_center_3d_conditioning.jl` shows the same thing quantitatively.

end


# `SymmetricField` and `ThetaPinchField` have no integration test block. Both are Poincaré-invariant
# fixtures: they ship `f_loop`/`f_surface` rather than `initial_conditions_*`, and the blocks that
# used to stand here called `guiding_center_3d_loop_ode` and `guiding_center_3d_surface_ode`, which
# this package no longer defines. Their constraints are exercised instead by the derivative tests in
# `test/structure_tests.jl`, which reach every module. Note that `b = e₃` in both, so `(g¹, g²)` is
# the only pair either of them can use.


@safetestset "Guiding Centre Dynamics in 3D: the constraint formulations agree                                    " begin

    using ChargedParticleDynamics.GuidingCenter3d
    using GeometricIntegrators
    using LinearAlgebra
    using Test
    using ..GuidingCenter3dTests

    # The substantive check on the three constraint pairs and the three formulations: where they are
    # all well conditioned they describe the same constrained system and must agree to the accuracy of
    # the nonlinear solve. `TokamakMediumCartesian` is the one equilibrium in the package for which no
    # component of `b` vanishes at the initial condition, so all three pairs are usable there at once.
    mx(ds) = maximum(abs(ds[i]) for i in eachindex(ds))

    M = GuidingCenter3d.TokamakMediumCartesian
    ic = M.initial_conditions_barely_passing()
    workload = (timespan = (0.0, 1E1), timestep = 0.1)

    problems = ((M.hodeproblem(ic; constraints = :g31, workload...)),
                (M.hodeproblem(ic; constraints = :g12, workload...)),
                (M.hodeproblem(ic; constraints = :g23, workload...)),
                (M.hodeproblem_canonical(ic; workload...)),
                (M.hodeproblem_compact(ic; workload...)),
                (M.hodeproblem_compact(ic; constraints = :g31, workload...)))

    # `options` is not exported by the helper module — the blocks above reach it through
    # `test_guiding_center_3d` rather than by name — so it is qualified here.
    solutions = [integrate(p, PartitionedGauss(2); GuidingCenter3dTests.options...)
                 for p in problems]

    reference = vcat(solutions[1].q[end], solutions[1].p[end])

    for (problem, solution) in zip(problems, solutions)
        y = vcat(solution.q[end], solution.p[end])
        @test all(isfinite, y)
        @test norm(y - reference) / norm(reference) < 1E-7

        # All three constraints vanish along the flow, whichever pair the problem retained.
        c = M.compute_constraints(solution)
        @test mx(c.g₁) < 1E-7
        @test mx(c.g₂) < 1E-7
        @test mx(c.g₃) < 1E-7

        _, e = M.compute_energy_error(solution.t, solution.q, solution.p, parameters(problem))
        @test mx(e) < 1E-7
    end

    # An unknown pair is a typo, not a silent fallback to the default.
    @test_throws MethodError M.hodeproblem(ic; constraints = :g13, workload...)

    # `:parallel` is the pair-independent regularisation of the compact form and means nothing for the
    # other two, which must say so rather than quietly picking something.
    @test_throws MethodError M.hodeproblem(ic; constraints = :parallel, workload...)
    @test_throws MethodError M.hodeproblem_canonical(ic; constraints = :parallel, workload...)
end
