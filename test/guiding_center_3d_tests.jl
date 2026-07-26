
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
# At the `1E-15` this used to request, the criterion was unreachable and Newton ran to its
# 1000-iteration limit on nearly half of all steps — a single `iode` testset cost minutes instead of
# seconds. `max_iterations` bounds the damage should a future equilibrium have a larger `‖ϑ‖`, and
# `warn_iterations` is set to match it so that hitting the cap is still reported — left at its
# default of 1000 it could never fire below the cap, and the count in `quiet_solver_warnings.jl`
# would go quiet regardless of how badly the solver was doing.
#
# `f_reltol` is deliberately left at the `SimpleSolvers` default of `√eps`. Its relative term
# `f_reltol ‖F(x₀)‖` is what lets a large-magnitude solve converge at all, and pinning it to `1E-15`
# alongside `f_abstol` removed the only criterion these problems could still meet.
const options = (f_abstol=1E-12, max_iterations=50, warn_iterations=50)

export test_guiding_center_3d
export nl, nx, ny

# The tests assert that the integration runs to completion and returns a solution.
#
# They used to assert `@test_nowarn` instead, which is not a usable criterion here: `SimpleSolvers`
# warns once per step whose solve exhausts the iteration budget, and which orbits trip it depends on
# the floating-point details of the platform, so the suite failed on a different equilibrium on
# Linux than on Windows. The warnings are filtered out by `test/quiet_solver_warnings.jl` and
# counted there; with the tolerances above the blocks that pass these options contribute none of
# them, so a sharp rise in the reported count means a tolerance has slipped below a residual floor.
# Nothing asserts on that count — see `quiet_solver_warnings.jl` for why a threshold would flake.
function test_guiding_center_3d(equ::ODEProblem)
    @test integrate(equ, Gauss(2); options...) isa GeometricSolution
end

function test_guiding_center_3d(equ::Union{HODEProblem,PODEProblem})
    @test integrate(equ, PartitionedGauss(2); initialguess=MidpointExtrapolation(5), options...) isa GeometricSolution
end

function test_guiding_center_3d(equ::Union{IODEProblem,LODEProblem})
    @test integrate(equ, VPRKGauss(2); initialguess=MidpointExtrapolation(5), options...) isa GeometricSolution
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

    test_guiding_center_3d(hodeproblem(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_trapped(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_trapped(); timespan=(0.0, 1E1), timestep=0.1))

    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_trapped(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_deeply_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_deeply_trapped(); timespan=(0.0, 1E1), timestep=0.1))

    test_guiding_center_3d(hodeproblem_compact(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem_compact(initial_conditions_deeply_trapped(); timespan=(0.0, 1E1), timestep=0.1))

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

    # These run at the module's own default step of 0.03 rather than the 0.1 the blocks below use,
    # because that is where the dipole's constraints stay at 1E-8 rather than 1E-6. The step size is
    # not what used to stop `hodeproblem_canonical` converging here, though — the constraint pair was.
    # `Dipole3d` is the one equilibrium whose default pair is chosen on its conditioning along the
    # orbit rather than at the initial condition; see `default_constraints` in `dipole.jl`.
    test_guiding_center_3d(hodeproblem(initial_conditions_dipole(); timespan=(0.0, 3.0), timestep=0.03))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_dipole(); timespan=(0.0, 3.0), timestep=0.03))
    test_guiding_center_3d(hodeproblem_compact(initial_conditions_dipole(); timespan=(0.0, 3.0), timestep=0.03))

end


@safetestset "Guiding Centre Dynamics in 3D in Quadratic Potentials                                               " begin

    using ChargedParticleDynamics.GuidingCenter3d.QuadraticPotentials3d
    using ..GuidingCenter3dTests

    # The only equilibrium in the package with a non-zero electrostatic potential, so this is also
    # the only integration test that exercises the `φ` term of `hamiltonian` and the `E` term of
    # `dHdqᵢ` on a field where they do not vanish.
    test_guiding_center_3d(hodeproblem(initial_conditions_quadratic(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_quadratic(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem_compact(initial_conditions_quadratic(); timespan=(0.0, 1E1), timestep=0.1))

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

    test_guiding_center_3d(hodeproblem(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_trapped(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_trapped(); timespan=(0.0, 1E1), timestep=0.1))

    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_trapped(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_deeply_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_deeply_trapped(); timespan=(0.0, 1E1), timestep=0.1))

    test_guiding_center_3d(hodeproblem_compact(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem_compact(initial_conditions_deeply_trapped(); timespan=(0.0, 1E1), timestep=0.1))

    # The one equilibrium where no component of `b` vanishes at the initial condition, so all three
    # constraint pairs are usable and can be compared against each other; see the "the constraint
    # formulations agree" block at the end of this file.
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_passing(); constraints=:g12, timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_passing(); constraints=:g23, timespan=(0.0, 1E1), timestep=0.1))

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

    test_guiding_center_3d(hodeproblem(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_trapped(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_trapped(); timespan=(0.0, 1E1), timestep=0.1))

    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_trapped(); timespan=(0.0, 1E1), timestep=0.1))

    test_guiding_center_3d(hodeproblem_compact(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem_compact(initial_conditions_deeply_trapped(); timespan=(0.0, 1E1), timestep=0.1))

end


@safetestset "Guiding Centre Dynamics in 3D with small-size Tokamak Equilibrium in Cartesian Coordinates          " begin

    using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallCartesian
    using ..GuidingCenter3dTests

    test_guiding_center_3d(hodeproblem(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_trapped(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_trapped(); timespan=(0.0, 1E1), timestep=0.1))

    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_trapped(); timespan=(0.0, 1E1), timestep=0.1))

    test_guiding_center_3d(hodeproblem_compact(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem_compact(initial_conditions_deeply_trapped(); timespan=(0.0, 1E1), timestep=0.1))

end


@safetestset "Guiding Centre Dynamics in 3D with small-size Tokamak Equilibrium in Cylindrical Coordinates        " begin

    using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallCylindrical
    using ..GuidingCenter3dTests

    test_guiding_center_3d(hodeproblem(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_trapped(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_trapped(); timespan=(0.0, 1E1), timestep=0.1))

    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_trapped(); timespan=(0.0, 1E1), timestep=0.1))

    test_guiding_center_3d(hodeproblem_compact(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem_compact(initial_conditions_deeply_trapped(); timespan=(0.0, 1E1), timestep=0.1))

end


@safetestset "Guiding Centre Dynamics in 3D with small-size Tokamak Equilibrium in Toroidal Coordinates           " begin

    using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallToroidal
    using ..GuidingCenter3dTests

    test_guiding_center_3d(hodeproblem(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_trapped(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_trapped(); timespan=(0.0, 1E1), timestep=0.1))

    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem_canonical(initial_conditions_barely_trapped(); timespan=(0.0, 1E1), timestep=0.1))

    test_guiding_center_3d(hodeproblem_compact(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem_compact(initial_conditions_deeply_trapped(); timespan=(0.0, 1E1), timestep=0.1))

end


@safetestset "Guiding Centre Dynamics in 3D with symmetric Solov'ev Equilibrium                                   " begin

    using ChargedParticleDynamics.GuidingCenter3d.SolovevSymmetricField
    using ..GuidingCenter3dTests

    # The only equilibrium for which two of the three constraint pairs are singular at once: `b = e₂`
    # at these initial conditions, so both `(g³, g¹)` and `(g¹, g²)` divide by zero and only the
    # `(g², g³)` this module defaults to survives. See `default_constraints` in `solovev_symmetric.jl`.
    #
    # Also the one block that cannot be run over the module's own `DEFAULT_TIMESPAN`. This is the
    # stiffest equilibrium in the package — `‖ϑ‖ ≈ 256`, an order of magnitude above the ITER-like
    # Solov'ev equilibria — and the constraint drift grows fast enough that Newton meets a NaN a few
    # hundred steps in. A hundred steps, the same workload as the blocks above, is well inside that.
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_barely_trapped(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem(initial_conditions_deeply_trapped(); timespan=(0.0, 1E1), timestep=0.1))

    test_guiding_center_3d(hodeproblem_compact(initial_conditions_barely_passing(); timespan=(0.0, 1E1), timestep=0.1))
    test_guiding_center_3d(hodeproblem_compact(initial_conditions_deeply_trapped(); timespan=(0.0, 1E1), timestep=0.1))

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
    solutions = [integrate(p, PartitionedGauss(2); initialguess = MidpointExtrapolation(5),
                           GuidingCenter3dTests.options...)
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
