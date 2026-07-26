
using SafeTestsets

module GuidingCenter4dTests

using GeometricIntegrators
using Test

const nl = 100
const nx = 10
const ny = 10

# See `guiding_center_3d_tests.jl` for why `f_abstol` has to stay above the residual's round-off
# floor and why `f_reltol` is left at its default.
const options = (f_abstol=1E-12, max_iterations=50, warn_iterations=50)

export test_guiding_center_4d
export nl, nx, ny

# Asserts that the integration returns a solution rather than `@test_nowarn`; see the comment on
# the corresponding functions in `guiding_center_3d_tests.jl` for why.
function test_guiding_center_4d(equ::ODEProblem)
    @test integrate(equ, Gauss(2); options...) isa GeometricSolution
end

function test_guiding_center_4d(equ::Union{IODEProblem,LODEProblem})
    @test integrate(equ, VPRKGauss(2); options...) isa GeometricSolution
end

end



@safetestset "Guiding Centre Dynamics in 4D with ITER-like Solov'ev Equilibrium with X-Point                      " begin

    using ChargedParticleDynamics.GuidingCenter4d.SolovevIterXpoint
    using ..GuidingCenter4dTests

    test_guiding_center_4d(odeproblem(initial_conditions_trapped(); timestep=1E4, timespan=(0, 1E6)))
    test_guiding_center_4d(odeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_trapped()))
    # test_guiding_center_4d(loop_odeproblem(nl), Δt=1.)
    # test_guiding_center_4d(surface_odeproblem(nx, ny), Δt=1.)

    test_guiding_center_4d(iodeproblem(initial_conditions_trapped(); timestep=1E4, timespan=(0, 1E6)))
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_passing(); timestep=1E-1, timespan=(0, 1E1)))
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_trapped(); timestep=1E-1, timespan=(0, 1E1)))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_passing(); timestep=1E-1, timespan=(0, 1E1)))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_trapped(); timestep=1E-1, timespan=(0, 1E1)))
    # test_guiding_center_4d(loop_iodeproblem(nl), Δt=1.)
    # test_guiding_center_4d(surface_iodeproblem(nx, ny), Δt=1.)

end


@safetestset "Guiding Centre Dynamics in 4D with medium-size Tokamak Equilibrium in Cartesian Coordinates         " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCartesian
    using ..GuidingCenter4dTests

    test_guiding_center_4d(odeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_trapped()))
    # test_guiding_center_4d(loop_odeproblem(nl))
    # test_guiding_center_4d(surface_odeproblem(nx, ny))

    test_guiding_center_4d(iodeproblem(initial_conditions_barely_passing(); timestep=1E-1, timespan=(0, 1E1)))
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_trapped(); timestep=1E-1, timespan=(0, 1E1)))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_passing(); timestep=1E-1, timespan=(0, 1E1)))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_trapped(); timestep=1E-1, timespan=(0, 1E1)))
    # test_guiding_center_4d(loop_iodeproblem(nl))
    # test_guiding_center_4d(surface_iodeproblem(nx, ny))

end


@safetestset "Guiding Centre Dynamics in 4D with medium-size Tokamak Equilibrium in Cylindrical Coordinates       " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCylindrical
    using ..GuidingCenter4dTests

    test_guiding_center_4d(odeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_trapped()))
    # test_guiding_center_4d(loop_odeproblem(nl))
    # test_guiding_center_4d(surface_odeproblem(nx, ny))

    test_guiding_center_4d(iodeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_trapped()))
    # test_guiding_center_4d(loop_iodeproblem(nl))
    # test_guiding_center_4d(surface_iodeproblem(nx, ny))

end


@safetestset "Guiding Centre Dynamics in 4D with small-size Tokamak Equilibrium in Cartesian Coordinates          " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCartesian
    using ..GuidingCenter4dTests

    test_guiding_center_4d(odeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_trapped()))
    # test_guiding_center_4d(loop_odeproblem(nl))
    # test_guiding_center_4d(surface_odeproblem(nx, ny))

    # the variational integrators do not converge at the default time step of Δt = 500
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_passing(); timestep=10.0, timespan=(0, 1E3)))
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_trapped(); timestep=10.0, timespan=(0, 1E3)))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_passing(); timestep=10.0, timespan=(0, 1E3)))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_trapped(); timestep=10.0, timespan=(0, 1E3)))
    # test_guiding_center_4d(loop_iodeproblem(nl; timestep = 10., timespan = [0, 1E3]))
    # test_guiding_center_4d(surface_iodeproblem(nx, ny; timestep = 10., timespan = [0, 1E3]))

end


@safetestset "Guiding Centre Dynamics in 4D with small-size Tokamak Equilibrium in Cylindrical Coordinates        " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCylindrical
    using ..GuidingCenter4dTests

    test_guiding_center_4d(odeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_trapped()))
    # test_guiding_center_4d(loop_odeproblem(nl))
    # test_guiding_center_4d(surface_odeproblem(nx, ny))

    test_guiding_center_4d(iodeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_trapped()))
    # test_guiding_center_4d(loop_iodeproblem(nl; timestep = 10., timespan = [0, 1E3]))
    # test_guiding_center_4d(surface_iodeproblem(nx, ny; timestep = 10., timespan = [0, 1E3]))

    # The κ-dependent "dg" formulation. Its `ḡ` is checked against finite differences in
    # `structure_tests.jl`, but until now it had never been integrated — and `κ = 0`, which is the
    # default, reduces it to the plain `iodeproblem` and so would not exercise the κ terms at all.
    test_guiding_center_4d(iodeproblem_dg(initial_conditions_barely_passing(); κ=0.0, timespan=(0.0, 1E3), timestep=10.0))
    test_guiding_center_4d(iodeproblem_dg(initial_conditions_barely_passing(); κ=0.5, timespan=(0.0, 1E3), timestep=10.0))
    test_guiding_center_4d(iodeproblem_dg(initial_conditions_barely_passing(); κ=1.0, timespan=(0.0, 1E3), timestep=10.0))

end


@safetestset "Guiding Centre Dynamics in 4D with small-size Tokamak Equilibrium in Toroidal Coordinates           " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakSmallToroidal
    using ..GuidingCenter4dTests

    test_guiding_center_4d(odeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_trapped()))
    # test_guiding_center_4d(loop_odeproblem(nl))
    # test_guiding_center_4d(surface_odeproblem(nx, ny))

    test_guiding_center_4d(iodeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_trapped()))
    # test_guiding_center_4d(loop_iodeproblem(nl))
    # test_guiding_center_4d(surface_iodeproblem(nx, ny))

end


@safetestset "Guiding Centre Dynamics in 4D with symmetric Solov'ev Equilibrium                                   " begin

    using ChargedParticleDynamics.GuidingCenter4d.SolovevSymmetricField
    using ..GuidingCenter4dTests

    test_guiding_center_4d(odeproblem(initial_conditions_barely_passing(); timestep=1E2, timespan=(0, 1E4)))
    # the barely trapped orbit needs a smaller time step for the solver to converge
    test_guiding_center_4d(odeproblem(initial_conditions_barely_trapped(); timestep=1E1, timespan=(0, 1E4)))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_trapped(); timestep=1E1, timespan=(0, 1E3)))

    test_guiding_center_4d(iodeproblem(initial_conditions_barely_passing(); timestep=1E0, timespan=(0, 1E2)))
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_trapped(); timestep=1E0, timespan=(0, 1E2)))

end


# These two equilibria carry only a Poincaré loop and surface parameterisation, no point initial
# condition, so what there is to test is the loop and surface problems. Both blocks used to be
# entirely commented out and ran zero tests: their calls passed the sampling counts positionally —
# `loop_odeproblem(nl)` — which the `PoincareInvariants` 0.5 rewrite replaced with keyword-only
# constructors that take the sample count at the invariant rather than at the problem.
@safetestset "Guiding Centre Dynamics in 4D with symmetric Equilibrium                                            " begin

    using ChargedParticleDynamics.GuidingCenter4d.SymmetricField
    using ..GuidingCenter4dTests

    test_guiding_center_4d(loop_odeproblem(; timespan=(0.0, 1E2), timestep=1.0))
    test_guiding_center_4d(surface_odeproblem(; timespan=(0.0, 1E2), timestep=1.0))

    test_guiding_center_4d(loop_iodeproblem(; timespan=(0.0, 1E2), timestep=1.0))
    test_guiding_center_4d(surface_iodeproblem(; timespan=(0.0, 1E2), timestep=1.0))

end


@safetestset "Guiding Centre Dynamics in 4D with Theta Pinch Equilibrium                                          " begin

    using ChargedParticleDynamics.GuidingCenter4d.ThetaPinchField
    using ..GuidingCenter4dTests

    # The theta pinch defines a loop but no surface.
    test_guiding_center_4d(loop_odeproblem(; timespan=(0.0, 1E2), timestep=1.0))
    test_guiding_center_4d(loop_iodeproblem(; timespan=(0.0, 1E2), timestep=1.0))

end
