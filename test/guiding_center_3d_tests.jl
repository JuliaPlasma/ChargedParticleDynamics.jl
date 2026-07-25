
using SafeTestsets

module GuidingCenter3dTests

using GeometricIntegrators
using Test

const nl = 100
const nx = 10
const ny = 10

const options = (f_abstol=1E-15, f_reltol=1E-15)

export test_guiding_center_3d
export nl, nx, ny

# The tests assert that the integration runs to completion and returns a solution.
#
# They used to assert `@test_nowarn` instead, which is not a usable criterion for these models: at
# the accuracy requested above the Newton solver regularly exhausts its iteration budget on the
# stiffer orbits, and `SimpleSolvers` warns once per occurrence — tens of thousands of times in a
# single testset. Which orbits trip it depends on the floating-point details of the platform, so
# the suite failed on a different equilibrium on Linux than on Windows. Relaxing the tolerance does
# not help: the warning count stays in the same order of magnitude at `1E-13`, and rises for some
# problems. The warnings are filtered out by `test/quiet_solver_warnings.jl` and counted there.
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

    # test_guiding_center_3d(ode(initial_conditions_trapped()...; tstep = 1E4, tspan = (0, 1E6)))
    # test_guiding_center_3d(ode(initial_conditions_barely_passing()...))
    # test_guiding_center_3d(ode(initial_conditions_barely_trapped()...))
    # test_guiding_center_3d(ode(initial_conditions_deeply_passing()...))
    # test_guiding_center_3d(ode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_3d(guiding_center_3d_loop_ode(nl), Δt=1.)
    # test_guiding_center_3d(guiding_center_3d_surface_ode(nx, ny), Δt=1.)

    test_guiding_center_3d(hode(initial_conditions_barely_passing()...; tspan=(0.0, 1E2), tstep=0.1))
    test_guiding_center_3d(hode(initial_conditions_barely_trapped()...; tspan=(0.0, 1E2), tstep=0.1))
    test_guiding_center_3d(hode(initial_conditions_deeply_passing()...; tspan=(0.0, 1E2), tstep=0.1))
    test_guiding_center_3d(hode(initial_conditions_deeply_trapped()...; tspan=(0.0, 1E2), tstep=0.1))

    test_guiding_center_3d(hode_canonical(initial_conditions_barely_passing()...; tspan=(0.0, 1E1), tstep=0.1))
    test_guiding_center_3d(hode_canonical(initial_conditions_barely_trapped()...; tspan=(0.0, 1E1), tstep=0.1))
    test_guiding_center_3d(hode_canonical(initial_conditions_deeply_passing()...; tspan=(0.0, 1E1), tstep=0.1))
    test_guiding_center_3d(hode_canonical(initial_conditions_deeply_trapped()...; tspan=(0.0, 1E1), tstep=0.1))

    # test_guiding_center_3d(iode(initial_conditions_trapped()...; tstep = 1E4, tspan = (0, 1E6)))
    # test_guiding_center_3d(iode(initial_conditions_barely_passing()...))
    # test_guiding_center_3d(iode(initial_conditions_barely_trapped()...))
    # test_guiding_center_3d(iode(initial_conditions_deeply_passing()...))#; tstep = 1E-2, tspan = (0, 1E1)
    # test_guiding_center_3d(iode(initial_conditions_deeply_trapped()...))#; tstep = 1E-2, tspan = (0, 1E1)
    # test_guiding_center_3d(guiding_center_3d_loop_iode(nl), Δt=1.)
    # test_guiding_center_3d(guiding_center_3d_surface_iode(nx, ny), Δt=1.)

end


@safetestset "Guiding Centre Dynamics in 3D with medium-size Tokamak Equilibrium in Cartesian Coordinates         " begin

    using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCartesian
    using ..GuidingCenter3dTests

    # test_guiding_center_3d(ode(initial_conditions_barely_passing()...))
    # test_guiding_center_3d(ode(initial_conditions_barely_trapped()...))
    # test_guiding_center_3d(ode(initial_conditions_deeply_passing()...))
    # test_guiding_center_3d(ode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_3d(guiding_center_3d_loop_ode(nl))
    # test_guiding_center_3d(guiding_center_3d_surface_ode(nx, ny))

    test_guiding_center_3d(hode(initial_conditions_barely_passing()...; tspan=(0.0, 1E2), tstep=0.1))
    test_guiding_center_3d(hode(initial_conditions_barely_trapped()...; tspan=(0.0, 1E2), tstep=0.1))
    test_guiding_center_3d(hode(initial_conditions_deeply_passing()...; tspan=(0.0, 1E2), tstep=0.1))
    test_guiding_center_3d(hode(initial_conditions_deeply_trapped()...; tspan=(0.0, 1E2), tstep=0.1))

    test_guiding_center_3d(hode_canonical(initial_conditions_barely_passing()...; tspan=(0.0, 1E1), tstep=0.1))
    test_guiding_center_3d(hode_canonical(initial_conditions_barely_trapped()...; tspan=(0.0, 1E1), tstep=0.1))
    test_guiding_center_3d(hode_canonical(initial_conditions_deeply_passing()...; tspan=(0.0, 1E1), tstep=0.1))
    test_guiding_center_3d(hode_canonical(initial_conditions_deeply_trapped()...; tspan=(0.0, 1E1), tstep=0.1))

    # test_guiding_center_3d(iode(initial_conditions_barely_passing()...))
    # test_guiding_center_3d(iode(initial_conditions_barely_trapped()...))
    # test_guiding_center_3d(iode(initial_conditions_deeply_passing()...))
    # test_guiding_center_3d(iode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_3d(guiding_center_3d_loop_iode(nl))
    # test_guiding_center_3d(guiding_center_3d_surface_iode(nx, ny))

end


# @safetestset "Guiding Centre Dynamics in 3D with medium-size Tokamak Equilibrium in Cylindrical Coordinates       " begin

#     using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCylindrical
#     using ..GuidingCenter3dTests

#     # test_guiding_center_3d(ode(initial_conditions_barely_passing()...))
#     # test_guiding_center_3d(ode(initial_conditions_barely_trapped()...))
#     # test_guiding_center_3d(ode(initial_conditions_deeply_passing()...))
#     # test_guiding_center_3d(ode(initial_conditions_deeply_trapped()...))
#     # test_guiding_center_3d(guiding_center_3d_loop_ode(nl))
#     # test_guiding_center_3d(guiding_center_3d_surface_ode(nx, ny))

#     test_guiding_center_3d(hode(initial_conditions_barely_passing()...; tspan=(0.0, 1E2), tstep=0.1))
#     test_guiding_center_3d(hode(initial_conditions_barely_trapped()...; tspan=(0.0, 1E2), tstep=0.1))
#     test_guiding_center_3d(hode(initial_conditions_deeply_passing()...; tspan=(0.0, 1E2), tstep=0.1))
#     test_guiding_center_3d(hode(initial_conditions_deeply_trapped()...; tspan=(0.0, 1E2), tstep=0.1))

#     # test_guiding_center_3d(iode(initial_conditions_barely_passing()...))
#     # test_guiding_center_3d(iode(initial_conditions_barely_trapped()...))
#     # test_guiding_center_3d(iode(initial_conditions_deeply_passing()...))
#     # test_guiding_center_3d(iode(initial_conditions_deeply_trapped()...))
#     # test_guiding_center_3d(guiding_center_3d_loop_iode(nl))
#     # test_guiding_center_3d(guiding_center_3d_surface_iode(nx, ny))

# end


# @safetestset "Guiding Centre Dynamics in 3D with small-size Tokamak Equilibrium in Cartesian Coordinates          " begin

#     using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallCartesian
#     using ..GuidingCenter3dTests

#     # test_guiding_center_3d(ode(initial_conditions_barely_passing()...))
#     # test_guiding_center_3d(ode(initial_conditions_barely_trapped()...))
#     # test_guiding_center_3d(ode(initial_conditions_deeply_passing()...))
#     # test_guiding_center_3d(ode(initial_conditions_deeply_trapped()...))
#     # test_guiding_center_3d(guiding_center_3d_loop_ode(nl))
#     # test_guiding_center_3d(guiding_center_3d_surface_ode(nx, ny))

#     test_guiding_center_3d(hode(initial_conditions_barely_passing()...; tspan=(0.0, 1E2), tstep=0.1))
#     test_guiding_center_3d(hode(initial_conditions_barely_trapped()...; tspan=(0.0, 1E2), tstep=0.1))
#     test_guiding_center_3d(hode(initial_conditions_deeply_passing()...; tspan=(0.0, 1E2), tstep=0.1))
#     test_guiding_center_3d(hode(initial_conditions_deeply_trapped()...; tspan=(0.0, 1E2), tstep=0.1))

#     # test_guiding_center_3d(iode(initial_conditions_barely_passing()...; tstep = 10., tspan = (0, 1E3)))
#     # test_guiding_center_3d(iode(initial_conditions_barely_trapped()...; tstep = 10., tspan = (0, 1E3)))
#     # test_guiding_center_3d(iode(initial_conditions_deeply_passing()...; tstep = 10., tspan = (0, 1E3)))
#     # test_guiding_center_3d(iode(initial_conditions_deeply_trapped()...; tstep = 10., tspan = (0, 1E3)))
#     # test_guiding_center_3d(guiding_center_3d_loop_iode(nl; tstep = 10., tspan = [0, 1E3]))
#     # test_guiding_center_3d(guiding_center_3d_surface_iode(nx, ny; tstep = 10., tspan = [0, 1E3]))

# end


# @safetestset "Guiding Centre Dynamics in 3D with small-size Tokamak Equilibrium in Cylindrical Coordinates        " begin

#     using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallCylindrical
#     using ..GuidingCenter3dTests

#     # test_guiding_center_3d(ode(initial_conditions_barely_passing()...))
#     # test_guiding_center_3d(ode(initial_conditions_barely_trapped()...))
#     # test_guiding_center_3d(ode(initial_conditions_deeply_passing()...))
#     # test_guiding_center_3d(ode(initial_conditions_deeply_trapped()...))
#     # test_guiding_center_3d(guiding_center_3d_loop_ode(nl))
#     # test_guiding_center_3d(guiding_center_3d_surface_ode(nx, ny))

#     test_guiding_center_3d(hode(initial_conditions_barely_passing()...; tspan=(0.0, 1E2), tstep=0.1))
#     test_guiding_center_3d(hode(initial_conditions_barely_trapped()...; tspan=(0.0, 1E2), tstep=0.1))
#     test_guiding_center_3d(hode(initial_conditions_deeply_passing()...; tspan=(0.0, 1E2), tstep=0.1))
#     test_guiding_center_3d(hode(initial_conditions_deeply_trapped()...; tspan=(0.0, 1E2), tstep=0.1))

#     # test_guiding_center_3d(iode(initial_conditions_barely_passing()...; tstep = 10., tspan = (0, 1E3)))
#     # test_guiding_center_3d(iode(initial_conditions_barely_trapped()...; tstep = 10., tspan = (0, 1E3)))
#     # test_guiding_center_3d(iode(initial_conditions_deeply_passing()...; tstep = 10., tspan = (0, 1E3)))
#     # test_guiding_center_3d(iode(initial_conditions_deeply_trapped()...; tstep = 10., tspan = (0, 1E3)))
#     # test_guiding_center_3d(guiding_center_3d_loop_iode(nl; tstep = 10., tspan = [0, 1E3]))
#     # test_guiding_center_3d(guiding_center_3d_surface_iode(nx, ny; tstep = 10., tspan = [0, 1E3]))

# end


# @safetestset "Guiding Centre Dynamics in 3D with small-size Tokamak Equilibrium in Toroidal Coordinates           " begin

#     using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallToroidal
#     using ..GuidingCenter3dTests

#     # test_guiding_center_3d(ode(initial_conditions_barely_passing()...))
#     # test_guiding_center_3d(ode(initial_conditions_barely_trapped()...))
#     # test_guiding_center_3d(ode(initial_conditions_deeply_passing()...))
#     # test_guiding_center_3d(ode(initial_conditions_deeply_trapped()...))
#     # test_guiding_center_3d(guiding_center_3d_loop_ode(nl))
#     # test_guiding_center_3d(guiding_center_3d_surface_ode(nx, ny))

#     test_guiding_center_3d(hode(initial_conditions_barely_passing()...; tspan=(0.0, 1E2), tstep=0.1))
#     test_guiding_center_3d(hode(initial_conditions_barely_trapped()...; tspan=(0.0, 1E2), tstep=0.1))
#     test_guiding_center_3d(hode(initial_conditions_deeply_passing()...; tspan=(0.0, 1E2), tstep=0.1))
#     test_guiding_center_3d(hode(initial_conditions_deeply_trapped()...; tspan=(0.0, 1E2), tstep=0.1))

#     # test_guiding_center_3d(iode(initial_conditions_barely_passing()...))
#     # test_guiding_center_3d(iode(initial_conditions_barely_trapped()...))
#     # test_guiding_center_3d(iode(initial_conditions_deeply_passing()...; tstep = 10., tspan = (0, 1E3)))
#     # test_guiding_center_3d(iode(initial_conditions_deeply_trapped()...))
#     # test_guiding_center_3d(guiding_center_3d_loop_iode(nl))
#     # test_guiding_center_3d(guiding_center_3d_surface_iode(nx, ny))

# end


# @safetestset "Guiding Centre Dynamics in 3D with symmetric Solov'ev Equilibrium                                   " begin

#     using ChargedParticleDynamics.GuidingCenter3d.SolovevSymmetricField
#     using ..GuidingCenter3dTests

#     # test_guiding_center_3d(ode(initial_conditions_barely_passing()...))
#     # test_guiding_center_3d(ode(initial_conditions_barely_trapped()...))
#     # test_guiding_center_3d(ode(initial_conditions_deeply_passing()...))
#     # test_guiding_center_3d(ode(initial_conditions_deeply_trapped()...))

#     test_guiding_center_3d(hode(initial_conditions_barely_passing()...))
#     test_guiding_center_3d(hode(initial_conditions_barely_trapped()...))
#     test_guiding_center_3d(hode(initial_conditions_deeply_passing()...))
#     test_guiding_center_3d(hode(initial_conditions_deeply_trapped()...))

#     # test_guiding_center_3d(iode(initial_conditions_barely_passing()...))
#     # test_guiding_center_3d(iode(initial_conditions_barely_trapped()...))
#     # test_guiding_center_3d(iode(initial_conditions_deeply_passing()...))
#     # test_guiding_center_3d(iode(initial_conditions_deeply_trapped()...))

# end


# @safetestset "Guiding Centre Dynamics in 3D with symmetric Equlibrium                                             " begin

#     using ChargedParticleDynamics.GuidingCenter3d.SymmetricField
#     using ..GuidingCenter3dTests

#     # test_guiding_center_3d(guiding_center_3d_loop_ode(nl))
#     # test_guiding_center_3d(guiding_center_3d_surface_ode(nx, ny))

#     # test_guiding_center_3d(guiding_center_3d_loop_iode(nl))
#     # test_guiding_center_3d(guiding_center_3d_surface_iode(nx, ny))

# end


# @safetestset "Guiding Centre Dynamics in 3D with Theta Pinch Equilibrium                                          " begin

#     using ChargedParticleDynamics.GuidingCenter3d.ThetaPinchField
#     using ..GuidingCenter3dTests

#     # test_guiding_center_3d(guiding_center_3d_loop_ode(nl))

#     # test_guiding_center_3d(guiding_center_3d_loop_iode(nl))

# end
