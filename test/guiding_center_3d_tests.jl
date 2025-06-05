
using SafeTestsets

module GuidingCenter3dTests

using GeometricIntegrators
using SimpleSolvers: Options
using Test

const nl = 100
const nx = 10
const ny = 10

const options = Options(x_reltol=1E-14, f_abstol=1E-14, f_reltol=1E-14)

export test_guiding_center_3d_glrk, test_guiding_center_3d_pglrk, test_guiding_center_3d_vpglrk
export nl, nx, ny

function test_guiding_center_3d_glrk(equ)
    @test_nowarn integrate(equ, Gauss(2); options=options)
end

function test_guiding_center_3d_pglrk(equ)
    @test_nowarn integrate(equ, PartitionedGauss(2); options=options)
end

function test_guiding_center_3d_vpglrk(equ)
    @test_nowarn integrate(equ, VPRKGauss(2); options=options)
end

end



@safetestset "Guiding Centre Dynamics in 3D with ITER-like Solov'ev Equilibrium with X-Point                      " begin

    using ChargedParticleDynamics.GuidingCenter3d.SolovevIterXpoint
    using ..GuidingCenter3dTests

    # test_guiding_center_3d_glrk(ode(initial_conditions_trapped()...; tstep = 1E4, tspan = (0, 1E6)))
    # test_guiding_center_3d_glrk(ode(initial_conditions_barely_passing()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_barely_trapped()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_deeply_passing()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_3d_glrk(guiding_center_3d_loop_ode(nl), Δt=1.)
    # test_guiding_center_3d_glrk(guiding_center_3d_surface_ode(nx, ny), Δt=1.)

    # test_guiding_center_3d_pglrk(hode(initial_conditions_barely_passing()...))
    # test_guiding_center_3d_pglrk(hode(initial_conditions_barely_trapped()...))
    # test_guiding_center_3d_pglrk(hode(initial_conditions_deeply_passing()...))
    # test_guiding_center_3d_pglrk(hode(initial_conditions_deeply_trapped()...))

    # test_guiding_center_3d_vpglrk(iode(initial_conditions_trapped()...; tstep = 1E4, tspan = (0, 1E6)))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_barely_passing()...))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_barely_trapped()...))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_deeply_passing()...; tstep = 1E-2, tspan = (0, 1E1)))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_deeply_trapped()...; tstep = 1E-2, tspan = (0, 1E1)))
    # test_guiding_center_3d_vpglrk(guiding_center_3d_loop_iode(nl), Δt=1.)
    # test_guiding_center_3d_vpglrk(guiding_center_3d_surface_iode(nx, ny), Δt=1.)

end


@safetestset "Guiding Centre Dynamics in 3D with medium-size Tokamak Equilibrium in Cartesian Coordinates         " begin

    using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCartesian
    using ..GuidingCenter3dTests

    # test_guiding_center_3d_glrk(ode(initial_conditions_barely_passing()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_barely_trapped()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_deeply_passing()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_3d_glrk(guiding_center_3d_loop_ode(nl))
    # test_guiding_center_3d_glrk(guiding_center_3d_surface_ode(nx, ny))

    test_guiding_center_3d_pglrk(hode(initial_conditions_barely_passing()...))
    test_guiding_center_3d_pglrk(hode(initial_conditions_barely_trapped()...))
    test_guiding_center_3d_pglrk(hode(initial_conditions_deeply_passing()...))
    test_guiding_center_3d_pglrk(hode(initial_conditions_deeply_trapped()...))

    # test_guiding_center_3d_vpglrk(iode(initial_conditions_barely_passing()...))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_barely_trapped()...))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_deeply_passing()...))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_3d_vpglrk(guiding_center_3d_loop_iode(nl))
    # test_guiding_center_3d_vpglrk(guiding_center_3d_surface_iode(nx, ny))

end


@safetestset "Guiding Centre Dynamics in 3D with medium-size Tokamak Equilibrium in Cylindrical Coordinates       " begin

    using ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCylindrical
    using ..GuidingCenter3dTests

    # test_guiding_center_3d_glrk(ode(initial_conditions_barely_passing()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_barely_trapped()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_deeply_passing()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_3d_glrk(guiding_center_3d_loop_ode(nl))
    # test_guiding_center_3d_glrk(guiding_center_3d_surface_ode(nx, ny))

    test_guiding_center_3d_pglrk(hode(initial_conditions_barely_passing()...))
    test_guiding_center_3d_pglrk(hode(initial_conditions_barely_trapped()...))
    test_guiding_center_3d_pglrk(hode(initial_conditions_deeply_passing()...))
    test_guiding_center_3d_pglrk(hode(initial_conditions_deeply_trapped()...))

    # test_guiding_center_3d_vpglrk(iode(initial_conditions_barely_passing()...))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_barely_trapped()...))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_deeply_passing()...))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_3d_vpglrk(guiding_center_3d_loop_iode(nl))
    # test_guiding_center_3d_vpglrk(guiding_center_3d_surface_iode(nx, ny))

end


@safetestset "Guiding Centre Dynamics in 3D with small-size Tokamak Equilibrium in Cartesian Coordinates          " begin

    using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallCartesian
    using ..GuidingCenter3dTests

    # test_guiding_center_3d_glrk(ode(initial_conditions_barely_passing()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_barely_trapped()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_deeply_passing()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_3d_glrk(guiding_center_3d_loop_ode(nl))
    # test_guiding_center_3d_glrk(guiding_center_3d_surface_ode(nx, ny))

    test_guiding_center_3d_pglrk(hode(initial_conditions_barely_passing()...))
    test_guiding_center_3d_pglrk(hode(initial_conditions_barely_trapped()...))
    test_guiding_center_3d_pglrk(hode(initial_conditions_deeply_passing()...))
    test_guiding_center_3d_pglrk(hode(initial_conditions_deeply_trapped()...))

    # test_guiding_center_3d_vpglrk(iode(initial_conditions_barely_passing()...; tstep = 10., tspan = (0, 1E3)))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_barely_trapped()...; tstep = 10., tspan = (0, 1E3)))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_deeply_passing()...; tstep = 10., tspan = (0, 1E3)))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_deeply_trapped()...; tstep = 10., tspan = (0, 1E3)))
    # test_guiding_center_3d_vpglrk(guiding_center_3d_loop_iode(nl; tstep = 10., tspan = [0, 1E3]))
    # test_guiding_center_3d_vpglrk(guiding_center_3d_surface_iode(nx, ny; tstep = 10., tspan = [0, 1E3]))

end


@safetestset "Guiding Centre Dynamics in 3D with small-size Tokamak Equilibrium in Cylindrical Coordinates        " begin

    using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallCylindrical
    using ..GuidingCenter3dTests

    # test_guiding_center_3d_glrk(ode(initial_conditions_barely_passing()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_barely_trapped()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_deeply_passing()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_3d_glrk(guiding_center_3d_loop_ode(nl))
    # test_guiding_center_3d_glrk(guiding_center_3d_surface_ode(nx, ny))

    test_guiding_center_3d_pglrk(hode(initial_conditions_barely_passing()...))
    test_guiding_center_3d_pglrk(hode(initial_conditions_barely_trapped()...))
    test_guiding_center_3d_pglrk(hode(initial_conditions_deeply_passing()...))
    test_guiding_center_3d_pglrk(hode(initial_conditions_deeply_trapped()...))

    # test_guiding_center_3d_vpglrk(iode(initial_conditions_barely_passing()...; tstep = 10., tspan = (0, 1E3)))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_barely_trapped()...; tstep = 10., tspan = (0, 1E3)))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_deeply_passing()...; tstep = 10., tspan = (0, 1E3)))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_deeply_trapped()...; tstep = 10., tspan = (0, 1E3)))
    # test_guiding_center_3d_vpglrk(guiding_center_3d_loop_iode(nl; tstep = 10., tspan = [0, 1E3]))
    # test_guiding_center_3d_vpglrk(guiding_center_3d_surface_iode(nx, ny; tstep = 10., tspan = [0, 1E3]))

end


@safetestset "Guiding Centre Dynamics in 3D with small-size Tokamak Equilibrium in Toroidal Coordinates           " begin

    using ChargedParticleDynamics.GuidingCenter3d.TokamakSmallToroidal
    using ..GuidingCenter3dTests

    # test_guiding_center_3d_glrk(ode(initial_conditions_barely_passing()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_barely_trapped()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_deeply_passing()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_3d_glrk(guiding_center_3d_loop_ode(nl))
    # test_guiding_center_3d_glrk(guiding_center_3d_surface_ode(nx, ny))

    test_guiding_center_3d_pglrk(hode(initial_conditions_barely_passing()...))
    test_guiding_center_3d_pglrk(hode(initial_conditions_barely_trapped()...))
    test_guiding_center_3d_pglrk(hode(initial_conditions_deeply_passing()...))
    test_guiding_center_3d_pglrk(hode(initial_conditions_deeply_trapped()...))

    # test_guiding_center_3d_vpglrk(iode(initial_conditions_barely_passing()...))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_barely_trapped()...))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_deeply_passing()...; tstep = 10., tspan = (0, 1E3)))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_3d_vpglrk(guiding_center_3d_loop_iode(nl))
    # test_guiding_center_3d_vpglrk(guiding_center_3d_surface_iode(nx, ny))

end


@safetestset "Guiding Centre Dynamics in 3D with symmetric Solov'ev Equilibrium                                   " begin

    using ChargedParticleDynamics.GuidingCenter3d.SolovevSymmetricField
    using ..GuidingCenter3dTests

    # test_guiding_center_3d_glrk(ode(initial_conditions_barely_passing()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_barely_trapped()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_deeply_passing()...))
    # test_guiding_center_3d_glrk(ode(initial_conditions_deeply_trapped()...))

    test_guiding_center_3d_pglrk(hode(initial_conditions_barely_passing()...))
    test_guiding_center_3d_pglrk(hode(initial_conditions_barely_trapped()...))
    test_guiding_center_3d_pglrk(hode(initial_conditions_deeply_passing()...))
    test_guiding_center_3d_pglrk(hode(initial_conditions_deeply_trapped()...))

    # test_guiding_center_3d_vpglrk(iode(initial_conditions_barely_passing()...))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_barely_trapped()...))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_deeply_passing()...))
    # test_guiding_center_3d_vpglrk(iode(initial_conditions_deeply_trapped()...))

end


@safetestset "Guiding Centre Dynamics in 3D with symmetric Equlibrium                                             " begin

    using ChargedParticleDynamics.GuidingCenter3d.SymmetricField
    using ..GuidingCenter3dTests

    # test_guiding_center_3d_glrk(guiding_center_3d_loop_ode(nl))
    # test_guiding_center_3d_glrk(guiding_center_3d_surface_ode(nx, ny))

    # test_guiding_center_3d_vpglrk(guiding_center_3d_loop_iode(nl))
    # test_guiding_center_3d_vpglrk(guiding_center_3d_surface_iode(nx, ny))

end


@safetestset "Guiding Centre Dynamics in 3D with Theta Pinch Equilibrium                                          " begin

    using ChargedParticleDynamics.GuidingCenter3d.ThetaPinchField
    using ..GuidingCenter3dTests

    # test_guiding_center_3d_glrk(guiding_center_3d_loop_ode(nl))

    # test_guiding_center_3d_vpglrk(guiding_center_3d_loop_iode(nl))

end
