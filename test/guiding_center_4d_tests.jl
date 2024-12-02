
using SafeTestsets

module GuidingCenter4dTests

    using Test
    using GeometricIntegrators

    const nl = 100
    const nx = 10
    const ny = 10

    export test_guiding_center_4d_glrk, test_guiding_center_4d_vpglrk
    export nl, nx, ny

    function test_guiding_center_4d_glrk(ode)
        sol = integrate(ode, Gauss(2))
        @test true
    end

    function test_guiding_center_4d_vpglrk(ode)
        sol = integrate(ode, VPRKGauss(2))
        @test true
    end

end



@safetestset "Guiding Centre Dynamics in 4D with ITER-like Solov'ev Equilibrium with X-Point                      " begin

    using ChargedParticleDynamics.GuidingCenter4d.SolovevIterXpoint
    using ..GuidingCenter4dTests

    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_trapped()...; tstep = 1E4, tspan = (0, 1E6)))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_passing()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_trapped()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_passing()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_4d_glrk(guiding_center_4d_loop_ode(nl), Δt=1.)
    # test_guiding_center_4d_glrk(guiding_center_4d_surface_ode(nx, ny), Δt=1.)

    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_trapped()...; tstep = 1E4, tspan = (0, 1E6)))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_passing()...))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_trapped()...))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_passing()...; tstep = 1E-2, tspan = (0, 1E1)))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_trapped()...; tstep = 1E-2, tspan = (0, 1E1)))
    # test_guiding_center_4d_vpglrk(guiding_center_4d_loop_iode(nl), Δt=1.)
    # test_guiding_center_4d_vpglrk(guiding_center_4d_surface_iode(nx, ny), Δt=1.)

end


@safetestset "Guiding Centre Dynamics in 4D with medium-size Tokamak Equilibrium in Cartesian Coordinates         " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCartesian
    using ..GuidingCenter4dTests

    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_passing()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_trapped()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_passing()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_4d_glrk(guiding_center_4d_loop_ode(nl))
    # test_guiding_center_4d_glrk(guiding_center_4d_surface_ode(nx, ny))

    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_passing()...))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_trapped()...))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_passing()...))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_4d_vpglrk(guiding_center_4d_loop_iode(nl))
    # test_guiding_center_4d_vpglrk(guiding_center_4d_surface_iode(nx, ny))

end


@safetestset "Guiding Centre Dynamics in 4D with medium-size Tokamak Equilibrium in Cylindrical Coordinates       " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCylindrical
    using ..GuidingCenter4dTests

    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_passing()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_trapped()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_passing()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_4d_glrk(guiding_center_4d_loop_ode(nl))
    # test_guiding_center_4d_glrk(guiding_center_4d_surface_ode(nx, ny))

    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_passing()...))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_trapped()...))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_passing()...))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_4d_vpglrk(guiding_center_4d_loop_iode(nl))
    # test_guiding_center_4d_vpglrk(guiding_center_4d_surface_iode(nx, ny))

end


@safetestset "Guiding Centre Dynamics in 4D with small-size Tokamak Equilibrium in Cartesian Coordinates          " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCartesian
    using ..GuidingCenter4dTests

    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_passing()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_trapped()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_passing()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_4d_glrk(guiding_center_4d_loop_ode(nl))
    # test_guiding_center_4d_glrk(guiding_center_4d_surface_ode(nx, ny))

    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_passing()...; tstep = 10., tspan = (0, 1E3)))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_trapped()...; tstep = 10., tspan = (0, 1E3)))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_passing()...; tstep = 10., tspan = (0, 1E3)))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_trapped()...; tstep = 10., tspan = (0, 1E3)))
    # test_guiding_center_4d_vpglrk(guiding_center_4d_loop_iode(nl; tstep = 10., tspan = [0, 1E3]))
    # test_guiding_center_4d_vpglrk(guiding_center_4d_surface_iode(nx, ny; tstep = 10., tspan = [0, 1E3]))

end


@safetestset "Guiding Centre Dynamics in 4D with small-size Tokamak Equilibrium in Cylindrical Coordinates        " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCylindrical
    using ..GuidingCenter4dTests

    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_passing()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_trapped()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_passing()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_4d_glrk(guiding_center_4d_loop_ode(nl))
    # test_guiding_center_4d_glrk(guiding_center_4d_surface_ode(nx, ny))

    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_passing()...; tstep = 10., tspan = (0, 1E3)))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_trapped()...; tstep = 10., tspan = (0, 1E3)))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_passing()...; tstep = 10., tspan = (0, 1E3)))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_trapped()...; tstep = 10., tspan = (0, 1E3)))
    # test_guiding_center_4d_vpglrk(guiding_center_4d_loop_iode(nl; tstep = 10., tspan = [0, 1E3]))
    # test_guiding_center_4d_vpglrk(guiding_center_4d_surface_iode(nx, ny; tstep = 10., tspan = [0, 1E3]))

end


@safetestset "Guiding Centre Dynamics in 4D with small-size Tokamak Equilibrium in Toroidal Coordinates           " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakSmallToroidal
    using ..GuidingCenter4dTests

    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_passing()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_trapped()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_passing()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_4d_glrk(guiding_center_4d_loop_ode(nl))
    # test_guiding_center_4d_glrk(guiding_center_4d_surface_ode(nx, ny))

    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_passing()...))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_trapped()...))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_passing()...; tstep = 10., tspan = (0, 1E3)))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_trapped()...))
    # test_guiding_center_4d_vpglrk(guiding_center_4d_loop_iode(nl))
    # test_guiding_center_4d_vpglrk(guiding_center_4d_surface_iode(nx, ny))

end


@safetestset "Guiding Centre Dynamics in 4D with symmetric Solov'ev Equilibrium                                   " begin

    using ChargedParticleDynamics.GuidingCenter4d.SolovevSymmetricField
    using ..GuidingCenter4dTests

    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_passing()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_trapped()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_passing()...))
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_trapped()...))

    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_passing()...))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_trapped()...))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_passing()...))
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_trapped()...))

end


@safetestset "Guiding Centre Dynamics in 4D with symmetric Equlibrium                                             " begin

    using ChargedParticleDynamics.GuidingCenter4d.SymmetricField
    using ..GuidingCenter4dTests

    # test_guiding_center_4d_glrk(guiding_center_4d_loop_ode(nl))
    # test_guiding_center_4d_glrk(guiding_center_4d_surface_ode(nx, ny))

    # test_guiding_center_4d_vpglrk(guiding_center_4d_loop_iode(nl))
    # test_guiding_center_4d_vpglrk(guiding_center_4d_surface_iode(nx, ny))

end


@safetestset "Guiding Centre Dynamics in 4D with Theta Pinch Equilibrium                                          " begin

    using ChargedParticleDynamics.GuidingCenter4d.ThetaPinchField
    using ..GuidingCenter4dTests

    # test_guiding_center_4d_glrk(guiding_center_4d_loop_ode(nl))

    # test_guiding_center_4d_vpglrk(guiding_center_4d_loop_iode(nl))

end
