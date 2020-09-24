
using SafeTestsets

module GuidingCenter4dTests

    using Test
    using GeometricIntegrators

    const Δt = 1.
    const nt = 1

    const nl = 100
    const nx = 10
    const ny = 10

    export test_guiding_center_4d_glrk, test_guiding_center_4d_vpglrk
    export nl, nx, ny

    function test_guiding_center_4d_glrk(ode; Δt = Δt)
        int = Integrator(ode, getTableauGLRK(1), Δt)
        sol = integrate(ode, int, nt)
        @test true
    end

    function test_guiding_center_4d_vpglrk(ode; Δt = Δt)
        int = Integrator(ode, getTableauVPGLRK(1), Δt)
        sol = integrate(ode, int, nt)
        @test true
    end

end



@safetestset "Guiding Centre Dynamics in 4D with ITER-like Solov'ev Equilibrium with X-Point                      " begin

    using ChargedParticleDynamics.GuidingCenter4d.SolovevIterXpoint
    using ..GuidingCenter4dTests

    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_trapped()...), Δt=10000.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_passing()...), Δt=1.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_trapped()...), Δt=1.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_passing()...), Δt=1.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_trapped()...), Δt=1.)
    # test_guiding_center_4d_glrk(guiding_center_4d_loop_ode(nl), Δt=1.)
    # test_guiding_center_4d_glrk(guiding_center_4d_surface_ode(nx, ny), Δt=1.)

    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_trapped()...), Δt=10000.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_passing()...), Δt=1.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_trapped()...), Δt=1.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_passing()...), Δt=1.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_trapped()...), Δt=1.)
    # test_guiding_center_4d_vpglrk(guiding_center_4d_loop_iode(nl), Δt=1.)
    # test_guiding_center_4d_vpglrk(guiding_center_4d_surface_iode(nx, ny), Δt=1.)

end


@safetestset "Guiding Centre Dynamics in 4D with medium-size Tokamak Equilibrium in Cartesian Coordinates         " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCartesian
    using ..GuidingCenter4dTests

    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_passing()...), Δt=1.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_trapped()...), Δt=1.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_passing()...), Δt=1.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_trapped()...), Δt=1.)
    test_guiding_center_4d_glrk(guiding_center_4d_loop_ode(nl), Δt=1.)
    test_guiding_center_4d_glrk(guiding_center_4d_surface_ode(nx, ny), Δt=1.)

    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_passing()...), Δt=1.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_trapped()...), Δt=1.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_passing()...), Δt=1.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_trapped()...), Δt=1.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_loop_iode(nl), Δt=.1)
    test_guiding_center_4d_vpglrk(guiding_center_4d_surface_iode(nx, ny), Δt=.1)

end


@safetestset "Guiding Centre Dynamics in 4D with medium-size Tokamak Equilibrium in Cylindrical Coordinates       " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCylindrical
    using ..GuidingCenter4dTests

    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_passing()...), Δt=1.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_trapped()...), Δt=1.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_passing()...), Δt=1.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_trapped()...), Δt=1.)
    test_guiding_center_4d_glrk(guiding_center_4d_loop_ode(nl), Δt=1.)
    test_guiding_center_4d_glrk(guiding_center_4d_surface_ode(nx, ny), Δt=1.)

    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_passing()...), Δt=1.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_trapped()...), Δt=1.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_passing()...), Δt=1.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_trapped()...), Δt=1.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_loop_iode(nl), Δt=.1)
    test_guiding_center_4d_vpglrk(guiding_center_4d_surface_iode(nx, ny), Δt=.1)

end


@safetestset "Guiding Centre Dynamics in 4D with small-size Tokamak Equilibrium in Cartesian Coordinates          " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCartesian
    using ..GuidingCenter4dTests

    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_passing()...), Δt=400.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_trapped()...), Δt=400.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_passing()...), Δt=400.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_trapped()...), Δt=400.)
    test_guiding_center_4d_glrk(guiding_center_4d_loop_ode(nl), Δt=400.)
    test_guiding_center_4d_glrk(guiding_center_4d_surface_ode(nx, ny), Δt=400.)

    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_passing()...), Δt=10.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_trapped()...), Δt=10.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_passing()...), Δt=10.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_trapped()...), Δt=10.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_loop_iode(nl), Δt=10.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_surface_iode(nx, ny), Δt=10.)

end


@safetestset "Guiding Centre Dynamics in 4D with small-size Tokamak Equilibrium in Circular Coordinates           " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCircular
    using ..GuidingCenter4dTests

    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_passing()...), Δt=400.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_trapped()...), Δt=400.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_passing()...), Δt=400.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_trapped()...), Δt=400.)
    test_guiding_center_4d_glrk(guiding_center_4d_loop_ode(nl), Δt=400.)
    test_guiding_center_4d_glrk(guiding_center_4d_surface_ode(nx, ny), Δt=400.)

    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_passing()...), Δt=400.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_trapped()...), Δt=400.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_passing()...), Δt=400.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_trapped()...), Δt=400.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_loop_iode(nl), Δt=400.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_surface_iode(nx, ny), Δt=400.)

end


@safetestset "Guiding Centre Dynamics in 4D with small-size Tokamak Equilibrium in Cylindrical Coordinates        " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCylindrical
    using ..GuidingCenter4dTests

    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_passing()...), Δt=400.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_trapped()...), Δt=400.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_passing()...), Δt=400.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_trapped()...), Δt=400.)
    test_guiding_center_4d_glrk(guiding_center_4d_loop_ode(nl), Δt=400.)
    test_guiding_center_4d_glrk(guiding_center_4d_surface_ode(nx, ny), Δt=400.)

    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_passing()...), Δt=10.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_trapped()...), Δt=10.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_passing()...), Δt=10.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_trapped()...), Δt=10.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_loop_iode(nl), Δt=10.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_surface_iode(nx, ny), Δt=10.)

end


@safetestset "Guiding Centre Dynamics in 4D with quadratic Solov'ev Equilibrium                                   " begin

    using ChargedParticleDynamics.GuidingCenter4d.SolovevQuadraticField
    using ..GuidingCenter4dTests

    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_passing()...), Δt=1.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_trapped()...), Δt=1.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_passing()...), Δt=1.)
    test_guiding_center_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_trapped()...), Δt=1.)

    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_passing()...), Δt=1.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_barely_trapped()...), Δt=1.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_passing()...), Δt=1.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_iode(initial_conditions_deeply_trapped()...), Δt=1.)

end


@safetestset "Guiding Centre Dynamics in 4D with symmetric quadratic Equlibrium                                   " begin

    using ChargedParticleDynamics.GuidingCenter4d.SymmetricField
    using ..GuidingCenter4dTests

    test_guiding_center_4d_glrk(guiding_center_4d_loop_ode(nl), Δt=1.)
    test_guiding_center_4d_glrk(guiding_center_4d_surface_ode(nx, ny), Δt=1.)

    test_guiding_center_4d_vpglrk(guiding_center_4d_loop_iode(nl), Δt=1.)
    test_guiding_center_4d_vpglrk(guiding_center_4d_surface_iode(nx, ny), Δt=1.)

end


@safetestset "Guiding Centre Dynamics in 4D with Theta Pinch Equilibrium                                          " begin

    using ChargedParticleDynamics.GuidingCenter4d.ThetaPinchField
    using ..GuidingCenter4dTests

    test_guiding_center_4d_glrk(guiding_center_4d_loop_ode(nl), Δt=1.)

    test_guiding_center_4d_vpglrk(guiding_center_4d_loop_iode(nl), Δt=1.)

end
