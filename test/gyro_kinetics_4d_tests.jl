
using SafeTestsets

module GyroKinetics4dTests

    using Test
    using GeometricIntegrators

    set_config(:nls_stol_break, Inf)

    const Δt = 1.
    const nt = 1

    const nl = 100
    const nx = 10
    const ny = 10

    export test_gyro_kinetics_4d_glrk, test_gyro_kinetics_4d_strang
    export nl, nx, ny

    function test_gyro_kinetics_4d_glrk(ode; Δt = Δt)
        int = Integrator(ode, getTableauGLRK(1), Δt)
        sol = integrate(ode, int, nt)
        @test true
    end

    function test_gyro_kinetics_4d_strang(ode; Δt = Δt)
        tab = getTableauStrang()
        mpi = Tuple(IntegratorConstructor(eltype(ode.q₀), ndims(ode), getTableauGLRK(1)) for i in 1:6)
        int = IntegratorComposition(ode, mpi, tab, Δt)
        sol = integrate(ode, int, nt)
        @test true
    end

end



@safetestset "Gyrokinetic GC Model in 4D with ITER-like Solov'ev Equilibrium with X-Point                         " begin

    using ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dSolovevIterXpoint
    using ..GyroKinetics4dTests

    test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_trapped()...), Δt=100.)
    # test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_passing()...), Δt=1.)
    # test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_trapped()...), Δt=1.)
    # test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_passing()...), Δt=1.)
    # test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_trapped()...), Δt=1.)

    test_gyro_kinetics_4d_strang(guiding_center_4d_sode(initial_conditions_trapped()...), Δt=100.)
    # test_gyro_kinetics_4d_strang(guiding_center_4d_sode(initial_conditions_barely_passing()...), Δt=1.)
    # test_gyro_kinetics_4d_strang(guiding_center_4d_sode(initial_conditions_barely_trapped()...), Δt=1.)
    # test_gyro_kinetics_4d_strang(guiding_center_4d_sode(initial_conditions_deeply_passing()...), Δt=1.)
    # test_gyro_kinetics_4d_strang(guiding_center_4d_sode(initial_conditions_deeply_trapped()...), Δt=1.)

end
