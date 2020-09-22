
using SafeTestsets

module GyroKinetics4dTests

    using Test
    using GeometricIntegrators
    using ChargedParticleDynamics.GyroKinetics4d

    const Δt = 1.
    const nt = 1

    const nl = 100
    const nx = 10
    const ny = 10

    export test_gyro_kinetics_4d_erk4, test_gyro_kinetics_4d_glrk, test_gyro_kinetics_4d_strang
    export nl, nx, ny

    function test_gyro_kinetics_4d_erk4(ode, ωabs; Δt = Δt)
        int = Integrator(ode, getTableauERK4(), Δt)
        sol = integrate(ode, int, nt)
        @test true
    end

    function test_gyro_kinetics_4d_glrk(ode, ωabs; Δt = Δt)
        int = Integrator(ode, getTableauGLRK(1), Δt)
        sol = integrate(ode, int, nt)
        @test true
    end

    function test_gyro_kinetics_4d_strang(ode, ωabs; Δt = Δt)
        tab = getTableauStrang()
        mpi = Tuple(IntegratorConstructor(eltype(ode.q₀), ndims(ode), getTableauGLRK(1)) for i in 1:6)
        int = IntegratorComposition(ode, mpi, tab, Δt)
        sol = integrate(ode, int, nt)
        @test true
    end

end


@safetestset "Gyrokinetic GC Model in 4D with ITER-like Solov'ev Equilibrium with X-Point                         " begin

    using GeometricIntegrators.Solvers
    using ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dSolovevIterXpoint
    # using ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dSolovevIterXpoint: scaling_factor, transform_q̃_to_q_jacobian!, transform_q̃_to_q_rhs!
    using ..GyroKinetics4dTests

    # test coordinate transformation
    q̃, params = initial_conditions_trapped()
    t  = 0.
    @test transform_q_to_q̃(t, transform_q̃_to_q(t, q̃, params), params) ≈ q̃ atol=1E-6

    # test ODE
    test_gyro_kinetics_4d_erk4(guiding_center_4d_ode(initial_conditions_trapped()...), ωabs; Δt=1.)
    test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_trapped()...), ωabs; Δt=1.)
    # test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_passing()...), Δt=1.)
    # test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_trapped()...), Δt=1.)
    # test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_passing()...), Δt=1.)
    # test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_trapped()...), Δt=1.)

    # test SODE
    test_gyro_kinetics_4d_strang(guiding_center_4d_sode(initial_conditions_trapped()...), ωabs; Δt=1.)
    # test_gyro_kinetics_4d_strang(guiding_center_4d_sode(initial_conditions_barely_passing()...), Δt=1.)
    # test_gyro_kinetics_4d_strang(guiding_center_4d_sode(initial_conditions_barely_trapped()...), Δt=1.)
    # test_gyro_kinetics_4d_strang(guiding_center_4d_sode(initial_conditions_deeply_passing()...), Δt=1.)
    # test_gyro_kinetics_4d_strang(guiding_center_4d_sode(initial_conditions_deeply_trapped()...), Δt=1.)

end
