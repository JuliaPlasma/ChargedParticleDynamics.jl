
using SafeTestsets

module GyroKinetics4dTests

    using Test
    using GeometricIntegrators
    using ChargedParticleDynamics.GyroKinetics4d

    const nl = 100
    const nx = 10
    const ny = 10

    export test_gyro_kinetics_4d_erk4, test_gyro_kinetics_4d_glrk, test_gyro_kinetics_4d_strang
    export nl, nx, ny

    function test_gyro_kinetics_4d_erk4(ode)
        sol = integrate(ode, RK4())
        @test true
    end

    function test_gyro_kinetics_4d_glrk(ode)
        sol = integrate(ode, Gauss(1))
        @test true
    end

    function test_gyro_kinetics_4d_strang(ode::SODEProblem{DT}) where {DT}
        mpi = Tuple(Gauss(1) for _ in 1:6)
        # mpi = Tuple((v::Function, Δt::Number; kwargs...) -> Integrator{DT, ndims(ode)}(v, Gauss(1), Δt; kwargs...) for i in 1:6)
        # mpi = Tuple((v::Function, Δt::Number; kwargs...) -> IntegratorFIRKwCT{DT, ndims(ode)}(v, ωabs, ode.parameters, Gauss(1), Δt; kwargs...) for i in 1:6)
        sol = integrate(ode, Composition(mpi, Strang()))
        @test true
    end

end


@safetestset "Gyrokinetic GC Model in 4D with ITER-like Solov'ev Equilibrium with X-Point                         " begin

    using SimpleSolvers
    using ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dSolovevIterXpoint
    # using ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dSolovevIterXpoint: scaling_factor, transform_q̃_to_q_jacobian!, transform_q̃_to_q_rhs!
    using ..GyroKinetics4dTests

    # test coordinate transformation
    q̃, params = initial_conditions_trapped()
    # q̃[2] = 1.0
    # q̃[3] = π/4
    t  = 0.
    # q  = q̃ ./ ωabs(t, [params[:x₀]..., q̃[4]], params) ./ scaling_factor(params)
    # j1 = zeros(eltype(q), length(q), length(q))
    # j2 = zeros(eltype(q), length(q), length(q))

    # transform_q̃_to_q_jacobian!(q, j1, t, params)
    # computeJacobianAD(q, j2, (q,b) -> transform_q̃_to_q_rhs!(q, b, t, q̃, params))
    # @test j1 ≈ j2 atol=1e-5

    @test transform_q_to_q̃(t, transform_q̃_to_q(t, q̃, params), params) ≈ q̃ atol=1E-6

    # test ODE
    test_gyro_kinetics_4d_erk4(guiding_center_4d_ode(initial_conditions_trapped()...; tstep = 1.))
    test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_trapped()...; tstep = 1.))
    # test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_passing()...), Δt=1.)
    # test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_trapped()...), Δt=1.)
    # test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_passing()...), Δt=1.)
    # test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_trapped()...), Δt=1.)

    # test SODE
    test_gyro_kinetics_4d_strang(guiding_center_4d_sode(initial_conditions_trapped()...; tstep = 1.))
    # test_gyro_kinetics_4d_strang(guiding_center_4d_sode(initial_conditions_barely_passing()...), Δt=1.)
    # test_gyro_kinetics_4d_strang(guiding_center_4d_sode(initial_conditions_barely_trapped()...), Δt=1.)
    # test_gyro_kinetics_4d_strang(guiding_center_4d_sode(initial_conditions_deeply_passing()...), Δt=1.)
    # test_gyro_kinetics_4d_strang(guiding_center_4d_sode(initial_conditions_deeply_trapped()...), Δt=1.)

end
