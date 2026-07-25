
using SafeTestsets

module GyroKinetics4dTests

using ChargedParticleDynamics.GyroKinetics4d
using GeometricIntegrators
using Test

const options = (f_abstol=1E-15, f_reltol=1E-15)

export test_gyro_kinetics_4d_erk4, test_gyro_kinetics_4d_glrk, test_gyro_kinetics_4d_strang
export strang_composition

# Assert that the integration returns a solution rather than `@test_nowarn`; see the comment on the
# corresponding functions in `guiding_center_3d_tests.jl` for why.
function test_gyro_kinetics_4d_erk4(equ::ODEProblem)
    @test integrate(equ, RK4()) isa GeometricSolution
end

function test_gyro_kinetics_4d_glrk(equ::ODEProblem)
    @test integrate(equ, Gauss(1); options...) isa GeometricSolution
end

"Symmetric composition of symplectic integrators for the six subsystems of the splitting."
strang_composition() = Composition(Tuple(Gauss(1) for _ in 1:6), Strang())

function test_gyro_kinetics_4d_strang(ode::SODEProblem)
    @test integrate(ode, strang_composition()) isa GeometricSolution
end

end


@safetestset "Gyrokinetic GC Model in 4D with ITER-like Solov'ev Equilibrium with X-Point                         " begin

    using ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dSolovevIterXpoint
    using ..GyroKinetics4dTests
    using Test

    # The coordinate transformation is not applied by the equations of motion, but it is exported,
    # so its round trip is still checked.
    q, params = initial_conditions_trapped()
    @test transform_q̃_to_q(0.0, transform_q_to_q̃(0.0, q, params), params) ≈ q atol = 1E-6

    test_gyro_kinetics_4d_erk4(guiding_center_4d_ode(initial_conditions_trapped()...))
    test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_trapped()...))
    test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_passing()...))
    test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_barely_trapped()...))
    test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_passing()...))
    test_gyro_kinetics_4d_glrk(guiding_center_4d_ode(initial_conditions_deeply_trapped()...))

    test_gyro_kinetics_4d_strang(guiding_center_4d_sode(initial_conditions_trapped()...))
    test_gyro_kinetics_4d_strang(guiding_center_4d_sode(initial_conditions_barely_passing()...))
    test_gyro_kinetics_4d_strang(guiding_center_4d_sode(initial_conditions_deeply_passing()...))

end


@safetestset "Gyrokinetic GC Model: the vector field is the guiding centre one, rescaled by B*∥                    " begin

    using ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dSolovevIterXpoint
    using ChargedParticleDynamics.GuidingCenter4d.SolovevIterXpoint
    using Test

    # The notes absorb the phase-space Jacobian into the distribution function, f̃ = B*∥ f, which
    # multiplies the guiding centre characteristics through by B*∥. The two vector fields must
    # therefore be exactly proportional, with B*∥ = `ωabs` the factor — same orbits, traversed in
    # the reparametrised time dt = B*∥ ds.
    GK = GuidingCenter4dSolovevIterXpoint
    GC = SolovevIterXpoint

    params = (μ=1E-2,)

    for q in ([6.2, 0.3, 0.0, 3.4E-1], [5.5, -0.8, 1.1, -2.0E-1], [7.0, 0.5, 2.0, 5.0E-1])
        vgk = zeros(4)
        GK.v(vgk, 0.0, q, params)

        vgc = zeros(4)
        GC.guiding_center_4d_v(vgc, 0.0, q, params)

        @test vgk ≈ GK.ωabs(0.0, q) .* vgc
    end

    # The six subsystems of the splitting must sum to the full vector field.
    q = [6.2, 0.3, 0.0, 3.4E-1]
    vfull = zeros(4)
    GK.v(vfull, 0.0, q, params)

    vsum = zeros(4)
    for V in (GK.v₁, GK.v₂, GK.v₃, GK.v₄, GK.v₅, GK.v₆)
        vᵢ = zeros(4)
        V(vᵢ, 0.0, q, params)
        vsum .+= vᵢ
    end

    @test vsum ≈ vfull

end


@safetestset "Gyrokinetic GC Model: the splitting preserves phase-space volume                                     " begin

    using ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dSolovevIterXpoint
    using ..GyroKinetics4dTests
    using GeometricIntegrators
    using LinearAlgebra
    using Test

    # Each subsystem freezes two of the four variables and is symplectic in the other two, so a
    # symplectic method applied to each preserves phase-space volume *exactly*, at any step size.
    # A generic method preserves it only to its order. Build the Jacobian of the one-step map by
    # central differences and compare the determinant of the two.
    q₀, params = initial_conditions_deeply_passing()

    function jacdet(step, Δs; h=1E-5)
        J = zeros(4, 4)
        for j in 1:4
            qp = copy(q₀); qp[j] += h
            qm = copy(q₀); qm[j] -= h
            J[:, j] = (step(qp, Δs) .- step(qm, Δs)) ./ 2h
        end
        det(J)
    end

    split(q, Δs) = integrate(guiding_center_4d_sode(q, params; tstep=Δs, tspan=(0.0, Δs)),
                             strang_composition()).q[end]
    euler(q, Δs) = integrate(guiding_center_4d_ode(q, params; tstep=Δs, tspan=(0.0, Δs)),
                             ExplicitEuler()).q[end]

    # The central differences bottom out around 1E-11; the splitting sits there for every step size
    # while the volume error of the first-order method grows linearly with it.
    for Δs in (1E-3, 1E-2, 3E-2)
        @test abs(jacdet(split, Δs) - 1) < 1E-9
    end

    @test abs(jacdet(euler, 1E-2) - 1) > 1E-7
    @test abs(jacdet(euler, 3E-2) / jacdet(euler, 1E-2) - 1) > 0   # grows with the step

end
