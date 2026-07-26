
using SafeTestsets

module GyroKinetics4dTests

using ChargedParticleDynamics.GyroKinetics4d
using GeometricIntegrators
using Test

# See `guiding_center_3d_tests.jl` for why `f_abstol` has to stay above the residual's round-off
# floor and why `f_reltol` is left at its default.
const options = (f_abstol=1E-12, max_iterations=50, warn_iterations=50)

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

    test_gyro_kinetics_4d_erk4(odeproblem(initial_conditions_trapped()))
    test_gyro_kinetics_4d_glrk(odeproblem(initial_conditions_trapped()))
    test_gyro_kinetics_4d_glrk(odeproblem(initial_conditions_barely_passing()))
    test_gyro_kinetics_4d_glrk(odeproblem(initial_conditions_barely_trapped()))
    test_gyro_kinetics_4d_glrk(odeproblem(initial_conditions_deeply_passing()))
    test_gyro_kinetics_4d_glrk(odeproblem(initial_conditions_deeply_trapped()))

    test_gyro_kinetics_4d_strang(sodeproblem(initial_conditions_trapped()))
    test_gyro_kinetics_4d_strang(sodeproblem(initial_conditions_barely_passing()))
    test_gyro_kinetics_4d_strang(sodeproblem(initial_conditions_deeply_passing()))

end


@safetestset "Gyrokinetic GC Model: the vector field is the guiding centre one, rescaled by B*∥                    " begin

    using ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dSolovevIterXpoint
    using ChargedParticleDynamics.GuidingCenter4d.SolovevIterXpoint
    using Test

    # The notes absorb the phasespace Jacobian into the distribution function, f̃ = B*∥ f, which
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


@safetestset "Gyrokinetic GC Model: the splitting preserves phasespace volume                                      " begin

    using ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dSolovevIterXpoint
    using ..GyroKinetics4dTests
    using GeometricIntegrators
    using LinearAlgebra
    using Test

    # Each subsystem freezes two of the four variables and is symplectic in the other two, so a
    # symplectic method applied to each preserves phasespace volume *exactly*, at any step size.
    # A generic method preserves it only to its order. Build the Jacobian of the one-step map by
    # central differences and compare the determinant of the two.
    ics = initial_conditions_deeply_passing()
    q₀, params = ics.q, ics.params

    function jacdet(step, Δs; h=1E-5)
        J = zeros(4, 4)
        for j in 1:4
            qp = copy(q₀); qp[j] += h
            qm = copy(q₀); qm[j] -= h
            J[:, j] = (step(qp, Δs) .- step(qm, Δs)) ./ 2h
        end
        det(J)
    end

    split(q, Δs) = integrate(sodeproblem(q; parameters=params, timestep=Δs, timespan=(0.0, Δs)),
                             strang_composition()).q[end]
    euler(q, Δs) = integrate(odeproblem(q; parameters=params, timestep=Δs, timespan=(0.0, Δs)),
                             ExplicitEuler()).q[end]

    # The central differences bottom out around 1E-11; the splitting sits there for every step size
    # while the volume error of the first-order method grows linearly with it.
    for Δs in (1E-3, 1E-2, 3E-2)
        @test abs(jacdet(split, Δs) - 1) < 1E-9
    end

    @test abs(jacdet(euler, 1E-2) - 1) > 1E-7
    @test abs(jacdet(euler, 3E-2) / jacdet(euler, 1E-2) - 1) > 0   # grows with the step

end


@safetestset "Gyrokinetic GC Model: every equilibrium integrates and preserves volume                " begin

    using ChargedParticleDynamics
    using ChargedParticleDynamics.GyroKinetics4d
    using ..GyroKinetics4dTests
    using GeometricIntegrators
    using LinearAlgebra
    using Test

    # Sweeps the whole set rather than the one equilibrium the blocks above cover. The point is the
    # curvilinear ones: `gc_common.jl` carries no metric, because the formulation is intrinsic —
    # ϑᵢ = Aᵢ + u bᵢ in covariant components and Ω = dϑ need none — and the Liouville measure is
    # √det Ω = B*∥ in any chart. So the volume-preserving splitting must work unchanged in
    # cylindrical and toroidal coordinates, and this is the assertion that says so.
    #
    # Each module's time step is the 4D guiding centre's divided by B*∥ at its own initial
    # condition, since the independent variable is the rescaled time with dt = B*∥ ds. Getting that
    # factor wrong by ~2·10² is exactly the bug the audit found in the one shipped equilibrium.
    equilibria = sort(filter(n -> isa(getfield(GyroKinetics4d, n), Module) && n !== :GyroKinetics4d,
                             names(GyroKinetics4d, all = true)))

    @test length(equilibria) == 8

    for name in equilibria
        M = getfield(GyroKinetics4d, name)
        ics = M.initial_conditions_deeply_passing()
        q₀, params = ics.q, ics.params
        Δs = M.DEFAULT_TIMESTEP

        @test integrate(M.odeproblem(ics), Gauss(2)) isa GeometricSolution
        @test integrate(M.sodeproblem(ics), strang_composition()) isa GeometricSolution

        # The rescaled step really is the physical one divided by B*∥: multiplying it back by
        # B*∥ must recover the step of the corresponding 4D guiding centre module, which is the
        # same equilibrium in physical time. This is what catches a step copied from a sibling —
        # the audit found the shipped one wrong by a factor of B*∥ ≈ 2·10².
        short = Symbol(replace(string(name), "GuidingCenter4d" => ""))
        G = getfield(ChargedParticleDynamics.GuidingCenter4d, short)
        Bpar = M.ωabs(0.0, q₀, params)
        @test 0.1 < Δs * Bpar / G.DEFAULT_TIMESTEP < 10.0

        step(q) = integrate(M.sodeproblem(q; parameters = params, timestep = Δs, timespan = (0.0, Δs)),
                            strang_composition()).q[end]
        h = 1E-7 * max(1.0, maximum(abs, q₀))
        J = zeros(4, 4)
        for j in 1:4
            qp = copy(q₀); qp[j] += h
            qm = copy(q₀); qm[j] -= h
            J[:, j] = (step(qp) .- step(qm)) ./ 2h
        end
        @test abs(det(J) - 1) < 1E-7
    end

end
