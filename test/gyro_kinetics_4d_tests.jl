
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


@safetestset "Gyrokinetic GC Model: the vector field is the guiding centre one, rescaled by ωabs                   " begin

    using ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dSolovevIterXpoint
    using ChargedParticleDynamics.GuidingCenter4d.SolovevIterXpoint
    using GeometricIntegrators
    using Test

    # The notes absorb the phasespace Jacobian into the distribution function, f̃ = ωabs f, which
    # multiplies the guiding centre characteristics through by it. The two vector fields must
    # therefore be exactly proportional, with `ωabs` the factor — same orbits, traversed in the same
    # direction, in the reparametrised time dt = `ωabs` ds.
    #
    # `ωabs = J B*∥` is positive in every chart, including this left-handed Solov'ev one, because it
    # carries the chart's `ORIENTATION`; see `ωabs` in `gc_common.jl`. The equality below is what
    # pins the two together, and the "one orbit in three charts" block at the end of this file is
    # what pins the sign.
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

    # ... and it must be the field the SODE hands the library as `v̄`. That is what the `q̇` of every
    # solution and the backward extrapolation of the solution step's pre-history are computed with,
    # and `SODE` — unlike `ODE`, which defaults `v̄` to `v` — defaults it to a function that writes
    # nothing. Left at that default the pre-history came back equal to the initial condition
    # bit-for-bit, so the Hermite initial guess of every subsystem solve of the composition warned
    # and fell back to a constant. The trajectories were unaffected, which is why only the warning
    # gave it away; this is the assertion that fails if the keyword is ever dropped again.
    v̄ = zeros(4)
    initialguess(GK.sodeproblem(q; parameters=params)).v(v̄, 0.0, q, params)

    @test v̄ == vfull

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
    # √det Ω in any chart, which is `ωabs` as this module computes it. So the volume-preserving
    # splitting must work unchanged in cylindrical and toroidal coordinates, and this is the assertion
    # that says so. Applying the chart's `ORIENTATION` to the field does not disturb it: negating a
    # divergence-free field leaves it divergence-free.
    #
    # Each module's time step is the 4D guiding centre's divided by the rescaling factor at its own
    # initial condition, since the independent variable is the rescaled time with dt = `ωabs` ds.
    # Getting that factor wrong by ~2·10² is exactly the bug the audit found in the one shipped
    # equilibrium.
    #
    # `ωabs` is the phasespace Jacobian `J B*∥`, and it is **positive in every chart**. The bare
    # contraction it is built from is `det(DF) · B*∥` — `∂₂ϑ₃ - ∂₃ϑ₂` is the coordinate curl, which
    # carries the signed Jacobian — so each module declares an `ORIENTATION` that restores it; see
    # `ωabs` in `gc_common.jl`. Both facts are asserted below: that the declared constant really is
    # the handedness of the chart, and that the factor that reaches the dynamics is positive.
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

        # The rescaled step really is the physical one divided by the rescaling factor: multiplying
        # it back must recover the step of the corresponding 4D guiding centre module, which is the
        # same equilibrium in physical time. This is what catches a step copied from a sibling —
        # the audit found the shipped one wrong by a factor of ≈2·10².
        short = Symbol(replace(string(name), "GuidingCenter4d" => ""))
        G = getfield(ChargedParticleDynamics.GuidingCenter4d, short)
        ωfac = M.ωabs(0.0, q₀, params)
        @test 0.1 < abs(Δs * ωfac / G.DEFAULT_TIMESTEP) < 10.0

        # The factor that reaches the dynamics is positive, so `s` runs with `t`. This is the
        # assertion that would have caught the sign the coordinate curl used to carry into the
        # vector field, and it holds for every chart rather than per chart.
        @test ωfac > 0

        # And the module's declared `ORIENTATION` really is its chart's handedness, checked against
        # the bare contraction rather than restated: `ORIENTATION · (bare) = ωabs > 0` means the
        # constant and the chart agree. Declaring the wrong one would silently reverse the orbit.
        bare = M.ω₁(0.0, q₀) * M.dϑ₁dx₄(0.0, q₀) +
               M.ω₂(0.0, q₀) * M.dϑ₂dx₄(0.0, q₀) +
               M.ω₃(0.0, q₀) * M.dϑ₃dx₄(0.0, q₀)
        @test sign(bare) == M.ORIENTATION

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


@safetestset "Gyrokinetic GC Model: one orbit, three charts, one direction                           " begin

    using ChargedParticleDynamics
    using ChargedParticleDynamics.GyroKinetics4d
    using ChargedParticleDynamics.GuidingCenter4d
    using Test

    # The regression test for the sign of the rescaling factor.
    #
    # `ωabs` is built from a *coordinate* curl, which carries `det(DF)`, so before each module
    # declared its `ORIENTATION` the factor was negative in the six left-handed charts — and since
    # the vector field carries the same factor, the model integrated those orbits **backwards**.
    # That is not a reparametrisation of time but chart-dependence: the same physical particle in
    # the same tokamak circulated one way in the cartesian chart and the other way in the
    # cylindrical one. Nothing else in the suite noticed, because everything else is a magnitude —
    # the proportionality `v_gk = ωabs v_gc` holds with either sign, and `|det J| = 1` is
    # insensitive to direction.
    #
    # So this asserts the one thing that is not a magnitude. The small tokamak carries all three
    # charts of a single equilibrium; started from the same physical point with the same `u` and
    # `μ`, the *physical* toroidal velocity `(scale factor) × (coordinate derivative) / ωabs` must
    # be the same number in all three, and must equal the guiding centre's.
    GK = GyroKinetics4d
    G4 = GuidingCenter4d

    params = (μ = 2.448E-6,)
    u = 1.623E-3

    # (gyrokinetic module, guiding centre module, toroidal coordinate index, its scale factor)
    charts = ((GK.GuidingCenter4dTokamakSmallCartesian,   G4.TokamakSmallCartesian,   2, q -> 1.0),
              (GK.GuidingCenter4dTokamakSmallCylindrical, G4.TokamakSmallCylindrical, 3, q -> q[1]),
              (GK.GuidingCenter4dTokamakSmallToroidal,    G4.TokamakSmallToroidal,    3, q -> 1.0 + q[1] * cos(q[2])))

    reference = nothing

    for (Mk, Mg, i, scale) in charts
        q = [Mg.from_cartesian(0, [1.05, 0.0, 0.0])..., u]

        vgk = zeros(4); Mk.v(vgk, 0.0, q, params)
        vgc = zeros(4); Mg.guiding_center_4d_v(vgc, 0.0, q, params)
        ω = Mk.ωabs(0.0, q, params)

        # the gyrokinetic field, un-rescaled, is the guiding centre field — sign included
        @test ω > 0
        @test vgk ≈ ω .* vgc

        # and the physical toroidal velocity is chart-independent, which is the statement that
        # failed while the factor was signed
        toroidal = scale(q) * vgk[i] / ω
        @test toroidal ≈ scale(q) * vgc[i]

        reference === nothing && (reference = toroidal)
        @test toroidal ≈ reference
        @test toroidal > 0
    end

end
