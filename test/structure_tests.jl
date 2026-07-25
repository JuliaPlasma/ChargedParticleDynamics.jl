#
# Structural tests that do not integrate anything.
#
# These are cheap and cover the whole cross product of model families and equilibria, which the
# integration tests cannot: they run one or two equilibria per family because each integration is
# expensive. Every check here corresponds to a class of bug that was found by hand in the audit and
# that no test would have caught:
#
#   * `initial_conditions_*` that raise (a `merge` called with a keyword rather than a named tuple),
#   * problem constructors that raise (a callback referring to a name defined in a sibling file),
#   * analytic derivatives that disagree with the function they claim to differentiate,
#   * projection forces that are not linear in the multiplier, or ignore it altogether,
#   * `default_parameters` methods that are missing or do not honour their type argument.
#

using SafeTestsets


module StructureTestUtils

using ChargedParticleDynamics

export eachmodel, central_difference

const CPD = ChargedParticleDynamics

const FAMILIES = (:ChargedParticle3d, :GuidingCenter3d, :GuidingCenter4d,
                  :PauliParticle3d, :GyroKinetics4d)

"""
Apply `f(family, name, module)` to every equilibrium submodule of every model family.
"""
function eachmodel(f, families = FAMILIES)
    for fam in families
        F = getfield(CPD, fam)
        for mname in sort(collect(names(F; all = true)), by = string)
            startswith(string(mname), "#") && continue
            M = try
                getfield(F, mname)
            catch
                continue
            end
            (M isa Module && M !== F) || continue
            f(fam, mname, M)
        end
    end
end

"Central difference of `f` with respect to component `i` of `x`."
function central_difference(f, x, i; h = 1e-6)
    xp = copy(x); xp[i] += h
    xm = copy(x); xm[i] -= h
    (f(xp) - f(xm)) / 2h
end

end


@safetestset "Every initial_conditions_* constructs                                                               " begin
    using ..StructureTestUtils
    using Test

    n = 0
    eachmodel() do fam, mname, M
        for f in names(M; all = true)
            startswith(string(f), "initial_conditions_") || continue
            fn = try getfield(M, f) catch; continue end
            fn isa Function || continue
            @test (fn(); true)
            n += 1
        end
    end
    @test n > 100   # guards against the loop silently finding nothing
end


@safetestset "Every module has a type-respecting default_parameters                                               " begin
    using ..StructureTestUtils
    using Test

    n = 0
    eachmodel() do fam, mname, M
        @test isdefined(M, :default_parameters)
        isdefined(M, :default_parameters) || return
        p = M.default_parameters()
        @test p isa NamedTuple
        # The type argument must propagate, as it does in GeometricProblems.
        p32 = M.default_parameters(Float32)
        @test all(v -> v isa Float32, values(p32))
        n += 1
    end
    @test n > 40
end


@safetestset "Every problem constructor builds                                                                    " begin
    using ..StructureTestUtils
    using Test

    ctors = Dict(
        :GuidingCenter4d   => ["guiding_center_4d_ode", "guiding_center_4d_iode",
                               "guiding_center_4d_lode", "guiding_center_4d_iode_λ",
                               "guiding_center_4d_dg", "guiding_center_4d_formal_lagrangian"],
        :GuidingCenter3d   => ["hode", "hode_canonical"],
        :ChargedParticle3d => ["charged_particle_3d_ode", "charged_particle_3d_iode",
                               "charged_particle_3d_lode", "charged_particle_3d_pode"],
        :PauliParticle3d   => ["pauli_particle_3d_pode", "pauli_particle_3d_hode",
                               "pauli_particle_3d_iode"],
    )

    icsnames = (:initial_conditions_barely_passing, :initial_conditions_deeply_passing,
                :initial_conditions_default, :initial_conditions_dipole,
                :initial_conditions_quadratic, :initial_conditions_trapped,
                :initial_conditions_pauli)

    n = 0
    eachmodel(collect(keys(ctors))) do fam, mname, M
        ics = nothing
        for name in icsnames
            if isdefined(M, name)
                ics = getfield(M, name)()
                break
            end
        end
        ics === nothing && return
        args = ics isa NamedTuple ? values(ics) : ics
        for c in ctors[fam]
            isdefined(M, Symbol(c)) || continue
            @test (getfield(M, Symbol(c))(args...); true)
            n += 1
        end
    end
    @test n > 80
end


@safetestset "3D guiding centre: second derivatives of H match finite differences                                 " begin
    using ChargedParticleDynamics.GuidingCenter3d
    using ..StructureTestUtils
    using Test

    # The metric cross-term of ∂²H/∂qᵢ∂qⱼ carried a factor of two that is only correct on the
    # diagonal, which made the Hessian asymmetric — by 65 % on the ITER Solov'ev X-point — in every
    # curvilinear coordinate system. Both the symmetry and the value are checked here.
    for M in (GuidingCenter3d.SolovevIterXpoint,        # cylindrical
              GuidingCenter3d.TokamakMediumCartesian,   # cartesian
              GuidingCenter3d.TokamakSmallCylindrical,  # cylindrical
              GuidingCenter3d.TokamakSmallToroidal)     # toroidal
        ic = M.initial_conditions_barely_passing()
        t, q, p, par = 0.0, ic.q, ic.p, ic.params

        dH  = (M.dHdq₁, M.dHdq₂, M.dHdq₃)
        d2H = [M.d²Hdq₁dq₂ M.d²Hdq₁dq₃ M.d²Hdq₂dq₃]

        for (i, j, f) in ((1, 2, M.d²Hdq₁dq₂), (1, 3, M.d²Hdq₁dq₃), (2, 3, M.d²Hdq₂dq₃),
                          (1, 1, M.d²Hdq₁dq₁), (2, 2, M.d²Hdq₂dq₂), (3, 3, M.d²Hdq₃dq₃))
            ana = f(t, q, p, par)
            # ∂/∂qⱼ of ∂H/∂qᵢ and ∂/∂qᵢ of ∂H/∂qⱼ must both reproduce it
            @test ana ≈ central_difference(x -> dH[i](t, x, p, par), q, j) rtol = 1e-5
            @test ana ≈ central_difference(x -> dH[j](t, x, p, par), q, i) rtol = 1e-5
        end
    end
end


@safetestset "4D guiding centre: ḡ matches finite differences of the one-form gradient                            " begin
    using ChargedParticleDynamics.GuidingCenter4d
    using ..StructureTestUtils
    using Test

    # ḡ_k(t,q,v) = q_l ∂²ϑ_l/∂x_k∂x_j v_j, the term the κ-modified one-form contributes to the
    # force. Two of the four ∂²ϑ/∂x_k∂x_j blocks wrote to the wrong matrix row and used the wrong
    # vector-potential component, and the helpers called them with the arguments transposed.
    cases = ((GuidingCenter4d.SymmetricField,            [0.35, -0.22, 0.4, 0.3]),
             (GuidingCenter4d.TokamakMediumCylindrical,  [2.4, 0.15, 0.3, 0.2]),
             (GuidingCenter4d.TokamakSmallToroidal,      [0.06, 0.4, 0.3, 3e-4]))

    for (M, q0) in cases
        t = 0.0
        q = collect(float.(q0))
        v = [0.7, -0.3, 0.45, 0.2]
        buf = zeros(4, 4)

        dϑk_dot_q(x, k) = (M.dϑ(buf, t, x); sum(buf[l, k] * q[l] for l in 1:4))

        for (k, ḡ) in enumerate((M.g̅₁, M.g̅₂, M.g̅₃, M.g̅₄))
            ana = ḡ(t, q, v)
            num = sum(central_difference(x -> dϑk_dot_q(x, k), q, j) * v[j] for j in 1:4)
            # Absolute tolerance as well: ḡ₄ vanishes identically where b is constant, and the
            # central difference of an O(1) quantity at h = 1e-6 carries ~1e-10 of noise.
            @test isapprox(ana, num; rtol = 1e-5, atol = 1e-9)
        end
    end
end


@safetestset "Projection forces are linear in the multiplier                                                      " begin
    using ChargedParticleDynamics.PauliParticle3d
    using ChargedParticleDynamics.ChargedParticle3d
    using Test

    # g = λ ⋅ ∇ϑ must be linear in λ. The Pauli version was built from `v` and ignored `λ`
    # entirely; the canonical charged particle version was quadratic in λ.
    cases = ((PauliParticle3d.TokamakSmallCylindrical,
              PauliParticle3d.TokamakSmallCylindrical.pauli_particle_3d_iode_g,
              [1.05, 0.0, 0.0], (μ = 2.31e-6,)),
             (ChargedParticle3d.SolovevIterXpoint,
              ChargedParticle3d.SolovevIterXpoint.charged_particle_3d_iode_g,
              [2.5, 0.0, 0.0], NamedTuple()))

    for (M, g, q, par) in cases
        t = 0.0
        v = [1.3e-3, -4.0e-4, 7.0e-4]
        λ = [0.31, -0.17, 0.42]

        g1 = zeros(3); g(g1, t, q, v, λ, par)
        g2 = zeros(3); g(g2, t, q, v, 2 .* λ, par)
        g0 = zeros(3); g(g0, t, q, v, zero(λ), par)

        @test g2 ≈ 2 .* g1
        @test all(iszero, g0)
        # and it must actually depend on λ
        @test !all(iszero, g1)
    end
end


@safetestset "Toroidal momentum is conserved                                                                      " begin
    using ChargedParticleDynamics.GuidingCenter4d
    using GeometricIntegrators
    using Test

    # `toroidal_momentum` used to be R ϑ₃ everywhere, which is not conserved: the factor of R is
    # spurious, since ϑ₃ is already the covariant φ-component. On the small tokamak the relative
    # variation over 10³ time units was 3e-3 for R ϑ₃ against 2e-13 for ϑ₃. In cartesian
    # coordinates the third coordinate is z rather than an angle, so neither is right there and
    # the momentum is the generator of rotation about the z-axis instead.
    #
    # The tolerances below are far tighter than the old definition could ever reach, which is the
    # point of the test; they are loose relative to what the integrator achieves.
    cases = ((GuidingCenter4d.TokamakSmallCylindrical, (tstep = 10.0, tspan = (0.0, 1e3)), 1e-10),
             (GuidingCenter4d.TokamakSmallToroidal,    (tstep = 10.0, tspan = (0.0, 1e3)), 1e-10),
             (GuidingCenter4d.TokamakSmallCartesian,   (tstep = 10.0, tspan = (0.0, 1e3)), 1e-8),
             (GuidingCenter4d.SolovevIterXpoint,       (tstep = 1.0,  tspan = (0.0, 1e2)), 1e-10))

    for (M, kw, tol) in cases
        prob = M.guiding_center_4d_ode(M.initial_conditions_barely_passing()...; kw...)
        sol  = integrate(prob, Gauss(2))
        _, perr = M.compute_toroidal_momentum_error(sol)
        @test maximum(abs(perr[i]) for i in eachindex(perr)) < tol
    end
end


@safetestset "3D guiding centre diagnostics                                                                       " begin
    using ChargedParticleDynamics.GuidingCenter3d
    using GeometricIntegrators
    using Test

    # `guiding_center_3d_diagnostics.jl` was an empty file, so the 3D model had no diagnostics at
    # all, and the `toroidal_momentum` its equilibria defined indexed `q[4]` of a three-component
    # state. The constraints are the model-specific diagnostic: they vanish along the exact flow,
    # so their magnitude is the drift off the manifold on which p = ϑ.
    mx(ds) = maximum(abs(ds[i]) for i in eachindex(ds))

    for M in (GuidingCenter3d.SolovevIterXpoint, GuidingCenter3d.TokamakMediumCartesian)
        prob = M.hode(M.initial_conditions_barely_passing()...; tstep = 0.1, tspan = (0.0, 1e2))
        sol  = integrate(prob, PartitionedGauss(2); initialguess = MidpointExtrapolation(5))

        _, eerr = M.compute_energy_error(sol)
        @test mx(eerr) < 1e-8

        c = M.compute_constraints(sol)
        @test mx(c.g₁) < 1e-8
        @test mx(c.g₂) < 1e-8

        if isdefined(M, :toroidal_momentum)
            _, perr = M.compute_toroidal_momentum_error(sol)
            @test mx(perr) < 1e-10
        end
    end
end


@safetestset "Symplectic two-forms are antisymmetric                                                              " begin
    using ChargedParticleDynamics.GuidingCenter4d
    using ChargedParticleDynamics.ChargedParticle3d
    using Test

    let M = GuidingCenter4d.TokamakSmallCylindrical
        Ω = zeros(4, 4)
        M.ω(Ω, 0.0, [1.05, 0.0, 0.0, 4.3e-4])
        @test Ω ≈ -transpose(Ω)
    end

    let M = ChargedParticle3d.TokamakSmallNoncanonical
        Ω = zeros(6, 6)
        M.ω(Ω, 0.0, [1.05, 0.0, 0.0, 1.0e-3, 0.0, -4.0e-4], NamedTuple())
        @test Ω ≈ -transpose(Ω)
    end

    let M = ChargedParticle3d.SolovevIterXpoint
        Ω = zeros(6, 6)
        M.ω(Ω, 0.0, [2.5, 0.0, 0.0], NamedTuple())
        @test Ω ≈ -transpose(Ω)
    end
end
