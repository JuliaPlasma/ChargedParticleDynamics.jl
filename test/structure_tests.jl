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
        :GuidingCenter4d   => ["odeproblem", "iodeproblem",
                               "lodeproblem", "iodeproblem_λ",
                               "iodeproblem_dg", "lodeproblem_formal_lagrangian"],
        :GuidingCenter3d   => ["hodeproblem", "hodeproblem_canonical"],
        :ChargedParticle3d => ["odeproblem", "iodeproblem",
                               "lodeproblem", "podeproblem"],
        :PauliParticle3d   => ["podeproblem", "hodeproblem",
                               "iodeproblem"],
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

        # Every family now returns the same shape — `(q = …, params = …)`, plus `p` or `v` where
        # the model has one — and every constructor takes that named tuple directly. This used to
        # branch on `ics isa NamedTuple` and splat, because the 3D guiding centre returned a named
        # tuple while the other three returned a bare tuple.
        @test ics isa NamedTuple
        @test haskey(ics, :q) && haskey(ics, :params)

        for c in ctors[fam]
            isdefined(M, Symbol(c)) || continue
            @test (getfield(M, Symbol(c))(ics); true)
            n += 1
        end
    end
    @test n > 80
end


@safetestset "Noncanonical charged particle: constructors build and the splitting is consistent                   " begin
    using ChargedParticleDynamics.ChargedParticle3d
    using GeometricIntegrators
    using Test

    # The noncanonical modules define no `initial_conditions_*`, so the constructor sweep above
    # skips them entirely and none of their problems was covered. Every constructor here defaults
    # its argument to the module's `qᵢ`.
    cartesian = (ChargedParticle3d.SingularField,
                 ChargedParticle3d.SymmetricField,
                 ChargedParticle3d.ThetaPinchNoncanonical)
    curvilinear = (ChargedParticle3d.TokamakSmallNoncanonical,)

    for M in (cartesian..., curvilinear...)
        @test (M.odeproblem(); true)
        @test (M.iodeproblem(); true)
        @test (M.lodeproblem(); true)
    end

    # The Boris splitting is the splitting of the *cartesian* Lorentz force. In a curvilinear chart
    # the frozen-position kick is quadratic in v, so the Cayley transform does not solve it and the
    # splitting has no exact flow; the constructor refuses rather than integrating a different
    # model. This is the assertion that fails if that guard is ever dropped.
    for M in curvilinear
        @test !M.has_trivial_metric(0.0, collect(float.(M.qᵢ)))
        @test_throws ArgumentError M.sodeproblem()
    end

    for M in cartesian
        @test M.has_trivial_metric(0.0, collect(float.(M.qᵢ)))
        @test (M.sodeproblem(); true)

        # `sodeproblem` used to accept `parameters` and then not forward it to `SODEProblem`, so
        # the keyword — and anything the named-tuple form put into it — was silently dropped.
        @test parameters(M.sodeproblem(; parameters = (a = 1,))) == (a = 1,)

        # The two maps of the splitting must sum to the full vector field. To first order in h,
        #     (fv(q₀) - q₀)/h + (fx(q₀) - q₀)/h  →  charged_particle_3d_v(q₀) ,
        # which is the Lie-Trotter consistency condition, and is what pins the splitting to the
        # model rather than to some subset of its terms. The electric field used to be absent from
        # both maps; every shipped equilibrium has φ = 0, so this passes either way today, but it
        # is the assertion that fails the moment a noncanonical equilibrium with a potential is
        # added — which is how that omission should have been caught.
        q₀ = collect(float.(M.qᵢ))
        h = 1e-7
        qv = zero(q₀); M.charged_particle_3d_sode_fv(qv, h, q₀, 0.0, NamedTuple())
        qx = zero(q₀); M.charged_particle_3d_sode_fx(qx, h, q₀, 0.0, NamedTuple())
        vfull = zero(q₀); M.charged_particle_3d_v(vfull, 0.0, q₀, NamedTuple())

        @test isapprox((qv .- q₀) ./ h .+ (qx .- q₀) ./ h, vfull; rtol = 1e-5, atol = 1e-9)

        # and each solution map must agree with the vector field of its own substep
        vv = zero(q₀); M.charged_particle_3d_sode_vv(vv, 0.0, q₀, NamedTuple())
        vx = zero(q₀); M.charged_particle_3d_sode_vx(vx, 0.0, q₀, NamedTuple())

        @test isapprox((qv .- q₀) ./ h, vv; rtol = 1e-5, atol = 1e-9)
        @test isapprox((qx .- q₀) ./ h, vx; rtol = 1e-5, atol = 1e-9)

        # The SODE integrates through a composition of the two exact-solution maps.
        sol = integrate(M.sodeproblem(; timespan = (0.0, 1.0), timestep = 0.1),
                        Composition(Strang()))
        @test sol isa GeometricSolution
    end
end


@safetestset "Noncanonical charged particle: the vector field is the one its own ω and H generate            " begin
    using ChargedParticleDynamics.ChargedParticle3d
    using LinearAlgebra
    using Test

    # `charged_particle_3d_v` used to be the plain cartesian Lorentz force while `ω`, `dϑ`, `ϑ` and
    # `hamiltonian` beside it all carried the metric, so the toroidal equilibrium integrated a
    # system inconsistent with its own structure. It is now derived from Ω ż = -∇H, and these are
    # the two assertions that hold it there:
    #
    #   * `dH` really is the gradient of `hamiltonian` — it used to be the gradient of a different,
    #     metric-free Hamiltonian, and referred to a `dφdxᵢ` that the field-code generator does not
    #     even emit, so it raised an `UndefVarError` on every call;
    #   * the vector field really is -Ω⁻¹∇H for the module's own `ω`.
    for M in (ChargedParticle3d.SingularField,
              ChargedParticle3d.SymmetricField,
              ChargedParticle3d.ThetaPinchNoncanonical,
              ChargedParticle3d.TokamakSmallNoncanonical)

        # displaced off the initial condition, which sits on a symmetry axis where several of the
        # metric derivatives vanish and would not exercise the new terms
        q = collect(float.(M.qᵢ)) .+ 0.01 .* [1.0, 2.0, 3.0, 1.0, 2.0, 3.0]

        g = zeros(6); M.dH(g, 0.0, q)
        fd = [(M.hamiltonian(0.0, (h = zeros(6); h[i] = 1e-6; q .+ h)) -
               M.hamiltonian(0.0, (h = zeros(6); h[i] = 1e-6; q .- h))) / 2e-6 for i in 1:6]
        @test isapprox(g, fd; rtol = 1e-6, atol = 1e-10)

        Ω = zeros(6, 6); M.ω(Ω, 0.0, q, NamedTuple())
        v = zeros(6);    M.charged_particle_3d_v(v, 0.0, q, NamedTuple())
        @test isapprox(v, -Ω \ g; rtol = 1e-10, atol = 1e-14)
    end
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


@safetestset "Projection forces are the gradient of the one-form                                                  " begin
    using ChargedParticleDynamics.PauliParticle3d
    using ChargedParticleDynamics.ChargedParticle3d
    using ..StructureTestUtils
    using Test

    # g must be gᵢ = (∂ϑⱼ/∂zⁱ) λʲ for the model's own one-form ϑ on its own state z. Three separate
    # ways of getting that wrong have been found in this package, so all three are checked:
    #
    #   * the Pauli version was built from `v` and ignored `λ` entirely,
    #   * the canonical charged particle version was quadratic in λ,
    #   * the noncanonical charged particle version was linear in λ and *looked* fine, but was
    #     written for a one-form without the metric — no ∂gⱼⱼ/∂xⁱ vʲ terms in the position block and
    #     an identity velocity block where ∂ϑⱼ/∂vⁱ = gᵢᵢ δᵢⱼ. Linearity alone does not catch that,
    #     which is why the finite-difference comparison below is the real assertion.
    #
    # `z` is the state g is indexed by and `oneform(z, v)` evaluates ϑ there, both of length `nz`.
    # `λ` has the dimension of the state too — for the noncanonical model that is six, and its last
    # three components multiply the identically vanishing ϑ₄₋₆.
    let
        Mp = PauliParticle3d.TokamakSmallCylindrical
        Mc = ChargedParticle3d.SolovevIterXpoint
        Mn = ChargedParticle3d.TokamakSmallNoncanonical

        v = [1.3e-3, -4.0e-4, 7.0e-4]
        λ3 = [0.31, -0.17, 0.42]
        λ6 = [0.31, -0.17, 0.42, 0.23, -0.51, 0.11]

        # (module, g, state, oneform, nz, λ, params)
        cases = (
            (Mp, Mp.pauli_particle_3d_iode_g, [1.05, 0.0, 0.0],
             (z, vv) -> (θ = zeros(3); Mp.ϑ(θ, 0.0, z, vv); θ), 3, λ3, (μ = 2.31e-6,)),
            (Mc, Mc.charged_particle_3d_iode_g, [2.5, 0.0, 0.0],
             (z, vv) -> (θ = zeros(3); Mc.ϑ(θ, 0.0, z, vv); θ), 3, λ3, NamedTuple()),
            (Mn, Mn.charged_particle_3d_iode_g, collect(float.(Mn.qᵢ)),
             (z, vv) -> (θ = zeros(6); Mn.ϑ(θ, 0.0, z); θ), 6, λ6, NamedTuple()),
        )

        for (M, g, z, oneform, nz, λ, par) in cases
            t = 0.0

            ga = zeros(nz); g(ga, t, z, v, λ, par)

            # gᵢ = ∂/∂zⁱ (ϑ ⋅ λ), by central differences
            for i in 1:nz
                num = central_difference(x -> sum(oneform(x, v) .* λ), z, i)
                @test isapprox(ga[i], num; rtol = 1e-5, atol = 1e-9)
            end

            # and it must be linear in λ, which the finite difference cannot see
            g2 = zeros(nz); g(g2, t, z, v, 2 .* λ, par)
            g0 = zeros(nz); g(g0, t, z, v, zero(λ), par)

            @test g2 ≈ 2 .* ga
            @test all(iszero, g0)
            @test !all(iszero, ga)
        end
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
    cases = ((GuidingCenter4d.TokamakSmallCylindrical, (timestep = 10.0, timespan = (0.0, 1e3)), 1e-10),
             (GuidingCenter4d.TokamakSmallToroidal,    (timestep = 10.0, timespan = (0.0, 1e3)), 1e-10),
             (GuidingCenter4d.TokamakSmallCartesian,   (timestep = 10.0, timespan = (0.0, 1e3)), 1e-8),
             (GuidingCenter4d.SolovevIterXpoint,       (timestep = 1.0,  timespan = (0.0, 1e2)), 1e-10))

    for (M, kw, tol) in cases
        prob = M.odeproblem(M.initial_conditions_barely_passing(); kw...)
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

    # This block keeps its 1000 steps deliberately: it is where the long-time behaviour of the 3D
    # model is checked, now that the integration tests in `guiding_center_3d_tests.jl` run 100. It
    # also keeps the library's default tolerance, so that it measures the conservation an
    # unconfigured `integrate` delivers — but it caps the iteration count, because a handful of steps
    # otherwise run to the default limit of 1000 and account for most of the block's runtime (48 s of
    # 49 s on the Solov'ev case). Capping to 50 leaves every error below bit-for-bit unchanged.
    #
    # `f_abstol` has to be restated even though it is the library default:
    # `GeometricIntegratorsBase` replaces its whole `default_options` set as soon as the caller passes
    # any option, and `SimpleSolvers`' own default is `f_abstol = 0`, which no residual can meet.
    options = (f_abstol = 8eps(), max_iterations = 50, warn_iterations = 50)

    for M in (GuidingCenter3d.SolovevIterXpoint, GuidingCenter3d.TokamakMediumCartesian)
        prob = M.hodeproblem(M.initial_conditions_barely_passing(); timestep = 0.1, timespan = (0.0, 1e2))
        sol  = integrate(prob, PartitionedGauss(2); initialguess = MidpointExtrapolation(5), options...)

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


@safetestset "Symplectic two-forms are antisymmetric and correctly sized                                          " begin
    using ChargedParticleDynamics.GuidingCenter4d
    using ChargedParticleDynamics.ChargedParticle3d
    using ..StructureTestUtils
    using Test

    # The buffer is sized from the state rather than hard-coded, because that is precisely what an
    # integrator does: the `LODE` two-form is called with a `D × D` array, `D = length(q)`. The
    # canonical charged particle used to return the 6 × 6 canonical form on (q,p) for a model whose
    # `q` has three components, which a hard-coded `zeros(6, 6)` here hid.
    let M = GuidingCenter4d.TokamakSmallCylindrical
        q = [1.05, 0.0, 0.0, 4.3e-4]
        Ω = zeros(length(q), length(q))
        M.ω(Ω, 0.0, q)
        @test Ω ≈ -transpose(Ω)
    end

    let M = ChargedParticle3d.TokamakSmallNoncanonical
        q = collect(float.(M.qᵢ))
        Ω = zeros(length(q), length(q))
        M.ω(Ω, 0.0, q, NamedTuple())
        @test Ω ≈ -transpose(Ω)
    end

    let M = ChargedParticle3d.SolovevIterXpoint
        q = [2.5, 0.0, 0.0]
        v = [1.3e-3, -4.0e-4, 7.0e-4]
        Ω = zeros(length(q), length(q))
        M.ω(Ω, 0.0, q, v, NamedTuple())
        @test Ω ≈ -transpose(Ω)

        # Ωᵢⱼ = ∂ϑᵢ/∂qʲ - ∂ϑⱼ/∂qⁱ, against central differences of the one-form itself.
        ϑi(x, i) = (θ = zeros(3); M.ϑ(θ, 0.0, x, v); θ[i])
        for i in 1:3, j in 1:3
            num = central_difference(x -> ϑi(x, i), q, j) - central_difference(x -> ϑi(x, j), q, i)
            @test isapprox(Ω[i, j], num; rtol = 1e-5, atol = 1e-9)
        end
    end
end
