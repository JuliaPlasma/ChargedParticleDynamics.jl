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

An equilibrium module is one that holds a field as `FIELD`. That leaves out the submodules holding
a family's equations, `ChargedParticle3d.Canonical` and `ChargedParticle3d.Noncanonical`.
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
            (M isa Module && M !== F && isdefined(M, :FIELD)) || continue
            f(fam, mname, M)
        end
    end
end

"Central difference of `f` with respect to component `i` of `x`."
function central_difference(f, x, i; h = 1e-6)
    xp = copy(x)
    xp[i] += h
    xm = copy(x)
    xm[i] -= h
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
            fn = try
                getfield(M, f)
            catch
                continue
            end
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
        # Every problem reads its field from `params.field`, and it is the module's own.
        @test p.field === M.FIELD
        # The type argument must propagate to the physical parameters, as it does in
        # GeometricProblems. The field is a value of its own and is not converted.
        p32 = M.default_parameters(Float32)
        @test all(v -> v isa Float32, values(Base.structdiff(p32, (field = nothing,))))
        n += 1
    end
    @test n > 40
end

@safetestset "Every problem constructor builds                                                                    " begin
    using ..StructureTestUtils
    using Test

    ctors = Dict(
        :GuidingCenter4d => ["odeproblem", "iodeproblem",
            "lodeproblem", "iodeproblem_λ",
            "iodeproblem_dg", "lodeproblem_formal_lagrangian"],
        :GuidingCenter3d => ["hodeproblem", "hodeproblem_canonical", "hodeproblem_compact"],
        :ChargedParticle3d => ["odeproblem", "iodeproblem",
            "lodeproblem", "podeproblem"],
        :PauliParticle3d => ["podeproblem", "hodeproblem",
            "iodeproblem"]
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

    # The equations of the formulation, which every one of the modules above forwards to.
    N = ChargedParticle3d.Noncanonical

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
        q₀ = collect(float.(M.qᵢ))
        @test !N.has_trivial_metric(0.0, N.fieldpoint(M.FIELD, 0.0, q₀))
        @test_throws ArgumentError M.sodeproblem()
    end

    for M in cartesian
        params = M.default_parameters()
        q₀ = collect(float.(M.qᵢ))
        @test N.has_trivial_metric(0.0, N.fieldpoint(M.FIELD, 0.0, q₀))
        @test (M.sodeproblem(); true)

        # `sodeproblem` used to accept `parameters` and then not forward it to `SODEProblem`, so
        # the keyword — and anything the named-tuple form put into it — was silently dropped.
        @test parameters(M.sodeproblem(; parameters = (field = M.FIELD, a = 1))).a == 1

        # The two maps of the splitting must sum to the full vector field. To first order in h,
        #     (fv(q₀) - q₀)/h + (fx(q₀) - q₀)/h  →  charged_particle_3d_v(q₀) ,
        # which is the Lie-Trotter consistency condition, and is what pins the splitting to the
        # model rather than to some subset of its terms. The electric field used to be absent from
        # both maps; every shipped equilibrium has φ = 0, so this passes either way today, but it
        # is the assertion that fails the moment a noncanonical equilibrium with a potential is
        # added — which is how that omission should have been caught.
        h = 1e-7
        qv = zero(q₀)
        N.charged_particle_3d_sode_fv(qv, h, q₀, 0.0, params)
        qx = zero(q₀)
        N.charged_particle_3d_sode_fx(qx, h, q₀, 0.0, params)
        vfull = zero(q₀)
        N.charged_particle_3d_v(vfull, 0.0, q₀, params)

        @test isapprox((qv .- q₀) ./ h .+ (qx .- q₀) ./ h, vfull; rtol = 1e-5, atol = 1e-9)

        # and each solution map must agree with the vector field of its own substep
        vv = zero(q₀)
        N.charged_particle_3d_sode_vv(vv, 0.0, q₀, params)
        vx = zero(q₀)
        N.charged_particle_3d_sode_vx(vx, 0.0, q₀, params)

        @test isapprox((qv .- q₀) ./ h, vv; rtol = 1e-5, atol = 1e-9)
        @test isapprox((qx .- q₀) ./ h, vx; rtol = 1e-5, atol = 1e-9)

        # The SODE integrates through a composition of the two exact-solution maps.
        sol = integrate(M.sodeproblem(; timespan = (0.0, 1.0), timestep = 0.1),
            Composition(Strang()))
        @test sol isa GeometricSolution
    end
end

@safetestset "Noncanonical charged particle: the vector field is the one its own ω and H generate                 " begin
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
    N = ChargedParticle3d.Noncanonical

    for M in (ChargedParticle3d.SingularField,
        ChargedParticle3d.SymmetricField,
        ChargedParticle3d.ThetaPinchNoncanonical,
        ChargedParticle3d.TokamakSmallNoncanonical)
        params = M.default_parameters()

        # displaced off the initial condition, which sits on a symmetry axis where several of the
        # metric derivatives vanish and would not exercise the new terms
        q = collect(float.(M.qᵢ)) .+ 0.01 .* [1.0, 2.0, 3.0, 1.0, 2.0, 3.0]

        g = zeros(6)
        N.dH(g, 0.0, N.fieldpoint(M.FIELD, 0.0, q))
        fd = [(M.hamiltonian(0.0, (h = zeros(6); h[i] = 1e-6; q .+ h), params) -
               M.hamiltonian(0.0, (h = zeros(6); h[i] = 1e-6; q .- h), params)) / 2e-6
              for i in 1:6]
        @test isapprox(g, fd; rtol = 1e-6, atol = 1e-10)

        Ω = zeros(6, 6)
        N.ω(Ω, 0.0, q, params)
        v = zeros(6)
        N.charged_particle_3d_v(v, 0.0, q, params)
        @test isapprox(v, -Ω \ g; rtol = 1e-10, atol = 1e-14)
    end
end

@safetestset "3D guiding centre: the field accessors are the derivatives they name                             " begin
    using ChargedParticleDynamics.GuidingCenter3d
    using ..StructureTestUtils
    using Test

    # `guiding_center_3d_v`/`_f` evaluate every field tensor once into a `FieldPoint` and read single
    # entries out of it, by the component names the equations are written in — `dA₁dx₂` is
    # `DA♭[1, 2]`, `dg¹¹dx₂` is `Dg♯[1, 1, 2]`. A transposed index there would surface only as a
    # silently wrong trajectory: the right-hand side would still run, and every one of these
    # equilibria integrates happily with a corrupted `∂b/∂x`. So each derivative accessor is checked
    # against a finite difference of the accessor it differentiates, on one chart of each kind.
    G = GuidingCenter3d
    V = (Val(1), Val(2), Val(3))

    for M in (G.SolovevIterXpoint,        # cylindrical
        G.TokamakMediumCartesian,   # cartesian
        G.TokamakSmallToroidal)     # toroidal
        ic = M.initial_conditions_barely_passing()
        t, q = 0.0, ic.q
        F = G.fieldpoint(M.FIELD, t, q)
        S = G.fieldpoint²(M.FIELD, t, q)
        P(x) = G.fieldpoint²(M.FIELD, t, x)

        @test isconcretetype(typeof(F))
        @test isconcretetype(typeof(S))

        fd(f, j) = central_difference(x -> f(P(x)), q, j)
        agree(a, b) = isapprox(a, b; rtol = 1e-5, atol = 1e-8)

        for i in 1:3, j in 1:3

            @test agree(G.dAᵢdxⱼ(V[i], V[j], t, F), fd(x -> G.Aᵢ(V[i], t, x), j))
            @test agree(G.dbᵢdxⱼ(V[i], V[j], t, F), fd(x -> G.bᵢ(V[i], t, x), j))

            for k in 1:3
                @test agree(G.d²Aᵢdxⱼdxₖ(V[i], V[j], V[k], t, S),
                    fd(x -> G.dAᵢdxⱼ(V[i], V[j], t, x), k))
                @test agree(G.d²bᵢdxⱼdxₖ(V[i], V[j], V[k], t, S),
                    fd(x -> G.dbᵢdxⱼ(V[i], V[j], t, x), k))
            end
        end

        # The per-component names the Hamiltonian gradients and `u` call without an index accessor.
        g = (G.g¹¹, G.g²², G.g³³)
        dg = ((G.dg¹¹dx₁, G.dg¹¹dx₂, G.dg¹¹dx₃), (G.dg²²dx₁, G.dg²²dx₂, G.dg²²dx₃),
            (G.dg³³dx₁, G.dg³³dx₂, G.dg³³dx₃))
        d²g¹¹ = ((G.d²g¹¹dx₁dx₁, G.d²g¹¹dx₁dx₂, G.d²g¹¹dx₁dx₃),
            (G.d²g¹¹dx₂dx₁, G.d²g¹¹dx₂dx₂, G.d²g¹¹dx₂dx₃),
            (G.d²g¹¹dx₃dx₁, G.d²g¹¹dx₃dx₂, G.d²g¹¹dx₃dx₃))
        dB = (G.dBdx₁, G.dBdx₂, G.dBdx₃)
        d²B = ((G.d²Bdx₁dx₁, G.d²Bdx₁dx₂, G.d²Bdx₁dx₃),
            (G.d²Bdx₂dx₁, G.d²Bdx₂dx₂, G.d²Bdx₂dx₃),
            (G.d²Bdx₃dx₁, G.d²Bdx₃dx₂, G.d²Bdx₃dx₃))

        for i in 1:3, j in 1:3

            @test agree(dg[i][j](t, F), fd(x -> g[i](t, x), j))
            @test agree(d²g¹¹[i][j](t, S), fd(x -> dg[1][i](t, x), j))
            @test agree(d²B[i][j](t, S), fd(x -> dB[i](t, x), j))
        end
        for j in 1:3
            @test agree(dB[j](t, F), fd(x -> G.B(t, x), j))
        end

        # The larger point answers everything the smaller one does, from the same tensors.
        for i in 1:3
            @test G.Aᵢ(V[i], t, S) == G.Aᵢ(V[i], t, F)
            @test G.bᵢ(V[i], t, S) == G.bᵢ(V[i], t, F)
            @test g[i](t, S) == g[i](t, F)
            @test dB[i](t, S) == dB[i](t, F)
        end
    end
end

@safetestset "3D guiding centre: second derivatives of H match finite differences                                 " begin
    using ChargedParticleDynamics.GuidingCenter3d
    using ..StructureTestUtils
    using Test

    # The metric cross-term of ∂²H/∂qᵢ∂qⱼ carried a factor of two that is only correct on the
    # diagonal, which made the Hessian asymmetric — by 65 % on the ITER Solov'ev X-point — in every
    # curvilinear coordinate system. Both the symmetry and the value are checked here.
    G = GuidingCenter3d

    for M in (G.SolovevIterXpoint,        # cylindrical
        G.TokamakMediumCartesian,   # cartesian
        G.TokamakSmallCylindrical,  # cylindrical
        G.TokamakSmallToroidal)     # toroidal
        ic = M.initial_conditions_barely_passing()
        t, q, p, par = 0.0, ic.q, ic.p, ic.params
        P(x) = G.fieldpoint²(M.FIELD, t, x)

        dH = (G.dHdq₁, G.dHdq₂, G.dHdq₃)

        for (i, j, f) in ((1, 2, G.d²Hdq₁dq₂), (1, 3, G.d²Hdq₁dq₃), (2, 3, G.d²Hdq₂dq₃),
            (1, 1, G.d²Hdq₁dq₁), (2, 2, G.d²Hdq₂dq₂), (3, 3, G.d²Hdq₃dq₃))
            ana = f(t, P(q), p, par)
            # ∂/∂qⱼ of ∂H/∂qᵢ and ∂/∂qᵢ of ∂H/∂qⱼ must both reproduce it
            @test ana ≈ central_difference(x -> dH[i](t, P(x), p, par), q, j) rtol = 1e-5
            @test ana ≈ central_difference(x -> dH[j](t, P(x), p, par), q, i) rtol = 1e-5
        end
    end
end

@safetestset "3D guiding centre: constraint derivatives match finite differences                                  " begin
    using ChargedParticleDynamics.GuidingCenter3d
    using ..StructureTestUtils
    using Test

    # The three constraints and every one of their first and second derivatives come from one indexed
    # definition each in `guiding_center_3d_constraints.jl`, so a sign or an index that is wrong is
    # wrong for all three pairs at once. This checks each of them against a central difference of the
    # derivative one order below, for every index combination, in three charts, and on the two
    # Poincaré-invariant fixtures, whose `b = e₃` makes `b₁` and `b₂` vanish. The fixtures have no
    # named initial condition, so a point on their loop serves.
    G = GuidingCenter3d
    test_point(M) =
        if isdefined(M, :initial_conditions_barely_passing)
            M.initial_conditions_barely_passing()
        else
            M.initial_conditions(0.0, M.f_loop(0.1))
        end

    for M in (G.SolovevIterXpoint,        # cylindrical
        G.TokamakMediumCartesian,   # cartesian
        G.TokamakSmallToroidal,     # toroidal
        G.SymmetricField,           # b = e₃
        G.ThetaPinchField)          # b = e₃
        ic = test_point(M)
        t, q, p = 0.0, ic.q, ic.p
        P(x) = G.fieldpoint²(M.FIELD, t, x)
        Pq = P(q)

        # A central difference at h = 1e-6 of a quantity of size s carries about `eps * s / h ≈
        # 2E-10 s` of round-off, and several entries of these Hessians vanish exactly, so the
        # comparison needs an absolute floor scaled to the problem rather than a fixed one. `s` is
        # bounded by the momentum, which runs from 1E-3 on the toroidal tokamak to 15 on ITER.
        atol = 1e-8 * max(1.0, maximum(abs, p))

        for ki in 1:3, li in 1:3

            k, l = Val(ki), Val(li)

            # ∂gᵏ/∂qₗ and ∂gᵏ/∂pₗ against a difference of gᵏ itself
            @test isapprox(G.dgᵏdqₗ(k, l, t, Pq, p),
                central_difference(x -> G.gᵏ(k, t, P(x), p), q, li); rtol = 1e-5, atol = atol)
            @test isapprox(G.dgᵏdpₗ(k, l, t, Pq, p),
                central_difference(y -> G.gᵏ(k, t, Pq, y), p, li); rtol = 1e-5, atol = atol)

            for mi in 1:3
                m = Val(mi)

                # ∂²gᵏ/∂qₗ∂qₘ has to reproduce a difference of ∂gᵏ/∂qₗ, and to be symmetric in the
                # two indices.
                @test isapprox(G.d²gᵏdqₗdqₘ(k, l, m, t, Pq, p),
                    central_difference(x -> G.dgᵏdqₗ(k, l, t, P(x), p), q, mi);
                    rtol = 1e-4, atol = atol)
                @test isapprox(G.d²gᵏdqₗdqₘ(k, l, m, t, Pq, p),
                    G.d²gᵏdqₗdqₘ(k, m, l, t, Pq, p); rtol = 1e-12, atol = 1e-14)

                # ∂²gᵏ/∂qₗ∂pₘ likewise, differentiated from either side
                @test isapprox(G.d²gᵏdqₗdpₘ(k, l, m, t, Pq, p),
                    central_difference(y -> G.dgᵏdqₗ(k, l, t, Pq, y), p, mi);
                    rtol = 1e-4, atol = atol)
                @test isapprox(G.d²gᵏdqₗdpₘ(k, l, m, t, Pq, p),
                    central_difference(x -> G.dgᵏdpₗ(k, m, t, P(x), p), q, li);
                    rtol = 1e-4, atol = atol)
            end
        end

        # `{cᵢ, cⱼ}(m) = bₘ D` with the same `D` for all three `m`. This is what lets the compact form
        # replace the singular factor `bₘ / {cᵢ, cⱼ}` by `1/D`, and it is the one step of that
        # derivation that is not pure algebra — it has to hold in every chart, which is exactly what
        # cannot be read off the paper, since Eq. (29) is written in cartesian coordinates. Checked
        # here in all three chart types.
        D = G.compact_denominator(t, Pq, p)
        for m in (Val(1), Val(2), Val(3))
            b = G.bᵢ(m, t, Pq)
            iszero(b) && continue
            @test G.bracket_cc(m, t, Pq, p) / b ≈ D rtol = 1e-10
        end

        # b·(b×v) = 0 identically, which in the paper's labelling reads b₁g² - b₂g¹ + b₃g³ = 0. It is
        # the identity the compact form uses to remove the singular factor, and it must hold whether
        # or not the constraints themselves do.
        for y in (p, p .+ 0.1)
            @test isapprox(
                G.b₁(t, Pq) * G.g₂(t, Pq, y) - G.b₂(t, Pq) * G.g₁(t, Pq, y) +
                G.b₃(t, Pq) * G.g₃(t, Pq, y),
                0.0;
                atol = 1e-12)
        end
    end
end

@safetestset "3D guiding centre: multiplier derivatives match finite differences                                  " begin
    using ChargedParticleDynamics.GuidingCenter3d
    using ..StructureTestUtils
    using Test

    # `∂λ/∂q` and `∂λ/∂p` are the only *composed* derivatives in the model — everything else is read
    # straight off the field tensors — and they are what separates `hodeproblem_canonical`
    # from `hodeproblem`. `dλ₂` carried the `∂λₒ` term with the wrong sign, which this catches.
    #
    # The test point has to sit *off* the constraint manifold. In the right-hand side both `∂λ/∂q g`
    # terms are multiplied by a constraint, so at the shipped initial condition the sign of the `∂λₒ`
    # term is invisible: `λ₁` and `λ₂` themselves are still finite and correct there, and only their
    # derivatives are wrong. Displacing `q` breaks `v × b = 0` and makes the term measurable.
    G = GuidingCenter3d

    for M in (G.SolovevIterXpoint,        # cylindrical
        G.TokamakMediumCartesian,   # cartesian
        G.TokamakSmallToroidal)     # toroidal
        ic = M.initial_conditions_barely_passing()
        t, p, par = 0.0, ic.p, ic.params
        q = ic.q .+ 0.01 .* [1.0, 2.0, 3.0]
        c = G.constraint_pair(M.default_constraints())
        S = G.fieldpoint²(M.FIELD, t, q)

        # Same reasoning as the constraint block above: a central difference at h = 1e-6 carries about
        # `eps * s / h` of round-off, and the multipliers run over several orders of magnitude between
        # these equilibria, so the floor is scaled to the problem rather than fixed.
        atol = 1e-8 * max(1.0, maximum(abs, p))

        for (λ, dλdq, dλdp) in ((G.λ₁, G.dλ₁dqⱼ, G.dλ₁dpⱼ), (G.λ₂, G.dλ₂dqⱼ, G.dλ₂dpⱼ))
            for ji in 1:3
                j = Val(ji)
                @test isapprox(dλdq(j, t, S, p, par, c),
                    central_difference(x -> λ(t, x, p, par, c), q, ji);
                    rtol = 1e-5, atol = atol)
                @test isapprox(dλdp(j, t, S, p, par, c),
                    central_difference(y -> λ(t, q, y, par, c), p, ji);
                    rtol = 1e-5, atol = atol)
            end
        end
    end
end

@safetestset "3D guiding centre: the compact form reproduces the paper's Eq. (29)                                 " begin
    using ChargedParticleDynamics.GuidingCenter3d
    using LinearAlgebra
    using Test

    # `guiding_center_3d_compact.jl` does not port Eq. (29) of Li, Zhang & Liu — that equation is
    # written with cartesian vector identities and `H = ½Σ(pᵢ-Aᵢ)²`, and eight of the thirteen
    # equilibria here are curvilinear — but derives an equivalent from the model's own objects. This
    # pins that derivation against the paper by spelling Eq. (29) out directly from the field
    # tensors, which is only possible in the cartesian charts, where the metric is the identity.
    G = GuidingCenter3d
    using ElectromagneticFields: g♯, b♭, A♭, Db♭, DA♭, DB, E♭, B

    # The field is read here straight from `ElectromagneticFields`, not through the component names
    # the right-hand side is written in, so that the check is independent of them.
    for M in (G.Dipole3d, G.QuadraticPotentials3d, G.TokamakMediumCartesian)
        ic = isdefined(M, :initial_conditions_dipole) ? M.initial_conditions_dipole() :
             isdefined(M, :initial_conditions_quadratic) ?
             M.initial_conditions_quadratic() :
             M.initial_conditions_barely_passing()
        t, q, p, par = 0.0, ic.q, ic.p, ic.params
        field = par.field

        @test g♯(field, t, q) == I

        b = Vector(b♭(field, t, q))
        A = Vector(A♭(field, t, q))
        v = p .- A
        db = Matrix(Db♭(field, t, q))
        dA = Matrix(DA♭(field, t, q))
        ∇B = Vector(DB(field, t, q))
        # `E = -∇Φ`: `hamiltonian` carries `+φ` while `dHdqᵢ` subtracts `Eᵢ`.
        ∇Φ = -Vector(E♭(field, t, q))
        ∇xb = [db[3, 2] - db[2, 3], db[1, 3] - db[3, 1], db[2, 1] - db[1, 2]]

        # The scalar denominator of Eq. (24), which `compact_denominator` reaches without ever
        # dividing by a component of b.
        D = B(field, t, q) + dot(v, ∇xb)
        @test G.compact_denominator(t, G.fieldpoint(field, t, q), p) ≈ D rtol = 1e-12

        ξ = cross(v, db * v)                                                             # Eq. (26)
        Ẋ = v .+ (ξ .+ par.μ .* cross(b, ∇B) .+ cross(b, ∇Φ)) ./ D                       # Eq. (29a)
        u = (p[3] - A[3]) / b[3]                                                         # the (g¹,g²) pair
        ṗ = -par.μ .* ∇B .- ∇Φ .+ dA' * Ẋ .+
            db' * ((u .* ξ .+ cross(v, par.μ .* ∇B .+ ∇Φ)) ./ D)                         # Eq. (29b)

        for constraints in (:g12, :parallel)
            # The two branches differ only in how they evaluate the parallel velocity and the
            # denominator, which agree on the constraint manifold the initial condition sits on.
            m = G.compact_index(constraints)
            v̄ = zeros(3)
            f̄ = zeros(3)
            G.guiding_center_3d_compact_v(v̄, t, q, p, par, m)
            G.guiding_center_3d_compact_f(f̄, t, q, p, par, m)

            @test isapprox(v̄, Ẋ; rtol = 1e-10)
            @test isapprox(f̄, ṗ; rtol = 1e-10, atol = 1e-14)
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
    cases = ((GuidingCenter4d.SymmetricField, [0.35, -0.22, 0.4, 0.3]),
        (GuidingCenter4d.TokamakMediumCylindrical, [2.4, 0.15, 0.3, 0.2]),
        (GuidingCenter4d.TokamakSmallToroidal, [0.06, 0.4, 0.3, 3e-4]))

    G = GuidingCenter4d

    for (M, q0) in cases
        t = 0.0
        q = collect(float.(q0))
        v = [0.7, -0.3, 0.45, 0.2]
        buf = zeros(4, 4)
        params = M.default_parameters()

        dϑk_dot_q(x, k) = (M.dϑ(buf, t, x, params); sum(buf[l, k] * q[l] for l in 1:4))

        for (k, ḡ) in enumerate((G.g̅₁, G.g̅₂, G.g̅₃, G.g̅₄))
            ana = ḡ(t, G.fieldpoint²(M.FIELD, t, q), v)
            num = sum(central_difference(x -> dϑk_dot_q(x, k), q, j) * v[j] for j in 1:4)
            # Absolute tolerance as well: ḡ₄ vanishes identically where b is constant, and the
            # central difference of an O(1) quantity at h = 1e-6 carries ~1e-10 of noise.
            @test isapprox(ana, num; rtol = 1e-5, atol = 1e-9)
        end
    end
end

@safetestset "4D guiding centre: the one-form and its curl take a coordinate vector and params                    " begin
    using ChargedParticleDynamics.GuidingCenter4d
    using ElectromagneticFields: A♭, b♭
    using ..StructureTestUtils
    using Test

    # `ϑ₁`…`ϑ₄`, `β₁`…`β₃` and `ϑ(t, q, params, k)` are exported, and a caller holds a coordinate
    # vector rather than a `FieldPoint`. The expected values do not come from the code under test:
    # the one-form is the field's own `A♭ + u b♭`, and its curl `β = ∇ × ϑ` is taken by central
    # differences of the one-form.
    G = GuidingCenter4d

    for M in (G.SolovevIterXpoint, G.TokamakSmallCartesian, G.TokamakSmallToroidal)
        ic = M.initial_conditions_barely_passing()
        t, q, par = 0.0, ic.q, ic.params
        x = q[1:3]
        θ = A♭(par.field, t, x) .+ q[4] .* b♭(par.field, t, x)
        ϑs = (G.ϑ₁, G.ϑ₂, G.ϑ₃)

        for k in 1:3
            @test ϑs[k](t, q, par) ≈ θ[k] rtol = 1e-14
            @test G.ϑ(t, q, par, k) == ϑs[k](t, q, par)
        end
        @test G.ϑ₄(t, q, par) == 0
        @test G.ϑ(t, q, par, 4) == 0

        d(k, j) = central_difference(y -> ϑs[k](t, y, par), q, j)
        atol = 1e-8 * max(1.0, maximum(abs, θ))
        @test isapprox(G.β₁(t, q, par), d(3, 2) - d(2, 3); rtol = 1e-6, atol = atol)
        @test isapprox(G.β₂(t, q, par), d(1, 3) - d(3, 1); rtol = 1e-6, atol = atol)
        @test isapprox(G.β₃(t, q, par), d(2, 1) - d(1, 2); rtol = 1e-6, atol = atol)
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

        # the equations each of the three modules forwards to
        Fp = PauliParticle3d
        Fc = ChargedParticle3d.Canonical
        Fn = ChargedParticle3d.Noncanonical

        v = [1.3e-3, -4.0e-4, 7.0e-4]
        λ3 = [0.31, -0.17, 0.42]
        λ6 = [0.31, -0.17, 0.42, 0.23, -0.51, 0.11]

        # (module, g, state, oneform, nz, λ, params)
        pp = (field = Mp.FIELD, μ = 2.31e-6)
        pc = Mc.default_parameters()
        pn = Mn.default_parameters()
        cases = (
            (Mp, Fp.pauli_particle_3d_iode_g, [1.05, 0.0, 0.0],
                (z, vv) -> (
                    θ = zeros(3); Fp.pauli_particle_3d_iode_ϑ(θ, 0.0, z, vv, pp); θ),
                3, λ3, pp),
            (Mc, Fc.charged_particle_3d_iode_g, [2.5, 0.0, 0.0],
                (z, vv) -> (
                    θ = zeros(3); Fc.charged_particle_3d_iode_ϑ(θ, 0.0, z, vv, pc); θ),
                3, λ3, pc),
            (Mn, Fn.charged_particle_3d_iode_g, collect(float.(Mn.qᵢ)),
                (z, vv) -> (θ = zeros(6); Mn.ϑ(θ, 0.0, z, pn); θ), 6, λ6, pn)
        )

        for (M, g, z, oneform, nz, λ, par) in cases
            t = 0.0

            ga = zeros(nz)
            g(ga, t, z, v, λ, par)

            # gᵢ = ∂/∂zⁱ (ϑ ⋅ λ), by central differences
            for i in 1:nz
                num = central_difference(x -> sum(oneform(x, v) .* λ), z, i)
                @test isapprox(ga[i], num; rtol = 1e-5, atol = 1e-9)
            end

            # and it must be linear in λ, which the finite difference cannot see
            g2 = zeros(nz)
            g(g2, t, z, v, 2 .* λ, par)
            g0 = zeros(nz)
            g(g0, t, z, v, zero(λ), par)

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

    # `toroidal_momentum` is ϑ₃, the covariant φ-component of the one-form, and not R ϑ₃, which is
    # not conserved: on the small cylindrical tokamak the relative variation over 10³ time units is
    # 2e-14 for ϑ₃ and 4e-3 for R ϑ₃. In cartesian coordinates the third coordinate is z rather
    # than an angle, so neither is right there and the momentum is the generator of rotation about
    # the z-axis instead.
    #
    # The tolerances below are far tighter than R ϑ₃ can reach, which is the point of the test;
    # they are loose relative to what the integrator achieves.
    cases = (
        (GuidingCenter4d.TokamakSmallCylindrical,
            (timestep = 10.0, timespan = (0.0, 1e3)), 1e-10),
        (GuidingCenter4d.TokamakSmallToroidal,
            (timestep = 10.0, timespan = (0.0, 1e3)), 1e-10),
        (GuidingCenter4d.TokamakSmallCartesian,
            (timestep = 10.0, timespan = (0.0, 1e3)), 1e-8),
        (GuidingCenter4d.SolovevIterXpoint, (timestep = 1.0, timespan = (0.0, 1e2)), 1e-10))

    for (M, kw, tol) in cases
        prob = M.odeproblem(M.initial_conditions_barely_passing(); kw...)
        sol = integrate(prob, Gauss(2))
        _, perr = M.compute_toroidal_momentum_error(sol)
        @test maximum(abs(perr[i]) for i in eachindex(perr)) < tol
    end
end

@safetestset "A formula written in one chart refuses a field in another                                           " begin
    using ChargedParticleDynamics
    using ChargedParticleDynamics.GuidingCenter4d
    using Test

    # The right-hand sides take any field whose chart is orthogonal. Two things are narrower, and
    # both refuse rather than answer wrongly: `toroidal_momentum`, whose formula each module writes
    # in its own chart, and the metric accessors of `FieldPoints`, which read the diagonal only.
    G = GuidingCenter4d
    cart, cyl = G.TokamakSmallCartesian, G.TokamakSmallCylindrical
    ic = cyl.initial_conditions_barely_passing()
    foreign = (field = cart.FIELD, μ = ic.params.μ)

    @test cyl.toroidal_momentum(0.0, ic.q, ic.params) isa Float64
    @test_throws ArgumentError cyl.toroidal_momentum(0.0, ic.q, foreign)
    @test_throws ArgumentError cart.toroidal_momentum(0.0, ic.q, ic.params)

    # The same kind of equilibrium with other parameters is the same chart; `GuidingCenter3d`'s
    # `compute_toroidal_momentum` goes through this check.
    @test ChargedParticleDynamics.check_chart(
        cyl.FIELD, G.TokamakMediumCylindrical.FIELD) === nothing
    @test_throws ArgumentError ChargedParticleDynamics.check_chart(cart.FIELD, cyl.FIELD)

    FP = ChargedParticleDynamics.FieldPoints
    x = [1.0, 0.0, 0.0]
    diagonal = [1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 3.0]
    skew = [1.0 0.1 0.0; 0.1 2.0 0.0; 0.0 0.0 3.0]
    @test FP.check_orthogonal((g♭ = diagonal, g♯ = diagonal), x) === nothing
    @test FP.check_orthogonal((B = 1.0,), x) === nothing
    @test_throws ArgumentError FP.check_orthogonal((g♭ = skew,), x)
    @test_throws ArgumentError FP.check_orthogonal((g♯ = skew,), x)
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
    # model is checked, now that the integration tests in `guiding_center_3d_tests.jl` run 100. The
    # iteration cap is there because a handful of steps otherwise run to the library's limit of 1000
    # and accounted for most of the block's runtime (48 s of 49 s on the Solov'ev case).
    #
    # `f_abstol` has to be stated because the library default is below this problem's round-off
    # floor, not because passing options would otherwise disturb the defaults: since
    # `GeometricIntegratorsBase` 0.5.1 the caller's options are *merged* into the method's
    # `default_options` rather than replacing them, and that default is
    # `f_abstol = max(8, solversize(method, problem)) * eps(datatype(problem))` — `2.7E-15` for
    # `PartitionedGauss(2)` on this three-component model (`solversize = 12`). Below `‖ϑ‖ eps`, so
    # still not reachable here.
    #
    # It used to restate the library default of `8eps() = 1.8E-15` so as to measure the conservation
    # an unconfigured `integrate` delivers. That was the wrong criterion to hold this block to: the
    # package's own measurements put the round-off floor of the ITER-scale residual at
    # `‖ϑ‖ eps ≈ 3.5E-15`, so `8eps()` is unreachable there — `quiet_solver_warnings.jl` and
    # `study_solver_tolerances.jl` both say so — and four steps of the Solov'ev X-point ran to the
    # iteration cap chasing it. `1E-12` is what every other test module, script and example in the
    # package uses, for the reasons set out on `options` in `guiding_center_3d_tests.jl`.
    #
    # Nothing is given up by stating it. The library defaults are still exercised, by the calls that
    # genuinely pass no options at all — `integrate(prob, Gauss(2))` earlier in this file, and the
    # three in `plots_tests.jl` — which are silent, and which `runtests.jl` asserts to be.
    # The energy error here is identical at both tolerances (5.388E-15 on the X-point,
    # 3.226E-10 on the medium tokamak), the medium tokamak is bit-identical throughout, and the
    # X-point's constraint errors move by under a quarter of one order against five orders of slack
    # in the bound below. What it buys is four fewer capped steps, four fewer suppressed warnings,
    # and a block that no longer asks the solver for something it cannot deliver.
    options = (f_abstol = 1E-12, max_iterations = 50, warn_iterations = 50)

    for M in (GuidingCenter3d.SolovevIterXpoint, GuidingCenter3d.TokamakMediumCartesian)
        prob = M.hodeproblem(M.initial_conditions_barely_passing(); timestep = 0.1, timespan = (
            0.0, 1e2))
        # No `initialguess`: `PartitionedGauss` defaults to `HermiteExtrapolation`, and the
        # `MidpointExtrapolation(5)` this used to request was 70% of the block's runtime while
        # leaving the trajectory bit-identical. See `guiding_center_3d_tests.jl`.
        #
        # At the `8eps()` this block used to request, four of these thousand steps on the Solov'ev
        # X-point ran to the iteration cap: steps 261, 262, 263 and 271, where the orbit crosses
        # `b₁ = 0` — `b₁` runs 3.6E-5, 1.0E-5, -1.5E-5 against a typical -5.9E-3 — so `λₒ ∝ b₁`
        # collapses and the residual is badly scaled. They did not converge at *any* iteration
        # budget: 50, 100, 200 and 1000 all terminated at the cap and all four gave a bit-identical
        # final state, while the mean iteration count over the run rose from 1.25 to 5.05. That is
        # why the cap is 50 and not higher, and why the fix was the tolerance rather than the budget.
        # At `1E-12` nothing caps and every step converges in one iteration.
        sol = integrate(prob, PartitionedGauss(2); options...)

        _, eerr = M.compute_energy_error(sol)
        @test mx(eerr) < 1e-8

        # All three components of `v × b` vanish along the flow, not only the two the problem's
        # constraint pair retains, so all three are asserted on — but not to the same bound. The
        # identity `b₁g² - b₂g¹ + b₃g³ = 0` holds pointwise, so it determines the omitted constraint
        # from the retained two divided by `bₘ`, the component of `b` the pair's multipliers divide by:
        #
        #     |g_omitted| ≲ max |g_retained| / min|bₘ|
        #
        # with the minimum taken along the orbit, not at the initial condition. Both equilibria here
        # cross the surface where `bₘ` is small — `min|b₁|` is 1E-5 on the ITER Solov'ev X-point and
        # 3.5E-4 on the medium tokamak — and the omitted `g²` duly comes out three orders of magnitude
        # larger than the retained pair on the latter. It is the sharpest argument for picking the pair
        # whose `bₘ` is largest, and the reason `compute_constraints` reports all three.
        c = M.compute_constraints(sol)
        gs = (c.g₁, c.g₂, c.g₃)
        G = GuidingCenter3d
        bs = map(b -> (t, x) -> b(t, G.fieldpoint(M.FIELD, t, x)), (G.b₁, G.b₂, G.b₃))

        retained = map(G.unval, G.constraint_pair(M.default_constraints()))
        omitted = only(setdiff(1:3, retained))

        # `compact_index` is the component of `b` the pair divides by, which is not `omitted`: the
        # labelling of the `gᵏ` is not the antisymmetric one, so `(g³, g¹)` omits `g²` but divides by
        # `b₁`.
        m = G.unval(G.compact_index(M.default_constraints()))
        bmin = minimum(abs(bs[m](sol.t[i], sol.q[i])) for i in eachindex(sol.t))

        for k in retained
            @test mx(gs[k]) < 1e-8
        end

        # The bound is anchored to what the retained pair actually achieved rather than to the 1E-8 it
        # is merely required to stay under. Dividing the requirement by `bmin` gives 9E-4 on the ITER
        # Solov'ev X-point against an observed 1.1E-12 — nine orders of slack, which asserts nothing.
        # Anchoring closes that to five: the retained pair does 1.4E-13 there, five orders better than
        # required.
        #
        # The factor is for the relation above being an inequality with no constant attached. Measured,
        # the two sides are within 0.3 on the medium tokamak and 8E-5 on the X-point, so 10 leaves the
        # tighter of the two about 35× of headroom — enough for a solver or platform difference, little
        # enough that a real regression in the omitted constraint trips it. The floor is for the `:g12`
        # equilibria, where `bₘ` is the *large* component of `b` and the omitted constraint comes out
        # better conserved than the retained pair, so the ratio bound goes below round-off.
        retained_max = maximum(mx(gs[k]) for k in retained)
        @test mx(gs[omitted]) < max(1e-12, 10 * retained_max / bmin)

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
        M.ω(Ω, 0.0, q, M.default_parameters())
        @test Ω ≈ -transpose(Ω)
    end

    let M = ChargedParticle3d.TokamakSmallNoncanonical
        q = collect(float.(M.qᵢ))
        Ω = zeros(length(q), length(q))
        ChargedParticle3d.Noncanonical.ω(Ω, 0.0, q, M.default_parameters())
        @test Ω ≈ -transpose(Ω)
    end

    let M = ChargedParticle3d.SolovevIterXpoint, C = ChargedParticle3d.Canonical
        params = M.default_parameters()
        q = [2.5, 0.0, 0.0]
        v = [1.3e-3, -4.0e-4, 7.0e-4]
        Ω = zeros(length(q), length(q))
        C.ω(Ω, 0.0, q, v, params)
        @test Ω ≈ -transpose(Ω)

        # Ωᵢⱼ = ∂ϑᵢ/∂qʲ - ∂ϑⱼ/∂qⁱ, against central differences of the one-form itself.
        ϑi(x, i) = (θ = zeros(3); C.charged_particle_3d_iode_ϑ(θ, 0.0, x, v, params); θ[i])
        for i in 1:3, j in 1:3

            num = central_difference(x -> ϑi(x, i), q, j) -
                  central_difference(x -> ϑi(x, j), q, i)
            @test isapprox(Ω[i, j], num; rtol = 1e-5, atol = 1e-9)
        end
    end
end
