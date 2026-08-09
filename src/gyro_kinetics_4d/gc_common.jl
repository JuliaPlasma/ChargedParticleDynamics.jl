
using Parameters

export hamiltonian, ϑ, ω, ωabs, β, γ, v


@inline function u(t, q)
    q[4]
end


@doc raw"""
The equilibrium guiding centre Hamiltonian,

```math
H_{0} = \tfrac{1}{2} u^{2} + \mu \vert B \vert + \varphi ,
```

matching `GuidingCenter4d`. The notes this module follows give the zero Larmor radius
Hamiltonian ``H^{\mathrm{zlr}}``, which carries ``A_{\parallel}`` and
``\langle \tilde{A}_{\parallel} \rangle^{2}`` terms on top of this and would change ``β`` and ``γ``
as well; only ``H_{0}`` is implemented. See `TODO.md`.
"""
function hamiltonian(t, q, params)
    @unpack μ = params
    0.5 * u(t,q)^2 + μ*B(t,q) + φ(t,q)
end


# ∂φ/∂xᵢ = -Eᵢ, hence the minus sign; ∂²φ/∂xᵢ∂xⱼ = -∂Eⱼ/∂xᵢ below.
dHdx₁(t,q,μ) = μ * dBdx₁(t,q) - E₁(t,q)
dHdx₂(t,q,μ) = μ * dBdx₂(t,q) - E₂(t,q)
dHdx₃(t,q,μ) = μ * dBdx₃(t,q) - E₃(t,q)
dHdx₄(t,q,μ) = u(t,q)

function dH(dH, t, q, params)
    @unpack μ = params
    dH[1] = dHdx₁(t,q,μ)
    dH[2] = dHdx₂(t,q,μ)
    dH[3] = dHdx₃(t,q,μ)
    dH[4] = dHdx₄(t,q,μ)
    nothing
end


d²Hdx₁dx₁(t, q, μ) = μ * d²Bdx₁dx₁(t,q) - dE₁dx₁(t,q)
d²Hdx₁dx₂(t, q, μ) = μ * d²Bdx₁dx₂(t,q) - dE₂dx₁(t,q)
d²Hdx₁dx₃(t, q, μ) = μ * d²Bdx₁dx₃(t,q) - dE₃dx₁(t,q)
d²Hdx₁dx₄(t, q, μ) = zero(eltype(q))
d²Hdx₂dx₁(t, q, μ) = μ * d²Bdx₂dx₁(t,q) - dE₁dx₂(t,q)
d²Hdx₂dx₂(t, q, μ) = μ * d²Bdx₂dx₂(t,q) - dE₂dx₂(t,q)
d²Hdx₂dx₃(t, q, μ) = μ * d²Bdx₂dx₃(t,q) - dE₃dx₂(t,q)
d²Hdx₂dx₄(t, q, μ) = zero(eltype(q))
d²Hdx₃dx₁(t, q, μ) = μ * d²Bdx₃dx₁(t,q) - dE₁dx₃(t,q)
d²Hdx₃dx₂(t, q, μ) = μ * d²Bdx₃dx₂(t,q) - dE₂dx₃(t,q)
d²Hdx₃dx₃(t, q, μ) = μ * d²Bdx₃dx₃(t,q) - dE₃dx₃(t,q)
d²Hdx₃dx₄(t, q, μ) = zero(eltype(q))
d²Hdx₄dx₁(t, q, μ) = zero(eltype(q))
d²Hdx₄dx₂(t, q, μ) = zero(eltype(q))
d²Hdx₄dx₃(t, q, μ) = zero(eltype(q))
d²Hdx₄dx₄(t, q, μ) = one(eltype(q))


ϑ₁(t,q) = A₁(t,q) + u(t,q) * b₁(t,q)
ϑ₂(t,q) = A₂(t,q) + u(t,q) * b₂(t,q)
ϑ₃(t,q) = A₃(t,q) + u(t,q) * b₃(t,q)
ϑ₄(t,q) = zero(eltype(q))

function ϑ(ϑ, t, q, params)
    ϑ[1] = ϑ₁(t,q)
    ϑ[2] = ϑ₂(t,q)
    ϑ[3] = ϑ₃(t,q)
    ϑ[4] = ϑ₄(t,q)
    nothing
end


dϑ₁dx₁(t,q) = dA₁dx₁(t,q) + u(t,q) * db₁dx₁(t,q)
dϑ₁dx₂(t,q) = dA₁dx₂(t,q) + u(t,q) * db₁dx₂(t,q)
dϑ₁dx₃(t,q) = dA₁dx₃(t,q) + u(t,q) * db₁dx₃(t,q)
dϑ₁dx₄(t,q) = b₁(t,q)
dϑ₂dx₁(t,q) = dA₂dx₁(t,q) + u(t,q) * db₂dx₁(t,q)
dϑ₂dx₂(t,q) = dA₂dx₂(t,q) + u(t,q) * db₂dx₂(t,q)
dϑ₂dx₃(t,q) = dA₂dx₃(t,q) + u(t,q) * db₂dx₃(t,q)
dϑ₂dx₄(t,q) = b₂(t,q)
dϑ₃dx₁(t,q) = dA₃dx₁(t,q) + u(t,q) * db₃dx₁(t,q)
dϑ₃dx₂(t,q) = dA₃dx₂(t,q) + u(t,q) * db₃dx₂(t,q)
dϑ₃dx₃(t,q) = dA₃dx₃(t,q) + u(t,q) * db₃dx₃(t,q)
dϑ₃dx₄(t,q) = b₃(t,q)
dϑ₄dx₁(t,q) = zero(eltype(q))
dϑ₄dx₂(t,q) = zero(eltype(q))
dϑ₄dx₃(t,q) = zero(eltype(q))
dϑ₄dx₄(t,q) = zero(eltype(q))


d²ϑ₁dx₁dx₁(t,q) = d²A₁dx₁dx₁(t,q) + u(t,q) * d²b₁dx₁dx₁(t,q)
d²ϑ₁dx₁dx₂(t,q) = d²A₁dx₁dx₂(t,q) + u(t,q) * d²b₁dx₁dx₂(t,q)
d²ϑ₁dx₁dx₃(t,q) = d²A₁dx₁dx₃(t,q) + u(t,q) * d²b₁dx₁dx₃(t,q)
d²ϑ₁dx₁dx₄(t,q) = db₁dx₁(t,q)

d²ϑ₁dx₂dx₁(t,q) = d²A₁dx₂dx₁(t,q) + u(t,q) * d²b₁dx₂dx₁(t,q)
d²ϑ₁dx₂dx₂(t,q) = d²A₁dx₂dx₂(t,q) + u(t,q) * d²b₁dx₂dx₂(t,q)
d²ϑ₁dx₂dx₃(t,q) = d²A₁dx₂dx₃(t,q) + u(t,q) * d²b₁dx₂dx₃(t,q)
d²ϑ₁dx₂dx₄(t,q) = db₁dx₂(t,q)

d²ϑ₁dx₃dx₁(t,q) = d²A₁dx₃dx₁(t,q) + u(t,q) * d²b₁dx₃dx₁(t,q)
d²ϑ₁dx₃dx₂(t,q) = d²A₁dx₃dx₂(t,q) + u(t,q) * d²b₁dx₃dx₂(t,q)
d²ϑ₁dx₃dx₃(t,q) = d²A₁dx₃dx₃(t,q) + u(t,q) * d²b₁dx₃dx₃(t,q)
d²ϑ₁dx₃dx₄(t,q) = db₁dx₃(t,q)

d²ϑ₁dx₄dx₁(t,q) = db₁dx₁(t,q)
d²ϑ₁dx₄dx₂(t,q) = db₁dx₂(t,q)
d²ϑ₁dx₄dx₃(t,q) = db₁dx₃(t,q)
d²ϑ₁dx₄dx₄(t,q) = zero(eltype(q))


d²ϑ₂dx₁dx₁(t,q) = d²A₂dx₁dx₁(t,q) + u(t,q) * d²b₂dx₁dx₁(t,q)
d²ϑ₂dx₁dx₂(t,q) = d²A₂dx₁dx₂(t,q) + u(t,q) * d²b₂dx₁dx₂(t,q)
d²ϑ₂dx₁dx₃(t,q) = d²A₂dx₁dx₃(t,q) + u(t,q) * d²b₂dx₁dx₃(t,q)
d²ϑ₂dx₁dx₄(t,q) = db₂dx₁(t,q)

d²ϑ₂dx₂dx₁(t,q) = d²A₂dx₂dx₁(t,q) + u(t,q) * d²b₂dx₂dx₁(t,q)
d²ϑ₂dx₂dx₂(t,q) = d²A₂dx₂dx₂(t,q) + u(t,q) * d²b₂dx₂dx₂(t,q)
d²ϑ₂dx₂dx₃(t,q) = d²A₂dx₂dx₃(t,q) + u(t,q) * d²b₂dx₂dx₃(t,q)
d²ϑ₂dx₂dx₄(t,q) = db₂dx₂(t,q)

d²ϑ₂dx₃dx₁(t,q) = d²A₂dx₃dx₁(t,q) + u(t,q) * d²b₂dx₃dx₁(t,q)
d²ϑ₂dx₃dx₂(t,q) = d²A₂dx₃dx₂(t,q) + u(t,q) * d²b₂dx₃dx₂(t,q)
d²ϑ₂dx₃dx₃(t,q) = d²A₂dx₃dx₃(t,q) + u(t,q) * d²b₂dx₃dx₃(t,q)
d²ϑ₂dx₃dx₄(t,q) = db₂dx₃(t,q)

d²ϑ₂dx₄dx₁(t,q) = db₂dx₁(t,q)
d²ϑ₂dx₄dx₂(t,q) = db₂dx₂(t,q)
d²ϑ₂dx₄dx₃(t,q) = db₂dx₃(t,q)
d²ϑ₂dx₄dx₄(t,q) = zero(eltype(q))


d²ϑ₃dx₁dx₁(t,q) = d²A₃dx₁dx₁(t,q) + u(t,q) * d²b₃dx₁dx₁(t,q)
d²ϑ₃dx₁dx₂(t,q) = d²A₃dx₁dx₂(t,q) + u(t,q) * d²b₃dx₁dx₂(t,q)
d²ϑ₃dx₁dx₃(t,q) = d²A₃dx₁dx₃(t,q) + u(t,q) * d²b₃dx₁dx₃(t,q)
d²ϑ₃dx₁dx₄(t,q) = db₃dx₁(t,q)

d²ϑ₃dx₂dx₁(t,q) = d²A₃dx₂dx₁(t,q) + u(t,q) * d²b₃dx₂dx₁(t,q)
d²ϑ₃dx₂dx₂(t,q) = d²A₃dx₂dx₂(t,q) + u(t,q) * d²b₃dx₂dx₂(t,q)
d²ϑ₃dx₂dx₃(t,q) = d²A₃dx₂dx₃(t,q) + u(t,q) * d²b₃dx₂dx₃(t,q)
d²ϑ₃dx₂dx₄(t,q) = db₃dx₂(t,q)

d²ϑ₃dx₃dx₁(t,q) = d²A₃dx₃dx₁(t,q) + u(t,q) * d²b₃dx₃dx₁(t,q)
d²ϑ₃dx₃dx₂(t,q) = d²A₃dx₃dx₂(t,q) + u(t,q) * d²b₃dx₃dx₂(t,q)
d²ϑ₃dx₃dx₃(t,q) = d²A₃dx₃dx₃(t,q) + u(t,q) * d²b₃dx₃dx₃(t,q)
d²ϑ₃dx₃dx₄(t,q) = db₃dx₃(t,q)

d²ϑ₃dx₄dx₁(t,q) = db₃dx₁(t,q)
d²ϑ₃dx₄dx₂(t,q) = db₃dx₂(t,q)
d²ϑ₃dx₄dx₃(t,q) = db₃dx₃(t,q)
d²ϑ₃dx₄dx₄(t,q) = zero(eltype(q))


d²ϑ₄dx₁dx₁(t,q) = zero(eltype(q))
d²ϑ₄dx₁dx₂(t,q) = zero(eltype(q))
d²ϑ₄dx₁dx₃(t,q) = zero(eltype(q))
d²ϑ₄dx₁dx₄(t,q) = zero(eltype(q))

d²ϑ₄dx₂dx₁(t,q) = zero(eltype(q))
d²ϑ₄dx₂dx₂(t,q) = zero(eltype(q))
d²ϑ₄dx₂dx₃(t,q) = zero(eltype(q))
d²ϑ₄dx₂dx₄(t,q) = zero(eltype(q))

d²ϑ₄dx₃dx₁(t,q) = zero(eltype(q))
d²ϑ₄dx₃dx₂(t,q) = zero(eltype(q))
d²ϑ₄dx₃dx₃(t,q) = zero(eltype(q))
d²ϑ₄dx₃dx₄(t,q) = zero(eltype(q))

d²ϑ₄dx₄dx₁(t,q) = zero(eltype(q))
d²ϑ₄dx₄dx₂(t,q) = zero(eltype(q))
d²ϑ₄dx₄dx₃(t,q) = zero(eltype(q))
d²ϑ₄dx₄dx₄(t,q) = zero(eltype(q))



ω₁(t,q) = dϑ₃dx₂(t,q) - dϑ₂dx₃(t,q)
ω₂(t,q) = dϑ₁dx₃(t,q) - dϑ₃dx₁(t,q)
ω₃(t,q) = dϑ₂dx₁(t,q) - dϑ₁dx₂(t,q)

function ω(ω::Vector, t, q, params=NamedTuple())
    ω[1] = ω₁(t,q)
    ω[2] = ω₂(t,q)
    ω[3] = ω₃(t,q)
    nothing
end

function ω(Ω::Matrix, t, q, params=NamedTuple())
    Ω[1,1] = 0
    Ω[1,2] = dϑ₁dx₂(t,q) - dϑ₂dx₁(t,q)
    Ω[1,3] = dϑ₁dx₃(t,q) - dϑ₃dx₁(t,q)
    Ω[1,4] = dϑ₁dx₄(t,q) - dϑ₄dx₁(t,q)

    Ω[2,1] = dϑ₂dx₁(t,q) - dϑ₁dx₂(t,q)
    Ω[2,2] = 0
    Ω[2,3] = dϑ₂dx₃(t,q) - dϑ₃dx₂(t,q)
    Ω[2,4] = dϑ₂dx₄(t,q) - dϑ₄dx₂(t,q)

    Ω[3,1] = dϑ₃dx₁(t,q) - dϑ₁dx₃(t,q)
    Ω[3,2] = dϑ₃dx₂(t,q) - dϑ₂dx₃(t,q)
    Ω[3,3] = 0
    Ω[3,4] = dϑ₃dx₄(t,q) - dϑ₄dx₃(t,q)

    Ω[4,1] = dϑ₄dx₁(t,q) - dϑ₁dx₄(t,q)
    Ω[4,2] = dϑ₄dx₂(t,q) - dϑ₂dx₄(t,q)
    Ω[4,3] = dϑ₄dx₃(t,q) - dϑ₃dx₄(t,q)
    Ω[4,4] = 0

    nothing
end

@doc raw"""
    ωabs(t, q, params)

The factor by which this module's vector field is rescaled relative to the 4D guiding centre one,
``v_{gk} = \omega_{abs} \, v_{gc}``, so that the independent variable is the rescaled time with
``dt = \omega_{abs} \, ds``. It is the phasespace Jacobian ``\sqrt{\det \Omega} = J \, B^{\star}_{\parallel}``:
the contraction ``\omega \cdot \partial\vartheta/\partial u`` of the curl of the one-form against
``b``, taken with the chart's orientation so that it is **positive in every chart**.

!!! note "Why `orientation()` appears here"
    `ω₁ = ∂₂ϑ₃ - ∂₃ϑ₂` is the *coordinate* curl, which is `det(DF)` times the contravariant one —
    the signed Jacobian, not the volume element `J`. The bare contraction is therefore
    `det(DF) · B*∥ = orientation() · J · B*∥`, and in the four left-handed charts of
    `ElectromagneticFields` — the cylindrical, the two toroidal and every `Solovev*` other than
    `SolovevSymmetric` — it comes out **negative** while the physical ``B^{\star}_{\parallel}`` is
    positive. Multiplying by `orientation()` recovers the Liouville density `√det Ω = |Pf(Ω)|`,
    which is non-negative by construction — as the name of this function says it should be.

    `orientation()` is generated into this module by the `@code` call at its top, alongside `J`,
    `DF` and the metric: `ElectromagneticFields` 0.7.1 emits the sign its own generator already
    reads to build `det DF`. It is the one generated function taking no arguments, since the
    handedness of a chart depends on neither `t` nor `q`. Nothing here declares it — asking the
    field for its own orientation is what keeps the two from drifting apart.

    That sign is not free to drop. The factor is the phasespace Jacobian, absorbed into the
    distribution function as `f̃ = ωabs f`; a distribution rescaled by a negative number is not one.
    Carried into the vector field, it made the model **chart-dependent**: the small tokamak's
    cartesian chart circulated the particle one way round the torus and its cylindrical chart the
    other, for the same physical initial condition. `test/gyro_kinetics_4d_tests.jl` asserts that it
    no longer does.

    The magnitude stays chart-dependent through `J`, which is correct for a density and is why the
    three charts of the small tokamak give `0.9511`, `0.9986` and `0.0499` for the same `B*∥ = 0.9511`.
    Only the sign was wrong.

    `orientation()` multiplies here and in [`v`](@ref) and nowhere else. `ω`, `Ω`, [`β`](@ref) and
    [`γ`](@ref) stay the plain coordinate objects they are — `Ω = dϑ` in particular must not be
    touched — so the identity `v_gk = ωabs · v_gc` holds exactly with both factors oriented.

    Until `ElectromagneticFields` 0.7.0 none of this was visible: `b` was itself reversed in exactly
    those charts, so the two sign errors cancelled and the bare contraction came out positive
    everywhere.
"""
function ωabs(t, q, params=NamedTuple())
    orientation() * (ω₁(t,q) * dϑ₁dx₄(t,q) + ω₂(t,q) * dϑ₂dx₄(t,q) + ω₃(t,q) * dϑ₃dx₄(t,q))
end


dω₁dx₁(t,q) = d²ϑ₃dx₂dx₁(t,q) - d²ϑ₂dx₃dx₁(t,q)
dω₁dx₂(t,q) = d²ϑ₃dx₂dx₂(t,q) - d²ϑ₂dx₃dx₂(t,q)
dω₁dx₃(t,q) = d²ϑ₃dx₂dx₃(t,q) - d²ϑ₂dx₃dx₃(t,q)
dω₁dx₄(t,q) = d²ϑ₃dx₂dx₄(t,q) - d²ϑ₂dx₃dx₄(t,q)
dω₂dx₁(t,q) = d²ϑ₁dx₃dx₁(t,q) - d²ϑ₃dx₁dx₁(t,q)
dω₂dx₂(t,q) = d²ϑ₁dx₃dx₂(t,q) - d²ϑ₃dx₁dx₂(t,q)
dω₂dx₃(t,q) = d²ϑ₁dx₃dx₃(t,q) - d²ϑ₃dx₁dx₃(t,q)
dω₂dx₄(t,q) = d²ϑ₁dx₃dx₄(t,q) - d²ϑ₃dx₁dx₄(t,q)
dω₃dx₁(t,q) = d²ϑ₂dx₁dx₁(t,q) - d²ϑ₁dx₂dx₁(t,q)
dω₃dx₂(t,q) = d²ϑ₂dx₁dx₂(t,q) - d²ϑ₁dx₂dx₂(t,q)
dω₃dx₃(t,q) = d²ϑ₂dx₁dx₃(t,q) - d²ϑ₁dx₂dx₃(t,q)
dω₃dx₄(t,q) = d²ϑ₂dx₁dx₄(t,q) - d²ϑ₁dx₂dx₄(t,q)


β₁(t,q,μ) = dHdx₄(t,q,μ) * ϑ₁(t,q)
β₂(t,q,μ) = dHdx₄(t,q,μ) * ϑ₂(t,q)
β₃(t,q,μ) = dHdx₄(t,q,μ) * ϑ₃(t,q)

@doc raw"""
    β(β, t, q, params)

The first of the two vector potentials of the volume-preserving formulation,

```math
\beta = \frac{1}{m} \frac{\partial H}{\partial u} A^{\star} ,
```

which reduces to ``\beta = u \, A^{\star}`` for the equilibrium Hamiltonian implemented here.
Together with [`γ`](@ref) it writes the equations of motion in the manifestly divergence-free form
``\dot{R} = \partial_{u} \gamma + \nabla \times \beta``, ``\dot{u} = - \nabla \cdot \gamma``.
"""
function β(β, t, q, params)
    @unpack μ = params
    β[1] = β₁(t,q,μ)
    β[2] = β₂(t,q,μ)
    β[3] = β₃(t,q,μ)
    nothing
end


γ₁(t,q,μ) = ϑ₂(t,q) * dHdx₃(t,q,μ) - dHdx₂(t,q,μ) * ϑ₃(t,q)
γ₂(t,q,μ) = ϑ₃(t,q) * dHdx₁(t,q,μ) - dHdx₃(t,q,μ) * ϑ₁(t,q)
γ₃(t,q,μ) = ϑ₁(t,q) * dHdx₂(t,q,μ) - dHdx₁(t,q,μ) * ϑ₂(t,q)

@doc raw"""
    γ(γ, t, q, params)

The second vector potential of the volume-preserving formulation,

```math
\gamma = \frac{1}{m} A^{\star} \times \nabla H ,
```

which reduces to ``\gamma = (\mu/m) \, A^{\star} \times \nabla B_{0}`` for the equilibrium
Hamiltonian implemented here. See [`β`](@ref).
"""
function γ(γ, t, q, params)
    @unpack μ = params
    γ[1] = γ₁(t,q,μ)
    γ[2] = γ₂(t,q,μ)
    γ[3] = γ₃(t,q,μ)
    nothing
end


dβ₁dx₂(t,q,μ) = d²Hdx₄dx₂(t,q,μ) * ϑ₁(t,q) + dHdx₄(t,q,μ) * dϑ₁dx₂(t,q)
dβ₁dx₁(t,q,μ) = d²Hdx₄dx₁(t,q,μ) * ϑ₁(t,q) + dHdx₄(t,q,μ) * dϑ₁dx₁(t,q)
dβ₁dx₃(t,q,μ) = d²Hdx₄dx₃(t,q,μ) * ϑ₁(t,q) + dHdx₄(t,q,μ) * dϑ₁dx₃(t,q)
dβ₁dx₄(t,q,μ) = d²Hdx₄dx₄(t,q,μ) * ϑ₁(t,q) + dHdx₄(t,q,μ) * dϑ₁dx₄(t,q)
dβ₂dx₁(t,q,μ) = d²Hdx₄dx₁(t,q,μ) * ϑ₂(t,q) + dHdx₄(t,q,μ) * dϑ₂dx₁(t,q)
dβ₂dx₂(t,q,μ) = d²Hdx₄dx₂(t,q,μ) * ϑ₂(t,q) + dHdx₄(t,q,μ) * dϑ₂dx₂(t,q)
dβ₂dx₃(t,q,μ) = d²Hdx₄dx₃(t,q,μ) * ϑ₂(t,q) + dHdx₄(t,q,μ) * dϑ₂dx₃(t,q)
dβ₂dx₄(t,q,μ) = d²Hdx₄dx₄(t,q,μ) * ϑ₂(t,q) + dHdx₄(t,q,μ) * dϑ₂dx₄(t,q)
dβ₃dx₁(t,q,μ) = d²Hdx₄dx₁(t,q,μ) * ϑ₃(t,q) + dHdx₄(t,q,μ) * dϑ₃dx₁(t,q)
dβ₃dx₂(t,q,μ) = d²Hdx₄dx₂(t,q,μ) * ϑ₃(t,q) + dHdx₄(t,q,μ) * dϑ₃dx₂(t,q)
dβ₃dx₃(t,q,μ) = d²Hdx₄dx₃(t,q,μ) * ϑ₃(t,q) + dHdx₄(t,q,μ) * dϑ₃dx₃(t,q)
dβ₃dx₄(t,q,μ) = d²Hdx₄dx₄(t,q,μ) * ϑ₃(t,q) + dHdx₄(t,q,μ) * dϑ₃dx₄(t,q)

dγ₁dx₁(t,q,μ) = dϑ₂dx₁(t,q) * dHdx₃(t,q,μ) + ϑ₂(t,q) * d²Hdx₃dx₁(t,q,μ) - d²Hdx₂dx₁(t,q,μ) * ϑ₃(t,q) - dHdx₂(t,q,μ) * dϑ₃dx₁(t,q)
dγ₁dx₂(t,q,μ) = dϑ₂dx₂(t,q) * dHdx₃(t,q,μ) + ϑ₂(t,q) * d²Hdx₃dx₂(t,q,μ) - d²Hdx₂dx₂(t,q,μ) * ϑ₃(t,q) - dHdx₂(t,q,μ) * dϑ₃dx₂(t,q)
dγ₁dx₃(t,q,μ) = dϑ₂dx₃(t,q) * dHdx₃(t,q,μ) + ϑ₂(t,q) * d²Hdx₃dx₃(t,q,μ) - d²Hdx₂dx₃(t,q,μ) * ϑ₃(t,q) - dHdx₂(t,q,μ) * dϑ₃dx₃(t,q)
dγ₁dx₄(t,q,μ) = dϑ₂dx₄(t,q) * dHdx₃(t,q,μ) + ϑ₂(t,q) * d²Hdx₃dx₄(t,q,μ) - d²Hdx₂dx₄(t,q,μ) * ϑ₃(t,q) - dHdx₂(t,q,μ) * dϑ₃dx₄(t,q)
dγ₂dx₁(t,q,μ) = dϑ₃dx₁(t,q) * dHdx₁(t,q,μ) + ϑ₃(t,q) * d²Hdx₁dx₁(t,q,μ) - d²Hdx₃dx₁(t,q,μ) * ϑ₁(t,q) - dHdx₃(t,q,μ) * dϑ₁dx₁(t,q)
dγ₂dx₂(t,q,μ) = dϑ₃dx₂(t,q) * dHdx₁(t,q,μ) + ϑ₃(t,q) * d²Hdx₁dx₂(t,q,μ) - d²Hdx₃dx₂(t,q,μ) * ϑ₁(t,q) - dHdx₃(t,q,μ) * dϑ₁dx₂(t,q)
dγ₂dx₃(t,q,μ) = dϑ₃dx₃(t,q) * dHdx₁(t,q,μ) + ϑ₃(t,q) * d²Hdx₁dx₃(t,q,μ) - d²Hdx₃dx₃(t,q,μ) * ϑ₁(t,q) - dHdx₃(t,q,μ) * dϑ₁dx₃(t,q)
dγ₂dx₄(t,q,μ) = dϑ₃dx₄(t,q) * dHdx₁(t,q,μ) + ϑ₃(t,q) * d²Hdx₁dx₄(t,q,μ) - d²Hdx₃dx₄(t,q,μ) * ϑ₁(t,q) - dHdx₃(t,q,μ) * dϑ₁dx₄(t,q)
dγ₃dx₁(t,q,μ) = dϑ₁dx₁(t,q) * dHdx₂(t,q,μ) + ϑ₁(t,q) * d²Hdx₂dx₁(t,q,μ) - d²Hdx₁dx₁(t,q,μ) * ϑ₂(t,q) - dHdx₁(t,q,μ) * dϑ₂dx₁(t,q)
dγ₃dx₂(t,q,μ) = dϑ₁dx₂(t,q) * dHdx₂(t,q,μ) + ϑ₁(t,q) * d²Hdx₂dx₂(t,q,μ) - d²Hdx₁dx₂(t,q,μ) * ϑ₂(t,q) - dHdx₁(t,q,μ) * dϑ₂dx₂(t,q)
dγ₃dx₃(t,q,μ) = dϑ₁dx₃(t,q) * dHdx₂(t,q,μ) + ϑ₁(t,q) * d²Hdx₂dx₃(t,q,μ) - d²Hdx₁dx₃(t,q,μ) * ϑ₂(t,q) - dHdx₁(t,q,μ) * dϑ₂dx₃(t,q)
dγ₃dx₄(t,q,μ) = dϑ₁dx₄(t,q) * dHdx₂(t,q,μ) + ϑ₁(t,q) * d²Hdx₂dx₄(t,q,μ) - d²Hdx₁dx₄(t,q,μ) * ϑ₂(t,q) - dHdx₁(t,q,μ) * dϑ₂dx₄(t,q)


# `orientation()` multiplies every subsystem, for the reason given on `ωabs`: the coordinate curl
# below carries `det(DF)`, and the rescaling factor has to be the phasespace Jacobian `√det Ω`,
# which is positive. Negating a divergence-free field leaves it divergence-free and a Hamiltonian
# vector field Hamiltonian, so nothing about the splitting changes — each `vᵢ` is still symplectic
# in the two variables it moves, and the six still sum to `v`.

function v₁(v, t, q, params)
    @unpack μ = params
    v[1] = + orientation() * dβ₃dx₂(t,q,μ)
    v[2] = - orientation() * dβ₃dx₁(t,q,μ)
    v[3] = 0
    v[4] = 0
end

function v₂(v, t, q, params)
    @unpack μ = params
    v[1] = - orientation() * dβ₂dx₃(t,q,μ)
    v[2] = 0
    v[3] = + orientation() * dβ₂dx₁(t,q,μ)
    v[4] = 0
end

function v₃(v, t, q, params)
    @unpack μ = params
    v[1] = 0
    v[2] = + orientation() * dβ₁dx₃(t,q,μ)
    v[3] = - orientation() * dβ₁dx₂(t,q,μ)
    v[4] = 0
end

function v₄(v, t, q, params)
    @unpack μ = params
    v[1] = + orientation() * dγ₁dx₄(t,q,μ)
    v[2] = 0
    v[3] = 0
    v[4] = - orientation() * dγ₁dx₁(t,q,μ)
end

function v₅(v, t, q, params)
    @unpack μ = params
    v[1] = 0
    v[2] = + orientation() * dγ₂dx₄(t,q,μ)
    v[3] = 0
    v[4] = - orientation() * dγ₂dx₂(t,q,μ)
end

function v₆(v, t, q, params)
    @unpack μ = params
    v[1] = 0
    v[2] = 0
    v[3] = + orientation() * dγ₃dx₄(t,q,μ)
    v[4] = - orientation() * dγ₃dx₃(t,q,μ)
end

@doc raw"""
    v(v, t, q, params)

The gyrokinetic guiding centre vector field in the rescaled time,

```math
\dfrac{dR}{ds} = \sigma \left( \dfrac{\partial \gamma}{\partial u} + \nabla \times \beta \right) ,
\qquad
\dfrac{du}{ds} = - \sigma \, \nabla \cdot \gamma ,
```

built from the potentials [`β`](@ref) and [`γ`](@ref), with ``\sigma =`` `orientation()` the sign of
`det(DF)` for this module's chart, as `ElectromagneticFields` generates it. The ``\nabla \times``
above is the *coordinate* curl, which carries `det(DF)`; ``\sigma`` restores the orientation, so the
bracket is the properly oriented curl and the field is the guiding centre one multiplied by the
phasespace Jacobian
``\omega_{abs} = J \, B^{\star}_{\parallel} > 0``. The independent variable is therefore related to
the physical time by ``dt = \omega_{abs} \, ds``, with ``s`` running *with* ``t`` in every chart.
See [`ωabs`](@ref) for why the sign is not optional.

``\sigma`` does not disturb the structure: negating a divergence-free field leaves it
divergence-free, so the splitting is volume preserving exactly as before.

The sub-fields `v₁` … `v₆` are the six subsystems this splits into; they sum to `v`.
"""
function v(v, t, q, params)
    @unpack μ = params
    v[1] = orientation() * (  dγ₁dx₄(t,q,μ) + dβ₃dx₂(t,q,μ) - dβ₂dx₃(t,q,μ) )
    v[2] = orientation() * (  dγ₂dx₄(t,q,μ) + dβ₁dx₃(t,q,μ) - dβ₃dx₁(t,q,μ) )
    v[3] = orientation() * (  dγ₃dx₄(t,q,μ) + dβ₂dx₁(t,q,μ) - dβ₁dx₂(t,q,μ) )
    v[4] = orientation() * ( -dγ₁dx₁(t,q,μ) - dγ₂dx₂(t,q,μ) - dγ₃dx₃(t,q,μ) )
    nothing
end
