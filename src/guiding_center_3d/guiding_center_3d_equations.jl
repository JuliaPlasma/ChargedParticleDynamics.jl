using LinearAlgebra
using Parameters

import GeometricEquations: HODEProblem
using ..FieldPoints
using ..ChargedParticleDynamics: periodic_domain

export hamiltonian, hamiltonian_canonical

# The constraints and the index accessors every generic expression below is written with.
# `include` resolves relative to the file containing the call. The compact form is pulled in the
# same way at the foot of this file.
include("guiding_center_3d_constraints.jl")

# The one-form on the four-component state `Q = (x, u)`, read from a `FieldPoint` of `Q`; the field
# is evaluated at `Q[1:3]`.
ϑ₁(t, Q) = A₁(t, Q) + Q[4] * b₁(t, Q)
ϑ₂(t, Q) = A₂(t, Q) + Q[4] * b₂(t, Q)
ϑ₃(t, Q) = A₃(t, Q) + Q[4] * b₃(t, Q)

# The named `v` and its first derivatives, kept because the second derivatives of the Hamiltonian in
# `guiding_center_3d_canonical.jl` are written out per index and read better that way. They are
# aliases of the generic accessors rather than second definitions of them.
for i in 1:3
    @eval $(Symbol("v", INDEX_SUBSCRIPTS[i]))(t, q, p) = vᵢ($(Val(i)), t, q, p)

    for j in 1:3
        @eval $(Symbol("dv", INDEX_SUBSCRIPTS[i], "dq", INDEX_SUBSCRIPTS[j]))(t, q, p) = dvᵢdqⱼ(
            $(Val(i)), $(Val(j)), t, q, p)
        @eval $(Symbol("dv", INDEX_SUBSCRIPTS[i], "dp", INDEX_SUBSCRIPTS[j]))(t, q, p) = dvᵢdpⱼ(
            $(Val(i)), $(Val(j)), t, q, p)
    end
end

function u(t, q::FieldPoint, p)
    v₁(t, q, p) * g¹¹(t, q) * b₁(t, q) + v₂(t, q, p) * g²²(t, q) * b₂(t, q) +
    v₃(t, q, p) * g³³(t, q) * b₃(t, q)
end

"""
    u(t, q, p, params)

The parallel velocity `Σᵢ gⁱⁱ (pᵢ - Aᵢ) bᵢ` of the state `(q, p)`.
"""
u(t, q, p, params) = u(t, fieldpoint(params.field, t, q), p)

function initial_momentum(tᵢ, Qᵢ::AbstractArray{T}, params) where {T <: Number}
    P = fieldpoint(params.field, tᵢ, Qᵢ)
    pᵢ = zeros(T, 3)
    pᵢ[1] = ϑ₁(tᵢ, P)
    pᵢ[2] = ϑ₂(tᵢ, P)
    pᵢ[3] = ϑ₃(tᵢ, P)
    return pᵢ
end

function fix_initial_momentum(tᵢ, qᵢ::AbstractArray{T}, pᵢ::AbstractArray{T}, params) where {T <:
                                                                                             Number}
    # `u` raises the index with the inverse metric; spelling the contraction out here without it
    # made this disagree with `u(t, q, p)` above in every curvilinear coordinate system.
    initial_momentum(tᵢ, [qᵢ..., u(tᵢ, qᵢ, pᵢ, params)], params)
end

function initial_conditions(tᵢ, Qᵢ::AbstractArray{T}, params) where {T <: Number}
    qᵢ = Qᵢ[1:3]
    pᵢ = initial_momentum(tᵢ, Qᵢ, params)
    (q = qᵢ, p = pᵢ)
end

# Which coordinates are periodic, and on what range, is a property of the chart, answered by the
# field; see `periodic_domain`.
function guiding_center_3d_periodicity(::Type{T}, field, periodic = true) where {T}
    periodic ? periodic_domain(field, T, 3) : (fill(-T(Inf), 3), fill(+T(Inf), 3))
end

function guiding_center_3d_periodicity(
        ::AbstractVector{<:AbstractArray{T}}, field, periodic = true) where {T <: Number}
    guiding_center_3d_periodicity(T, field, periodic)
end
function guiding_center_3d_periodicity(
        ::AbstractArray{T}, field, periodic = true) where {T <: Number}
    guiding_center_3d_periodicity(T, field, periodic)
end

function hamiltonian(t, q, p, params)
    q = fieldpoint(params.field, t, q)
    g¹¹(t, q) * v₁(t, q, p)^2 / 2 + g²²(t, q) * v₂(t, q, p)^2 / 2 +
    g³³(t, q) * v₃(t, q, p)^2 / 2 + params.μ * B(t, q) + φ(t, q)
end
# The same Hamiltonian with the parallel velocity in place of the full kinetic energy. The two
# agree wherever the constraint `v × b = 0` holds, since `v` is then purely parallel and `|b| = 1`,
# so their difference along a trajectory is a measure of the drift off the constraint manifold —
# which is what `scripts/guiding_center_3d_*.jl` plot it for. It omitted `φ`, unlike `hamiltonian`
# beside it, which made the two incomparable for any equilibrium with a potential.
function hamiltonian_u(t, q, p, params)
    q = fieldpoint(params.field, t, q)
    u(t, q, p)^2 / 2 + params.μ * B(t, q) + φ(t, q)
end

function dHdq₁(t, q, p, params)
    v₁(t, q, p) * g¹¹(t, q) * dv₁dq₁(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * dv₂dq₁(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * dv₃dq₁(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₁(t, q) * v₁(t, q, p) / 2 +
    v₂(t, q, p) * dg²²dx₁(t, q) * v₂(t, q, p) / 2 +
    v₃(t, q, p) * dg³³dx₁(t, q) * v₃(t, q, p) / 2 +
    params.μ * dBdx₁(t, q) -
    E₁(t, q)
end

function dHdq₂(t, q, p, params)
    v₁(t, q, p) * g¹¹(t, q) * dv₁dq₂(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * dv₂dq₂(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * dv₃dq₂(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₂(t, q) * v₁(t, q, p) / 2 +
    v₂(t, q, p) * dg²²dx₂(t, q) * v₂(t, q, p) / 2 +
    v₃(t, q, p) * dg³³dx₂(t, q) * v₃(t, q, p) / 2 +
    params.μ * dBdx₂(t, q) -
    E₂(t, q)
end

function dHdq₃(t, q, p, params)
    v₁(t, q, p) * g¹¹(t, q) * dv₁dq₃(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * dv₂dq₃(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * dv₃dq₃(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₃(t, q) * v₁(t, q, p) / 2 +
    v₂(t, q, p) * dg²²dx₃(t, q) * v₂(t, q, p) / 2 +
    v₃(t, q, p) * dg³³dx₃(t, q) * v₃(t, q, p) / 2 +
    params.μ * dBdx₃(t, q) -
    E₃(t, q)
end

function dHdp₁(t, q, p, params)
    v₁(t, q, p) * g¹¹(t, q) * dv₁dp₁(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * dv₂dp₁(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * dv₃dp₁(t, q, p)
end

function dHdp₂(t, q, p, params)
    v₁(t, q, p) * g¹¹(t, q) * dv₁dp₂(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * dv₂dp₂(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * dv₃dp₂(t, q, p)
end

function dHdp₃(t, q, p, params)
    v₁(t, q, p) * g¹¹(t, q) * dv₁dp₃(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * dv₂dp₃(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * dv₃dp₃(t, q, p)
end

for i in 1:3
    @eval @inline dHdqᵢ(::Val{$i}, t, q, p, params) = $(Symbol("dHdq", INDEX_SUBSCRIPTS[i]))(t, q, p, params)
    @eval @inline dHdpᵢ(::Val{$i}, t, q, p, params) = $(Symbol("dHdp", INDEX_SUBSCRIPTS[i]))(t, q, p, params)
end

# The two Poisson brackets the Lagrange multipliers are built from,
# `{f, l} = Σᵢ (∂f/∂qᵢ ∂l/∂pᵢ - ∂f/∂pᵢ ∂l/∂qᵢ)`.
@inline bracket_gg(k₁::Val, k₂::Val, t, q, p) = contract(l -> dgᵏdqₗ(k₁, l, t, q, p) *
                                                              dgᵏdpₗ(k₂, l, t, q, p) -
                                                              dgᵏdpₗ(k₁, l, t, q, p) *
                                                              dgᵏdqₗ(k₂, l, t, q, p))

@inline bracket_gH(k::Val, t, q, p, params) = contract(l -> dgᵏdqₗ(k, l, t, q, p) *
                                                            dHdpᵢ(l, t, q, p, params) -
                                                            dgᵏdpₗ(k, l, t, q, p) *
                                                            dHdqᵢ(l, t, q, p, params))

"""
    λₒ(t, q, p, params, c)

The Poisson bracket `{g₁, g₂}` of the two constraints of the pair `c`, which divides both Lagrange
multipliers. It equals `±bₘ [B + (p-A)·(∇×b)]` with `m` the index of the constraint the pair omits, so
it is where the formulation becomes singular — see
`scripts/study_guiding_center_3d_conditioning.jl`.

The sign depends on the pair's ordering and is `+` for `:g31` and `:g12`, `-` for `:g23`; see
[`constraint_pair`](@ref). Only `λₒ = 0` matters for whether the pair is usable, so nothing turns on
it, but `λₒ` and `bₘ` do not always share a sign and the tabulated values reflect that.

The four-argument form `λₒ(t, q, p, c)` takes a `FieldPoint` in place of `q`.
"""
λₒ(t, q::FieldPoint, p, c) = bracket_gg(c[1], c[2], t, q, p)
λₒ(t, q, p, params, c) = λₒ(t, fieldpoint(params.field, t, q), p, c)

"""
    λ₁(t, q, p, params, c)
    λ₂(t, q, p, params, c)

The two Lagrange multipliers, `λ₁ = {g₂, H} / {g₁, g₂}` and `λ₂ = -{g₁, H} / {g₁, g₂}`, for the
constraint pair `c`.
"""
function λ₁(t, q, p, params, c)
    q = fieldpoint(params.field, t, q)
    +bracket_gH(c[2], t, q, p, params) / λₒ(t, q, p, c)
end
function λ₂(t, q, p, params, c)
    q = fieldpoint(params.field, t, q)
    -bracket_gH(c[1], t, q, p, params) / λₒ(t, q, p, c)
end

"""
    multipliers(t, q, p, params, c)

Both Lagrange multipliers as a tuple, sharing the one evaluation of `λₒ` that divides them.

`λ₁` and `λ₂` each recompute it, and the right-hand sides need both, so calling them separately
evaluated the twelve Poisson-bracket terms behind `{g₁, g₂}` four times per step-stage rather than
twice. Nothing else divides by `λₒ`, so this is the only place the sharing is worth spelling out.
"""
@inline function multipliers(t, q, p, params, c)
    q = fieldpoint(params.field, t, q)
    lo = λₒ(t, q, p, c)
    (+bracket_gH(c[2], t, q, p, params) / lo,
        -bracket_gH(c[1], t, q, p, params) / lo)
end

function hamiltonian_canonical(t, q, p, params, c)
    q = fieldpoint(params.field, t, q)
    hamiltonian(t, q, p, params) + λ₁(t, q, p, params, c) * gᵏ(c[1], t, q, p) +
    λ₂(t, q, p, params, c) * gᵏ(c[2], t, q, p)
end

# Fill or accumulate into the three components of a right-hand side from a function of the
# compile-time index.
@inline function components!(x, f)
    x[1] = f(Val(1))
    x[2] = f(Val(2))
    x[3] = f(Val(3))
    x
end

@inline function addcomponents!(x, f)
    x[1] += f(Val(1))
    x[2] += f(Val(2))
    x[3] += f(Val(3))
    x
end

# `fieldpoint` evaluates each field tensor once and the expressions below read the results out of
# it; see the header of `guiding_center_3d_constraints.jl` for why that is worth doing. `F` goes
# where `q` would: everything downstream is generic in that argument.
function guiding_center_3d_v(v, t, q, p, params, c)
    F = fieldpoint(params.field, t, q)
    l₁, l₂ = multipliers(t, F, p, params, c)

    components!(v,
        i -> dHdpᵢ(i, t, F, p, params) + l₁ * dgᵏdpₗ(c[1], i, t, F, p)
             + l₂ * dgᵏdpₗ(c[2], i, t, F, p))
    nothing
end

function guiding_center_3d_f(f, t, q, p, params, c)
    F = fieldpoint(params.field, t, q)
    l₁, l₂ = multipliers(t, F, p, params, c)

    components!(f,
        i -> -dHdqᵢ(i, t, F, p, params) - l₁ * dgᵏdqₗ(c[1], i, t, F, p)
             -
             l₂ * dgᵏdqₗ(c[2], i, t, F, p))
    nothing
end

"""
    hodeproblem(q₀, p₀; timespan, timestep, parameters, constraints, periodic = true)
    hodeproblem(x₀; timespan, timestep, parameters, constraints, periodic = true)
    hodeproblem(ics::NamedTuple; kwargs...)

The constrained canonical guiding centre system as an `HODEProblem` in the position and its
conjugate momentum — the Hamilton-Dirac form, with the Lagrange multipliers substituted.

The first form takes the position and momentum directly. The second takes the four-component state
``(x, u)`` and recovers the momentum from it through `initial_conditions`; this is the form
an equilibrium module's constant `qᵢ` is in. The third takes the named tuple that every
`initial_conditions_*` returns, so `hodeproblem(initial_conditions_barely_passing())` carries that
condition's own `μ`. Each equilibrium module defaults all of these to its own.

`constraints` selects which pair of the three constraints is retained; see
[`constraint_pair`](@ref). Each equilibrium module defaults it to its `default_constraints()`, a
pair that is regular at its own initial condition. Where the pair is singular the multipliers
are infinite and the problem cannot be integrated at all, so this is not a free choice.
"""
function hodeproblem(
        q₀::AbstractVector, p₀::AbstractVector; timespan, timestep, parameters,
        periodic = true, constraints)
    c = constraint_pair(constraints)

    _v(v, t, q, p, params) = guiding_center_3d_v(v, t, q, p, params, c)
    _f(f, t, q, p, params) = guiding_center_3d_f(f, t, q, p, params, c)

    HODEProblem(
        _v,
        _f,
        hamiltonian,
        timespan, timestep, q₀, p₀;
        parameters = parameters,
        periodicity = guiding_center_3d_periodicity(q₀, parameters.field, periodic))
end

function hodeproblem(x₀::AbstractVector; timespan, parameters, kwargs...)
    ics = initial_conditions(timespan[begin], x₀, parameters)
    hodeproblem(ics.q, ics.p; timespan = timespan, parameters = parameters, kwargs...)
end

function hodeproblem(ics::NamedTuple; kwargs...)
    hodeproblem(ics.q, ics.p; parameters = ics.params, kwargs...)
end

# The compact form, Eq. (29), which drops the constraint-proportional terms from the right-hand side.
# It needs the brackets, the Hamiltonian gradients and `u` from this file, and nothing at all from
# `guiding_center_3d_canonical.jl` — no second derivatives — so it is included here rather than there.
include("guiding_center_3d_compact.jl")
