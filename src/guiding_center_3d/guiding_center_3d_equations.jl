using LinearAlgebra
using Parameters

import GeometricEquations: HODEProblem
import GeometricSolutions: GeometricSolution, DataSeries, TimeSeries

export hamiltonian, hamiltonian_canonical
export hodeproblem, hodeproblem_canonical

# The constraints and the index accessors every generic expression below is written with. Included
# from here rather than from the equilibrium modules so that adding it did not need thirteen edits;
# `include` resolves relative to the file containing the call.
include("guiding_center_3d_constraints.jl")


@views ϑ₁(t, Q) = A₁(t, Q[1:3]) + Q[4] * b₁(t, Q[1:3])
@views ϑ₂(t, Q) = A₂(t, Q[1:3]) + Q[4] * b₂(t, Q[1:3])
@views ϑ₃(t, Q) = A₃(t, Q[1:3]) + Q[4] * b₃(t, Q[1:3])

# The named `v` and its first derivatives, kept because the second derivatives of the Hamiltonian in
# `guiding_center_3d_canonical.jl` are written out per index and read better that way. They are
# aliases of the generic accessors rather than second definitions of them.
for i in 1:3
    @eval $(Symbol("v", INDEX_SUBSCRIPTS[i]))(t, q, p) = vᵢ($(Val(i)), t, q, p)

    for j in 1:3
        @eval $(Symbol("dv", INDEX_SUBSCRIPTS[i], "dq", INDEX_SUBSCRIPTS[j]))(t, q, p) = dvᵢdqⱼ($(Val(i)), $(Val(j)), t, q, p)
        @eval $(Symbol("dv", INDEX_SUBSCRIPTS[i], "dp", INDEX_SUBSCRIPTS[j]))(t, q, p) = dvᵢdpⱼ($(Val(i)), $(Val(j)), t, q, p)
    end
end


u(t, q, p) = v₁(t, q, p) * g¹¹(t, q) * b₁(t, q) + v₂(t, q, p) * g²²(t, q) * b₂(t, q) + v₃(t, q, p) * g³³(t, q) * b₃(t, q)


function initial_momentum(tᵢ, Qᵢ::AbstractArray{T}) where {T<:Number}
    pᵢ = zeros(T, 3)
    pᵢ[1] = ϑ₁(tᵢ, Qᵢ)
    pᵢ[2] = ϑ₂(tᵢ, Qᵢ)
    pᵢ[3] = ϑ₃(tᵢ, Qᵢ)
    return pᵢ
end

function fix_initial_momentum(tᵢ, qᵢ::AbstractArray{T}, pᵢ::AbstractArray{T}) where {T<:Number}
    # `u` raises the index with the inverse metric; spelling the contraction out here without it
    # made this disagree with `u(t, q, p)` above in every curvilinear coordinate system.
    initial_momentum(tᵢ, [qᵢ..., u(tᵢ, qᵢ, pᵢ)])
end

function initial_conditions(tᵢ, Qᵢ::AbstractArray{T}) where {T<:Number}
    qᵢ = Qᵢ[1:3]
    pᵢ = initial_momentum(tᵢ, Qᵢ)
    (q=qᵢ, p=pᵢ)
end


# `rangemin`/`rangemax` take an evaluation point only for uniformity with the other generated field
# functions (`B(t,x)`, `A₁(t,x)`, …). The range of a coordinate is a property of the chart, not a
# field sampled at a point: `ElectromagneticFields` bakes `minx¹`…`maxx³` in as literals when it
# injects the field code, so the argument is discarded. Pass the origin rather than the `±Inf` that
# `xmin`/`xmax` are initialised to, which read as if the range were being asked for at the point at
# infinity.
function guiding_center_3d_periodicity(::Type{T}, periodic=true) where {T}
    xmin = -Inf * ones(T, 3)
    xmax = +Inf * ones(T, 3)

    if periodic
        xmin[1:3] .= rangemin(zeros(T, 3))
        xmax[1:3] .= rangemax(zeros(T, 3))
    end

    return (xmin, xmax)
end

guiding_center_3d_periodicity(::AbstractVector{<:AbstractArray{T}}, periodic=true) where {T<:Number} = guiding_center_3d_periodicity(T, periodic)
guiding_center_3d_periodicity(::AbstractArray{T}, periodic=true) where {T<:Number} = guiding_center_3d_periodicity(T, periodic)


hamiltonian(t, q, p, params) = g¹¹(t, q) * v₁(t, q, p)^2 / 2 + g²²(t, q) * v₂(t, q, p)^2 / 2 + g³³(t, q) * v₃(t, q, p)^2 / 2 + params.μ * B(t, q) + φ(t, q)
# The same Hamiltonian with the parallel velocity in place of the full kinetic energy. The two
# agree wherever the constraint `v × b = 0` holds, since `v` is then purely parallel and `|b| = 1`,
# so their difference along a trajectory is a measure of the drift off the constraint manifold —
# which is what `scripts/guiding_center_3d_*.jl` plot it for. It omitted `φ`, unlike `hamiltonian`
# beside it, which made the two incomparable for any equilibrium with a potential.
hamiltonian_u(t, q, p, params) = u(t, q, p)^2 / 2 + params.μ * B(t, q) + φ(t, q)


dHdq₁(t, q, p, params) = v₁(t, q, p) * g¹¹(t, q) * dv₁dq₁(t, q, p) +
                         v₂(t, q, p) * g²²(t, q) * dv₂dq₁(t, q, p) +
                         v₃(t, q, p) * g³³(t, q) * dv₃dq₁(t, q, p) +
                         v₁(t, q, p) * dg¹¹dx₁(t, q) * v₁(t, q, p) / 2 +
                         v₂(t, q, p) * dg²²dx₁(t, q) * v₂(t, q, p) / 2 +
                         v₃(t, q, p) * dg³³dx₁(t, q) * v₃(t, q, p) / 2 +
                         params.μ * dBdx₁(t, q) -
                         E₁(t, q)

dHdq₂(t, q, p, params) = v₁(t, q, p) * g¹¹(t, q) * dv₁dq₂(t, q, p) +
                         v₂(t, q, p) * g²²(t, q) * dv₂dq₂(t, q, p) +
                         v₃(t, q, p) * g³³(t, q) * dv₃dq₂(t, q, p) +
                         v₁(t, q, p) * dg¹¹dx₂(t, q) * v₁(t, q, p) / 2 +
                         v₂(t, q, p) * dg²²dx₂(t, q) * v₂(t, q, p) / 2 +
                         v₃(t, q, p) * dg³³dx₂(t, q) * v₃(t, q, p) / 2 +
                         params.μ * dBdx₂(t, q) -
                         E₂(t, q)

dHdq₃(t, q, p, params) = v₁(t, q, p) * g¹¹(t, q) * dv₁dq₃(t, q, p) +
                         v₂(t, q, p) * g²²(t, q) * dv₂dq₃(t, q, p) +
                         v₃(t, q, p) * g³³(t, q) * dv₃dq₃(t, q, p) +
                         v₁(t, q, p) * dg¹¹dx₃(t, q) * v₁(t, q, p) / 2 +
                         v₂(t, q, p) * dg²²dx₃(t, q) * v₂(t, q, p) / 2 +
                         v₃(t, q, p) * dg³³dx₃(t, q) * v₃(t, q, p) / 2 +
                         params.μ * dBdx₃(t, q) -
                         E₃(t, q)

dHdp₁(t, q, p, params) = v₁(t, q, p) * g¹¹(t, q) * dv₁dp₁(t, q, p) +
                         v₂(t, q, p) * g²²(t, q) * dv₂dp₁(t, q, p) +
                         v₃(t, q, p) * g³³(t, q) * dv₃dp₁(t, q, p)

dHdp₂(t, q, p, params) = v₁(t, q, p) * g¹¹(t, q) * dv₁dp₂(t, q, p) +
                         v₂(t, q, p) * g²²(t, q) * dv₂dp₂(t, q, p) +
                         v₃(t, q, p) * g³³(t, q) * dv₃dp₂(t, q, p)

dHdp₃(t, q, p, params) = v₁(t, q, p) * g¹¹(t, q) * dv₁dp₃(t, q, p) +
                         v₂(t, q, p) * g²²(t, q) * dv₂dp₃(t, q, p) +
                         v₃(t, q, p) * g³³(t, q) * dv₃dp₃(t, q, p)

for i in 1:3
    @eval @inline dHdqᵢ(::Val{$i}, t, q, p, params) = $(Symbol("dHdq", INDEX_SUBSCRIPTS[i]))(t, q, p, params)
    @eval @inline dHdpᵢ(::Val{$i}, t, q, p, params) = $(Symbol("dHdp", INDEX_SUBSCRIPTS[i]))(t, q, p, params)
end


# The two Poisson brackets the Lagrange multipliers are built from,
# `{f, l} = Σᵢ (∂f/∂qᵢ ∂l/∂pᵢ - ∂f/∂pᵢ ∂l/∂qᵢ)`.
@inline bracket_gg(k₁::Val, k₂::Val, t, q, p) = contract(l ->
    dgᵏdqₗ(k₁, l, t, q, p) * dgᵏdpₗ(k₂, l, t, q, p) -
    dgᵏdpₗ(k₁, l, t, q, p) * dgᵏdqₗ(k₂, l, t, q, p))

@inline bracket_gH(k::Val, t, q, p, params) = contract(l ->
    dgᵏdqₗ(k, l, t, q, p) * dHdpᵢ(l, t, q, p, params) -
    dgᵏdpₗ(k, l, t, q, p) * dHdqᵢ(l, t, q, p, params))


"""
    λₒ(t, q, p, c = default_constraint_pair())

The Poisson bracket `{g₁, g₂}` of the two constraints of the pair `c`, which divides both Lagrange
multipliers. It equals `bₘ [B + (p-A)·(∇×b)]` with `m` the index of the constraint the pair omits, so
it is where the formulation becomes singular — see
`scripts/study_guiding_center_3d_conditioning.jl`.
"""
λₒ(t, q, p, c = default_constraint_pair()) = bracket_gg(c[1], c[2], t, q, p)

"""
    λ₁(t, q, p, params, c = default_constraint_pair())
    λ₂(t, q, p, params, c = default_constraint_pair())

The two Lagrange multipliers, `λ₁ = {g₂, H} / {g₁, g₂}` and `λ₂ = -{g₁, H} / {g₁, g₂}`, for the
constraint pair `c`.
"""
λ₁(t, q, p, params, c = default_constraint_pair()) = +bracket_gH(c[2], t, q, p, params) / λₒ(t, q, p, c)
λ₂(t, q, p, params, c = default_constraint_pair()) = -bracket_gH(c[1], t, q, p, params) / λₒ(t, q, p, c)

hamiltonian_canonical(t, q, p, params, c = default_constraint_pair()) =
    hamiltonian(t, q, p, params) + λ₁(t, q, p, params, c) * gᵏ(c[1], t, q, p) +
                                   λ₂(t, q, p, params, c) * gᵏ(c[2], t, q, p)


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


function guiding_center_3d_v(v, t, q, p, params, c = default_constraint_pair())
    l₁ = λ₁(t, q, p, params, c)
    l₂ = λ₂(t, q, p, params, c)

    components!(v, i -> dHdpᵢ(i, t, q, p, params) + l₁ * dgᵏdpₗ(c[1], i, t, q, p)
                                                  + l₂ * dgᵏdpₗ(c[2], i, t, q, p))
    nothing
end

function guiding_center_3d_f(f, t, q, p, params, c = default_constraint_pair())
    l₁ = λ₁(t, q, p, params, c)
    l₂ = λ₂(t, q, p, params, c)

    components!(f, i -> -dHdqᵢ(i, t, q, p, params) - l₁ * dgᵏdqₗ(c[1], i, t, q, p)
                                                   - l₂ * dgᵏdqₗ(c[2], i, t, q, p))
    nothing
end


"""
    hodeproblem(q₀, p₀; kwargs...)
    hodeproblem(x₀ = qᵢ; kwargs...)
    hodeproblem(ics::NamedTuple; kwargs...)

The constrained canonical guiding centre system as an `HODEProblem` in the position and its
conjugate momentum — the Hamilton-Dirac form, with the Lagrange multipliers substituted.

The first form takes the position and momentum directly. The second takes the four-component state
``(x, u)`` and recovers the momentum from it through `initial_conditions`; this is the form
the module constant `qᵢ` is in. The third takes the named tuple that every `initial_conditions_*`
returns, so `hodeproblem(initial_conditions_barely_passing())` carries that condition's own `μ`.

`constraints` selects which pair of the three constraints is retained; see
[`constraint_pair`](@ref). It defaults to [`default_constraints`](@ref), which each equilibrium sets
to a pair that is regular at its own initial condition. Where the pair is singular the multipliers
are infinite and the problem cannot be integrated at all, so this is not a free choice.
"""
function hodeproblem(q₀::AbstractVector, p₀::AbstractVector; timespan = DEFAULT_TIMESPAN,
                     timestep = DEFAULT_TIMESTEP, parameters = default_parameters(), periodic = true,
                     constraints = default_constraints())
    c = constraint_pair(constraints)

    _v(v, t, q, p, params) = guiding_center_3d_v(v, t, q, p, params, c)
    _f(f, t, q, p, params) = guiding_center_3d_f(f, t, q, p, params, c)

    HODEProblem(
        _v,
        _f,
        hamiltonian,
        timespan, timestep, q₀, p₀;
        parameters=parameters,
        periodicity=guiding_center_3d_periodicity(q₀, periodic))
end

function hodeproblem(x₀::AbstractVector = qᵢ; timespan = DEFAULT_TIMESPAN, kwargs...)
    ics = initial_conditions(timespan[begin], x₀)
    hodeproblem(ics.q, ics.p; timespan = timespan, kwargs...)
end

hodeproblem(ics::NamedTuple; kwargs...) = hodeproblem(ics.q, ics.p; parameters = ics.params, kwargs...)
