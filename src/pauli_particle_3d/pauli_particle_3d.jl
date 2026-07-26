
using Parameters

using GeometricEquations: HODEProblem, IODEProblem, LODEProblem, PODEProblem


ϑ₁(t, q, v) = g₁₁(t, q) * v[1] + A₁(t, q)
ϑ₂(t, q, v) = g₂₂(t, q) * v[2] + A₂(t, q)
ϑ₃(t, q, v) = g₃₃(t, q) * v[3] + A₃(t, q)

function ϑ(θ, t, q, v)
    θ[1] = ϑ₁(t, q, v)
    θ[2] = ϑ₂(t, q, v)
    θ[3] = ϑ₃(t, q, v)
    nothing
end

v¹(t, q, p) = g¹¹(t, q) * (p[1] - A₁(t, q))
v²(t, q, p) = g²²(t, q) * (p[2] - A₂(t, q))
v³(t, q, p) = g³³(t, q) * (p[3] - A₃(t, q))

v(t, q, p) = [v¹(t, q, p), v²(t, q, p), v³(t, q, p)]


# ϕ₀(x) = E₀*sin(2π*x)
# ϕ(t,q) = ϕ₀(q[3])

# E₁(t,q) = zero(eltype(q))
# E₂(t,q) = zero(eltype(q))
# E₃(t,q) = - 2π*E₀*cos(2π*q[3])

function hamiltonian(t, q, p, params)
    @unpack μ = params
    0.5 * (g₁₁(t, q) * v¹(t, q, p)^2 + g₂₂(t, q) * v²(t, q, p)^2 + g₃₃(t, q) * v³(t, q, p)^2) + μ * B(t, q) + φ(t, q)
end


function initial_conditions(x₀, v₀)
    u₀ = v₀' * b(0, x₀)
    vpar = u₀ .* b⃗(0, x₀)
    vper = v₀ .- vpar
    μ = vper' * vper / 2 / B(0, x₀)

    (q = x₀, v = vpar, params = (μ = μ,))
end


function pauli_particle_3d_pᵢ(tᵢ, qᵢ, vᵢ)
    pᵢ = zero(qᵢ)
    ϑ(pᵢ, tᵢ, qᵢ, vᵢ)
    return pᵢ
end


function pauli_particle_3d_pode_v(v, t, q, p, params)
    v[1] = v¹(t, q, p)
    v[2] = v²(t, q, p)
    v[3] = v³(t, q, p)
    nothing
end

function pauli_particle_3d_pode_f(f, t, q, p, params)
    @unpack μ = params
    f[1] = dA₁dx₁(t, q) * v¹(t, q, p) + dA₂dx₁(t, q) * v²(t, q, p) + dA₃dx₁(t, q) * v³(t, q, p) + E₁(t, q) - μ * dBdx₁(t, q) +
           (dg₁₁dx₁(t, q) * v¹(t, q, p)^2 + dg₂₂dx₁(t, q) * v²(t, q, p)^2 + dg₃₃dx₁(t, q) * v³(t, q, p)^2) / 2
    f[2] = dA₁dx₂(t, q) * v¹(t, q, p) + dA₂dx₂(t, q) * v²(t, q, p) + dA₃dx₂(t, q) * v³(t, q, p) + E₂(t, q) - μ * dBdx₂(t, q) +
           (dg₁₁dx₂(t, q) * v¹(t, q, p)^2 + dg₂₂dx₂(t, q) * v²(t, q, p)^2 + dg₃₃dx₂(t, q) * v³(t, q, p)^2) / 2
    f[3] = dA₁dx₃(t, q) * v¹(t, q, p) + dA₂dx₃(t, q) * v²(t, q, p) + dA₃dx₃(t, q) * v³(t, q, p) + E₃(t, q) - μ * dBdx₃(t, q) +
           (dg₁₁dx₃(t, q) * v¹(t, q, p)^2 + dg₂₂dx₃(t, q) * v²(t, q, p)^2 + dg₃₃dx₃(t, q) * v³(t, q, p)^2) / 2
    nothing
end


pauli_particle_3d_iode_ϑ(θ, t, q, v, params) = ϑ(θ, t, q, v)

function pauli_particle_3d_iode_f(f, t, q, v, params)
    @unpack μ = params
    f[1] = dA₁dx₁(t, q) * v[1] + dA₂dx₁(t, q) * v[2] + dA₃dx₁(t, q) * v[3] + E₁(t, q) - μ * dBdx₁(t, q) +
           (dg₁₁dx₁(t, q) * v[1]^2 + dg₂₂dx₁(t, q) * v[2]^2 + dg₃₃dx₁(t, q) * v[3]^2) / 2
    f[2] = dA₁dx₂(t, q) * v[1] + dA₂dx₂(t, q) * v[2] + dA₃dx₂(t, q) * v[3] + E₂(t, q) - μ * dBdx₂(t, q) +
           (dg₁₁dx₂(t, q) * v[1]^2 + dg₂₂dx₂(t, q) * v[2]^2 + dg₃₃dx₂(t, q) * v[3]^2) / 2
    f[3] = dA₁dx₃(t, q) * v[1] + dA₂dx₃(t, q) * v[2] + dA₃dx₃(t, q) * v[3] + E₃(t, q) - μ * dBdx₃(t, q) +
           (dg₁₁dx₃(t, q) * v[1]^2 + dg₂₂dx₃(t, q) * v[2]^2 + dg₃₃dx₃(t, q) * v[3]^2) / 2
    nothing
end

# The projection force is gᵢ = (∂ϑⱼ/∂qⁱ) λʲ = (∂Aⱼ/∂qⁱ + ∂gⱼⱼ/∂qⁱ vʲ) λʲ. It was previously built
# from `v` alone, leaving the multiplier `λ` unused, so the projection did not depend on what it
# was projecting.
function pauli_particle_3d_iode_g(g, t, q, v, λ, params)
    g[1] = dA₁dx₁(t, q) * λ[1] + dA₂dx₁(t, q) * λ[2] + dA₃dx₁(t, q) * λ[3] +
           dg₁₁dx₁(t, q) * v[1] * λ[1] + dg₂₂dx₁(t, q) * v[2] * λ[2] + dg₃₃dx₁(t, q) * v[3] * λ[3]
    g[2] = dA₁dx₂(t, q) * λ[1] + dA₂dx₂(t, q) * λ[2] + dA₃dx₂(t, q) * λ[3] +
           dg₁₁dx₂(t, q) * v[1] * λ[1] + dg₂₂dx₂(t, q) * v[2] * λ[2] + dg₃₃dx₂(t, q) * v[3] * λ[3]
    g[3] = dA₁dx₃(t, q) * λ[1] + dA₂dx₃(t, q) * λ[2] + dA₃dx₃(t, q) * λ[3] +
           dg₁₁dx₃(t, q) * v[1] * λ[1] + dg₂₂dx₃(t, q) * v[2] * λ[2] + dg₃₃dx₃(t, q) * v[3] * λ[3]
    nothing
end


"""
    podeproblem(q₀, v₀; kwargs...)
    podeproblem(ics::NamedTuple; kwargs...)
    podeproblem(; kwargs...)

The Pauli particle as a `PODEProblem`.

!!! note "`v₀` is the parallel velocity"
    The state is the position and the *parallel* velocity, not the full one, and `μ` is the moment
    of the perpendicular part. `initial_conditions(x₀, v₀)` performs that split and returns both,
    which is why the no-argument form goes through it rather than handing the module's `vᵢ` to the
    constructor directly.

The second form takes the named tuple `initial_conditions` returns; the third splits the module's
own `(qᵢ, vᵢ)`.
"""
function podeproblem(q₀::AbstractVector, v₀::AbstractVector; timespan = DEFAULT_TIMESPAN,
                     timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
    PODEProblem(
        pauli_particle_3d_pode_v,
        pauli_particle_3d_pode_f,
        timespan, timestep, q₀, pauli_particle_3d_pᵢ(timespan[begin], q₀, v₀);
        parameters=parameters,
        invariants=(h=hamiltonian,)
    )
end


"""
    hodeproblem(q₀, v₀; kwargs...)
    hodeproblem(ics::NamedTuple; kwargs...)

The Pauli particle as an `HODEProblem`; see [`podeproblem`](@ref) for the argument forms.
"""
function hodeproblem(q₀::AbstractVector, v₀::AbstractVector; timespan = DEFAULT_TIMESPAN,
                     timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
    HODEProblem(
        pauli_particle_3d_pode_v,
        pauli_particle_3d_pode_f,
        hamiltonian,
        timespan, timestep, q₀, pauli_particle_3d_pᵢ(timespan[begin], q₀, v₀);
        parameters=parameters)
end


"""
    iodeproblem(q₀, v₀; kwargs...)
    iodeproblem(ics::NamedTuple; kwargs...)

The Pauli particle as an `IODEProblem`; see [`podeproblem`](@ref) for the argument forms.
"""
function iodeproblem(q₀::AbstractVector, v₀::AbstractVector; timespan = DEFAULT_TIMESPAN,
                     timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
    IODEProblem(
        pauli_particle_3d_iode_ϑ,
        pauli_particle_3d_iode_f,
        pauli_particle_3d_iode_g,
        timespan, timestep, q₀, pauli_particle_3d_pᵢ(timespan[begin], q₀, v₀);
        parameters=parameters,
        invariants=(h=hamiltonian,),
        v̄=pauli_particle_3d_pode_v)
end


for problem in (:podeproblem, :hodeproblem, :iodeproblem)
    # `μ` on its own, for callers that have the moment rather than a parameter tuple
    @eval $problem(q₀::AbstractVector, v₀::AbstractVector, μ::Real; kwargs...) =
        $problem(q₀, v₀; parameters = (μ = μ,), kwargs...)

    # the named tuple `initial_conditions` returns, whose `v` is already the parallel velocity
    @eval $problem(ics::NamedTuple; kwargs...) =
        $problem(ics.q, ics.v; parameters = ics.params, kwargs...)

    # no arguments: split the module's own `(qᵢ, vᵢ)`. This must go through `initial_conditions`
    # rather than defaulting `v₀ = vᵢ` in the signature above — `vᵢ` is the *full* velocity, and
    # handing it to the constructor as if it were the parallel one puts the particle on a
    # trajectory the solver cannot follow.
    @eval $problem(; kwargs...) = $problem(initial_conditions(qᵢ, vᵢ); kwargs...)
end


# function lodeproblem(q₀::AbstractVector, v₀::AbstractVector, parameters::NamedTuple; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP)
#     LODEProblem(
#         pauli_particle_3d_iode_ϑ,
#         pauli_particle_3d_iode_f,
#         pauli_particle_3d_iode_g,
#         pauli_particle_3d_ω, lagrangian,
#         timespan, step, q₀, pauli_particle_3d_pᵢ(q₀, v₀);
#         parameters = parameters,
#         invariants = (h = hamiltonian,))
# end

# lodeproblem(q₀::AbstractVector, v₀::AbstractVector, μ::Real) = lodeproblem(q₀, v₀, (μ=μ,))

# lodeproblem(qᵢ=qᵢ, vᵢ=vᵢ) = lodeproblem(initial_conditions(qᵢ, vᵢ)...)
