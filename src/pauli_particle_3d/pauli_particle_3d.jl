
using Parameters

using GeometricEquations: HODEProblem, IODEProblem, LODEProblem, PODEProblem


ϑ₁(t, q, v) = g₁₁(t,q) * v[1] + A₁(t, q)
ϑ₂(t, q, v) = g₂₂(t,q) * v[2] + A₂(t, q)
ϑ₃(t, q, v) = g₃₃(t,q) * v[3] + A₃(t, q)

function ϑ(θ, t, q, v)
    θ[1] = ϑ₁(t,q,v)
    θ[2] = ϑ₂(t,q,v)
    θ[3] = ϑ₃(t,q,v)
    nothing
end

v¹(t, q, p) = g¹¹(t, q) * (p[1] - A₁(t, q))
v²(t, q, p) = g²²(t, q) * (p[2] - A₂(t, q))
v³(t, q, p) = g³³(t, q) * (p[3] - A₃(t, q))

# ϕ₀(x) = E₀*sin(2π*x)
# ϕ(t,q) = ϕ₀(q[3])

# E₁(t,q) = zero(eltype(q))
# E₂(t,q) = zero(eltype(q))
# E₃(t,q) = - 2π*E₀*cos(2π*q[3])

function hamiltonian(t, q, p, params)
    @unpack μ = params
    0.5 * (g₁₁(t,q) * v¹(t,q,p)^2 + g₂₂(t,q) * v²(t,q,p)^2 + g₃₃(t,q) * v³(t,q,p)^2) + μ*B(t,q) + φ(t,q)
end


function initial_conditions(x₀, v₀)
    u₀ = v₀' * b(0, x₀)
    vpar = u₀ .* b⃗(0, x₀)
    vper = v₀ .- vpar
    μ = vper' * vper / 2 / B(0, x₀)
    
    (x₀, vpar, (μ = μ,))
end


function pauli_particle_3d_pᵢ(qᵢ, vᵢ, tᵢ=0)
    pᵢ = zero(qᵢ)
    ϑ(pᵢ, tᵢ, qᵢ, vᵢ)
    return pᵢ
end


function pauli_particle_3d_pode_v(v, t, q, p, params)
    v[1] = v¹(t,q,p)
    v[2] = v²(t,q,p)
    v[3] = v³(t,q,p)
    nothing
end

function pauli_particle_3d_pode_f(f, t, q, p, params)
    @unpack μ = params
    f[1] = dA₁dx₁(t,q) * v¹(t,q,p) + dA₂dx₁(t,q) * v²(t,q,p) + dA₃dx₁(t,q) * v³(t,q,p) + E₁(t,q) - μ * dBdx₁(t,q) +
           (dg₁₁dx₁(t,q) * v¹(t,q,p)^2 + dg₂₂dx₁(t,q) * v²(t,q,p)^2 + dg₃₃dx₁(t,q) * v³(t,q,p)^2) / 2
    f[2] = dA₁dx₂(t,q) * v¹(t,q,p) + dA₂dx₂(t,q) * v²(t,q,p) + dA₃dx₂(t,q) * v³(t,q,p) + E₂(t,q) - μ * dBdx₂(t,q) +
           (dg₁₁dx₂(t,q) * v¹(t,q,p)^2 + dg₂₂dx₂(t,q) * v²(t,q,p)^2 + dg₃₃dx₂(t,q) * v³(t,q,p)^2) / 2
    f[3] = dA₁dx₃(t,q) * v¹(t,q,p) + dA₂dx₃(t,q) * v²(t,q,p) + dA₃dx₃(t,q) * v³(t,q,p) + E₃(t,q) - μ * dBdx₃(t,q) +
           (dg₁₁dx₃(t,q) * v¹(t,q,p)^2 + dg₂₂dx₃(t,q) * v²(t,q,p)^2 + dg₃₃dx₃(t,q) * v³(t,q,p)^2) / 2
    nothing
end


pauli_particle_3d_iode_ϑ(θ, t, q, v, params) = ϑ(θ, t, q, v)

function pauli_particle_3d_iode_f(f, t, q, v, params)
    @unpack μ = params
    f[1] = dA₁dx₁(t,q) * v[1] + dA₂dx₁(t,q) * v[2] + dA₃dx₁(t,q) * v[3] + E₁(t,q) - μ * dBdx₁(t,q) +
           (dg₁₁dx₁(t,q) * v[1]^2 + dg₂₂dx₁(t,q) * v[2]^2 + dg₃₃dx₁(t,q) * v[3]^2) / 2
    f[2] = dA₁dx₂(t,q) * v[1] + dA₂dx₂(t,q) * v[2] + dA₃dx₂(t,q) * v[3] + E₂(t,q) - μ * dBdx₂(t,q) +
           (dg₁₁dx₂(t,q) * v[1]^2 + dg₂₂dx₂(t,q) * v[2]^2 + dg₃₃dx₂(t,q) * v[3]^2) / 2
    f[3] = dA₁dx₃(t,q) * v[1] + dA₂dx₃(t,q) * v[2] + dA₃dx₃(t,q) * v[3] + E₃(t,q) - μ * dBdx₃(t,q) +
           (dg₁₁dx₃(t,q) * v[1]^2 + dg₂₂dx₃(t,q) * v[2]^2 + dg₃₃dx₃(t,q) * v[3]^2) / 2
    nothing
end

function pauli_particle_3d_iode_g(g, t, q, v, λ, params)
    g[1] = dA₁dx₁(t,q) * v[1] + dA₂dx₁(t,q) * v[2] + dA₃dx₁(t,q) * v[3] +
           (dg₁₁dx₁(t,q) * v[1]^2 + dg₂₂dx₁(t,q) * v[2]^2 + dg₃₃dx₁(t,q) * v[3]^2) / 2
    g[2] = dA₁dx₂(t,q) * v[1] + dA₂dx₂(t,q) * v[2] + dA₃dx₂(t,q) * v[3] +
           (dg₁₁dx₂(t,q) * v[1]^2 + dg₂₂dx₂(t,q) * v[2]^2 + dg₃₃dx₂(t,q) * v[3]^2) / 2
    g[3] = dA₁dx₃(t,q) * v[1] + dA₂dx₃(t,q) * v[2] + dA₃dx₃(t,q) * v[3] +
           (dg₁₁dx₃(t,q) * v[1]^2 + dg₂₂dx₃(t,q) * v[2]^2 + dg₃₃dx₃(t,q) * v[3]^2) / 2
    nothing
end


function pauli_particle_3d_pode(q₀, v₀, parameters; tspan = tspan, tstep = Δt / 100)
    PODEProblem(
        pauli_particle_3d_pode_v,
        pauli_particle_3d_pode_f, 
        tspan, tstep, q₀, pauli_particle_3d_pᵢ(q₀, v₀);
        parameters = parameters,
        invariants = (h = hamiltonian,)
    )
end

pauli_particle_3d_pode(qᵢ=qᵢ, vᵢ=vᵢ; kwargs...) = pauli_particle_3d_pode(initial_conditions(qᵢ, vᵢ)...; kwargs...)


function pauli_particle_3d_hode(q₀, v₀, parameters; tspan = tspan, tstep = Δt / 100)
    HODEProblem(
        pauli_particle_3d_pode_v,
        pauli_particle_3d_pode_f,
        hamiltonian,
        tspan, tstep, q₀, pauli_particle_3d_pᵢ(q₀, v₀);
        parameters = parameters)
end

pauli_particle_3d_hode(qᵢ=qᵢ, vᵢ=vᵢ; kwargs...) = pauli_particle_3d_hode(initial_conditions(qᵢ, vᵢ)...; kwargs...)


function pauli_particle_3d_iode(q₀::AbstractVector, v₀::AbstractVector, parameters::NamedTuple; tspan = tspan, tstep = Δt)
    IODEProblem(
        pauli_particle_3d_iode_ϑ,
        pauli_particle_3d_iode_f,
        pauli_particle_3d_iode_g,
        tspan, tstep, q₀, pauli_particle_3d_pᵢ(q₀, v₀);
        parameters = parameters,
        invariants = (h = hamiltonian,))
end

pauli_particle_3d_iode(q₀::AbstractVector, v₀::AbstractVector, μ::Real; kwargs...) = pauli_particle_3d_iode(q₀, v₀, (μ=μ,); kwargs...)

pauli_particle_3d_iode(qᵢ=qᵢ, vᵢ=vᵢ; kwargs...) = pauli_particle_3d_iode(initial_conditions(qᵢ, vᵢ)...; kwargs...)


# function pauli_particle_3d_lode(q₀::AbstractVector, v₀::AbstractVector, parameters::NamedTuple; tspan = tspan, tstep = Δt)
#     LODEProblem(
#         pauli_particle_3d_iode_ϑ,
#         pauli_particle_3d_iode_f,
#         pauli_particle_3d_iode_g,
#         pauli_particle_3d_ω, lagrangian,
#         tspan, step, q₀, pauli_particle_3d_pᵢ(q₀, v₀);
#         parameters = parameters,
#         invariants = (h = hamiltonian,))
# end

# pauli_particle_3d_lode(q₀::AbstractVector, v₀::AbstractVector, μ::Real) = pauli_particle_3d_lode(q₀, v₀, (μ=μ,))

# pauli_particle_3d_lode(qᵢ=qᵢ, vᵢ=vᵢ) = pauli_particle_3d_lode(initial_conditions(qᵢ, vᵢ)...)
