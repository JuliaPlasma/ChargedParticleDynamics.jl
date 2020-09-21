
using Parameters
using GeometricIntegrators.Equations


ϑ₁(t, q, v) = g₁₁(t,q) * v[1] + A₁(t, q)
ϑ₂(t, q, v) = g₂₂(t,q) * v[2] + A₂(t, q)
ϑ₃(t, q, v) = g₃₃(t,q) * v[3] + A₃(t, q)

function ϑ(t, q, v, θ)
    θ[1] = ϑ₁(t,q,v)
    θ[2] = ϑ₂(t,q,v)
    θ[3] = ϑ₃(t,q,v)
    nothing
end

v¹(t, q, p) = g¹¹(t, q) * (p[1] - A₁(t, q))
v²(t, q, p) = g²²(t, q) * (p[2] - A₂(t, q))
v³(t, q, p) = g³³(t, q) * (p[3] - A₃(t, q))

ϕ₀(x) = E₀*sin(2π*x)
ϕ(t,q) = ϕ₀(q[3])

E₁(t,q) = zero(eltype(q))
E₂(t,q) = zero(eltype(q))
E₃(t,q) = - 2π*E₀*cos(2π*q[3])

function hamiltonian(t, q, p, params)
    @unpack μ = params
    0.5 * (g₁₁(t,q) * v¹(t,q,p)^2 + g₂₂(t,q) * v²(t,q,p)^2 + g₃₃(t,q) * v³(t,q,p)^2) + μ*B(t,q)
end


function initial_conditions(x₀, v₀)
    b(t,x) = [b₁(t,x), b₂(t,x), b₃(t,x)]
    b⃗(t,x) = [b¹(t,x), b²(t,x), b³(t,x)]

    u₀ = v₀' * b(0, x₀)
    vpar = u₀ .* b⃗(0, x₀)
    vper = v₀ .- vpar
    μ = vper' * vper / 2 / B(0, x₀)
    
    (x₀, vpar, (μ = μ,))
end


function pauli_particle_3d_pᵢ(qᵢ, vᵢ, tᵢ=0)
    pᵢ = zero(qᵢ)

    if ndims(qᵢ) == 1
        ϑ(tᵢ, qᵢ, vᵢ, pᵢ)
    else
        for i in 1:size(qᵢ,2)
            @views ϑ(tᵢ, qᵢ[:,i], vᵢ[:,i], pᵢ[:,i])
        end
    end
    pᵢ
end


function pauli_particle_3d_pode_v(t, q, p, v, params)
    v[1] = v¹(t,q,p)
    v[2] = v²(t,q,p)
    v[3] = v³(t,q,p)
    nothing
end

function pauli_particle_3d_pode_f(t, q, p, f, params)
    @unpack μ = params
    f[1] = dA₁dx₁(t,q) * v¹(t,q,p) + dA₂dx₁(t,q) * v²(t,q,p) + dA₃dx₁(t,q) * v³(t,q,p) + E₁(t,q) - μ * dBdx₁(t,q)
    f[2] = dA₁dx₂(t,q) * v¹(t,q,p) + dA₂dx₂(t,q) * v²(t,q,p) + dA₃dx₂(t,q) * v³(t,q,p) + E₂(t,q) - μ * dBdx₂(t,q)
    f[3] = dA₁dx₃(t,q) * v¹(t,q,p) + dA₂dx₃(t,q) * v²(t,q,p) + dA₃dx₃(t,q) * v³(t,q,p) + E₃(t,q) - μ * dBdx₃(t,q)
    nothing
end


pauli_particle_3d_iode_ϑ(t, q, v, θ, params) = ϑ(t, q, v, θ)

function pauli_particle_3d_iode_f(t, q, v, f, params)
    @unpack μ = params
    f[1] = dA₁dx₁(t,q) * v[1] + dA₂dx₁(t,q) * v[2] + dA₃dx₁(t,q) * v[3] + E₁(t,q) - μ * dBdx₁(t,q)
    f[2] = dA₁dx₂(t,q) * v[1] + dA₂dx₂(t,q) * v[2] + dA₃dx₂(t,q) * v[3] + E₂(t,q) - μ * dBdx₂(t,q)
    f[3] = dA₁dx₃(t,q) * v[1] + dA₂dx₃(t,q) * v[2] + dA₃dx₃(t,q) * v[3] + E₃(t,q) - μ * dBdx₃(t,q)
    nothing
end

function pauli_particle_3d_iode_g(t, q, v, f, params)
    f[1] = dA₁dx₁(t,q) * v[1] + dA₂dx₁(t,q) * v[2] + dA₃dx₁(t,q) * v[3]
    f[2] = dA₁dx₂(t,q) * v[1] + dA₂dx₂(t,q) * v[2] + dA₃dx₂(t,q) * v[3]
    f[3] = dA₁dx₃(t,q) * v[1] + dA₂dx₃(t,q) * v[2] + dA₃dx₃(t,q) * v[3]
    nothing
end


function pauli_particle_3d_pode(q₀, v₀, params)
    PODE(pauli_particle_3d_pode_v, pauli_particle_3d_pode_f, q₀, pauli_particle_3d_pᵢ(q₀, v₀);
         parameters=params, h=hamiltonian)
end

pauli_particle_3d_pode(qᵢ=qᵢ, vᵢ=vᵢ) = pauli_particle_3d_pode(initial_conditions(qᵢ, vᵢ)...)

function pauli_particle_3d_iode(q₀, v₀, params)
    IODE(pauli_particle_3d_iode_ϑ, pauli_particle_3d_iode_f, pauli_particle_3d_iode_g,
            q₀, pauli_particle_3d_pᵢ(q₀, v₀);
            parameters=params, h=hamiltonian, v = (t, q, v, params) -> nothing)
end

pauli_particle_3d_iode(qᵢ=qᵢ, vᵢ=vᵢ) = pauli_particle_3d_iode(initial_conditions(qᵢ, vᵢ)...)

