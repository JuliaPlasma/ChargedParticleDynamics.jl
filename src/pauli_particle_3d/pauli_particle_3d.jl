
using Parameters

ϑ₁(t, q, v) = v[1] + A₁(t, q)
ϑ₂(t, q, v) = v[2] + A₂(t, q)
ϑ₃(t, q, v) = v[3] + A₃(t, q)

v₁(t, q, p) = p[1] - A₁(t, q)
v₂(t, q, p) = p[2] - A₂(t, q)
v₃(t, q, p) = p[3] - A₃(t, q)

ϕ₀(x) = E₀*sin(2π*x)
ϕ(t,q) = ϕ₀(q[3])

E₁(t,q) = zero(eltype(q))
E₂(t,q) = zero(eltype(q))
E₃(t,q) = - 2π*E₀*cos(2π*q[3])

hamiltonian(t,q,p,μ) = 0.5 * (v₁(t,q,p)^2 + v₂(t,q,p)^2 + v₃(t,q,p)^2) + μ * B(t,q) + ϕ(t, q)


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
    v[1] = v₁(t,q,p)
    v[2] = v₂(t,q,p)
    v[3] = v₃(t,q,p)
    nothing
end

function pauli_particle_3d_pode_f(t, q, p, f, params)
    @unpack μ = params
    f[1] = dA₁dx₁(t,q) * v₁(t,q,p) + dA₂dx₁(t,q) * v₂(t,q,p) + dA₃dx₁(t,q) * v₃(t,q,p) + E₁(t,q) - μ * dBdx₁(t,q)
    f[2] = dA₁dx₂(t,q) * v₁(t,q,p) + dA₂dx₂(t,q) * v₂(t,q,p) + dA₃dx₂(t,q) * v₃(t,q,p) + E₂(t,q) - μ * dBdx₂(t,q)
    f[3] = dA₁dx₃(t,q) * v₁(t,q,p) + dA₂dx₃(t,q) * v₂(t,q,p) + dA₃dx₃(t,q) * v₃(t,q,p) + E₃(t,q) - μ * dBdx₃(t,q)
    nothing
end


function pauli_particle_3d_pode(q₀=qᵢ, v₀=vᵢ; params=parameters)
    PODE(pauli_particle_3d_pode_v, pauli_particle_3d_pode_f, q₀, pauli_particle_3d_pᵢ(q₀, v₀);
         parameters=params, h=hamiltonian)
end

