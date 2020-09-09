
using GeometricIntegrators.Equations


ϑ₁(t, q) = q[4] + A₁(t,q)
ϑ₂(t, q) = q[5] + A₂(t,q)
ϑ₃(t, q) = q[6] + A₃(t,q)

function ϑ(t, q, p)
    p[1] = ϑ₁(t,q)
    p[2] = ϑ₂(t,q)
    p[3] = ϑ₃(t,q)
    p[4] = zero(eltype(q))
    p[5] = zero(eltype(q))
    p[6] = zero(eltype(q))
    nothing
end


ϕ₀(x) = E₀*sin(2π*x)
ϕ(t,q) = ϕ₀(q[3])

E₁(t,q) = zero(eltype(q))
E₂(t,q) = zero(eltype(q))
E₃(t,q) = - 2π*E₀*cos(2π*q[3])

hamiltonian(t,q) = 0.5 * (q[4]^2 + q[5]^2 + q[6]^2) + ϕ(t, q)

angular_momentum(t,q) = q[1] * ϑ₂(t,q) - q[2] * ϑ₁(t,q)


function charged_particle_3d_pᵢ(qᵢ, tᵢ=0)
    pᵢ = zero(qᵢ)

    if ndims(qᵢ) == 1
        ϑ(tᵢ, qᵢ, pᵢ)
    else
        for i in 1:size(qᵢ,2)
            @views ϑ(tᵢ, qᵢ[:,i], pᵢ[:,i])
        end
    end
    pᵢ
end


charged_particle_3d_iode_ϑ(t, q, v, p) = ϑ(t, q, p)

function charged_particle_3d_iode_f(t, q, v, f)
    f[1] = dA₁dx₁(t,q) * v[1] + dA₂dx₁(t,q) * v[2] + dA₃dx₁(t,q) * v[3] + E₁(t,q)
    f[2] = dA₁dx₂(t,q) * v[1] + dA₂dx₂(t,q) * v[2] + dA₃dx₂(t,q) * v[3] + E₂(t,q)
    f[3] = dA₁dx₃(t,q) * v[1] + dA₂dx₃(t,q) * v[2] + dA₃dx₃(t,q) * v[3] + E₂(t,q)
    f[4] = v[1] - q[4]
    f[5] = v[2] - q[5]
    f[6] = v[3] - q[6]
    nothing
end

function charged_particle_3d_iode_g(t, q, λ, g)
    g[1] = dA₁dx₁(t,q) * λ[1] + dA₂dx₁(t,q) * λ[2] + dA₃dx₁(t,q) * λ[3]
    g[2] = dA₁dx₂(t,q) * λ[1] + dA₂dx₂(t,q) * λ[2] + dA₃dx₂(t,q) * λ[3]
    g[3] = dA₁dx₃(t,q) * λ[1] + dA₂dx₃(t,q) * λ[2] + dA₃dx₃(t,q) * λ[3]
    g[4] = λ[1]
    g[5] = λ[2]
    g[6] = λ[3]
    nothing
end

function charged_particle_3d_iode_v(t, q, v)
    v[1] = q[4]
    v[2] = q[5]
    v[3] = q[6]
    v[4] = E₁(t,q) + q[5] * B₃(t,q) - q[6] * B₂(t,q)
    v[5] = E₂(t,q) + q[6] * B₁(t,q) - q[4] * B₂(t,q)
    v[6] = E₃(t,q) + q[4] * B₂(t,q) - q[5] * B₁(t,q)
    nothing
end

function charged_particle_3d_iode(q₀=q₀)
    IODE(charged_particle_3d_iode_ϑ, charged_particle_3d_iode_f,
         charged_particle_3d_iode_g, q₀, charged_particle_3d_pᵢ(q₀);
         h=hamiltonian, v=charged_particle_3d_iode_v)
end
