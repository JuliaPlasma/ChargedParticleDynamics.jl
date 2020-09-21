
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

hamiltonian(t,q,p) = 0.5 * (g₁₁(t,q) * v¹(t,q,p)^2 + g₂₂(t,q) * v²(t,q,p)^2 + g₃₃(t,q) * v³(t,q,p)^2)
toroidal_momentum(t,q,p) = p[3]


function charged_particle_3d_pᵢ(qᵢ, vᵢ, tᵢ=0)
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


function charged_particle_3d_pode_v(t, q, p, v)
    v[1] = v¹(t,q,p)
    v[2] = v²(t,q,p)
    v[3] = v³(t,q,p)
    nothing
end

function charged_particle_3d_pode_f(t, q, p, f)
    f[1] = dA₁dx₁(t,q) * v¹(t,q,p) + dA₂dx₁(t,q) * v²(t,q,p) + dA₃dx₁(t,q) * v³(t,q,p) + E₁(t,q)
    f[2] = dA₁dx₂(t,q) * v¹(t,q,p) + dA₂dx₂(t,q) * v²(t,q,p) + dA₃dx₂(t,q) * v³(t,q,p) + E₂(t,q)
    f[3] = dA₁dx₃(t,q) * v¹(t,q,p) + dA₂dx₃(t,q) * v²(t,q,p) + dA₃dx₃(t,q) * v³(t,q,p) + E₃(t,q)
    nothing
end


charged_particle_3d_iode_ϑ(t, q, v, θ) = ϑ(t, q, v, θ)

function charged_particle_3d_iode_f(t, q, v, f)
    f[1] = dA₁dx₁(t,q) * v[1] + dA₂dx₁(t,q) * v[2] + dA₃dx₁(t,q) * v[3] + E₁(t,q)
    f[2] = dA₁dx₂(t,q) * v[1] + dA₂dx₂(t,q) * v[2] + dA₃dx₂(t,q) * v[3] + E₂(t,q)
    f[3] = dA₁dx₃(t,q) * v[1] + dA₂dx₃(t,q) * v[2] + dA₃dx₃(t,q) * v[3] + E₃(t,q)
    nothing
end

function charged_particle_3d_iode_g(t, q, v, f)
    f[1] = dA₁dx₁(t,q) * v[1] + dA₂dx₁(t,q) * v[2] + dA₃dx₁(t,q) * v[3]
    f[2] = dA₁dx₂(t,q) * v[1] + dA₂dx₂(t,q) * v[2] + dA₃dx₂(t,q) * v[3]
    f[3] = dA₁dx₃(t,q) * v[1] + dA₂dx₃(t,q) * v[2] + dA₃dx₃(t,q) * v[3]
    nothing
end


# function charged_particle_3d_pode(q₀=qᵢ, v₀=vᵢ)
#     PODE(charged_particle_3d_pode_v, charged_particle_3d_pode_f, q₀, charged_particle_3d_pᵢ(q₀, v₀))
# end
function charged_particle_3d_pode(q₀=qᵢ, p₀=pᵢ)
    PODE(charged_particle_3d_pode_v, charged_particle_3d_pode_f, q₀, p₀)
end


function charged_particle_3d_iode(q₀=qᵢ, p₀=pᵢ)
    IODE(charged_particle_3d_iode_ϑ, charged_particle_3d_iode_f, charged_particle_3d_iode_g,
            q₀, p₀;
            h=hamiltonian, v = (t, q, v) -> nothing)
end