
using Parameters

export hamiltonian


@inline function u(t, q)
    q[4]
end


function hamiltonian(t, q, params)
    @unpack μ = params
    0.5 * u(t,q)^2 + μ*B(t,q)
end


dHdx₁(t,q,μ) = μ * dBdx₁(t,q)
dHdx₂(t,q,μ) = μ * dBdx₂(t,q)
dHdx₃(t,q,μ) = μ * dBdx₃(t,q)
dHdx₄(t,q,μ) = u(t,q)

function dH(t, q, dH, params)
    @unpack μ = params
    dH[1] = dHdx₁(t,q,μ)
    dH[2] = dHdx₂(t,q,μ)
    dH[3] = dHdx₃(t,q,μ)
    dH[4] = dHdx₄(t,q,μ)
    nothing
end


ϑ₁(t,q) = A₁(t,q) + u(t,q) * b₁(t,q)
ϑ₂(t,q) = A₂(t,q) + u(t,q) * b₂(t,q)
ϑ₃(t,q) = A₃(t,q) + u(t,q) * b₃(t,q)
ϑ₄(t,q) = zero(eltype(q))

function ϑ(t, q, ϑ, params)
    ϑ[1] = ϑ₁(t,q)
    ϑ[2] = ϑ₂(t,q)
    ϑ[3] = ϑ₃(t,q)
    ϑ[4] = ϑ₄(t,q)
    nothing
end


β₁(t, q) = u(t,q) * ϑ₁(t,q)
β₂(t, q) = u(t,q) * ϑ₂(t,q)
β₃(t, q) = u(t,q) * ϑ₃(t,q)

function β(t, q, β, params)
    β[1] = β₁(t,q)
    β[2] = β₂(t,q)
    β[3] = β₃(t,q)
    nothing
end


γ₁(t, q, μ) = ϑ₂(t,q) * dHdx₃(t,q,μ) - dHdx₂(t,q,μ) * ϑ₃(t,q)
γ₂(t, q, μ) = ϑ₃(t,q) * dHdx₁(t,q,μ) - dHdx₃(t,q,μ) * ϑ₁(t,q)
γ₃(t, q, μ) = ϑ₁(t,q) * dHdx₂(t,q,μ) - dHdx₁(t,q,μ) * ϑ₂(t,q)

function γ(t, q, γ, params)
    @unpack μ = params
    γ[1] = γ₁(t,q,μ)
    γ[2] = γ₂(t,q,μ)
    γ[3] = γ₃(t,q,μ)
    nothing
end


d²Hdx₁dx₁(t, q, μ) = μ * d²Bdx₁dx₁(t,q)
d²Hdx₁dx₂(t, q, μ) = μ * d²Bdx₁dx₂(t,q)
d²Hdx₁dx₃(t, q, μ) = μ * d²Bdx₁dx₃(t,q)
d²Hdx₂dx₁(t, q, μ) = μ * d²Bdx₂dx₁(t,q)
d²Hdx₂dx₂(t, q, μ) = μ * d²Bdx₂dx₂(t,q)
d²Hdx₂dx₃(t, q, μ) = μ * d²Bdx₂dx₃(t,q)
d²Hdx₃dx₁(t, q, μ) = μ * d²Bdx₃dx₁(t,q)
d²Hdx₃dx₂(t, q, μ) = μ * d²Bdx₃dx₂(t,q)
d²Hdx₃dx₃(t, q, μ) = μ * d²Bdx₃dx₃(t,q)

dϑ₁dx₁(t, q) = dA₁dx₁(t,q) + u(t,q) * db₁dx₁(t,q)
dϑ₁dx₂(t, q) = dA₁dx₂(t,q) + u(t,q) * db₁dx₂(t,q)
dϑ₁dx₃(t, q) = dA₁dx₃(t,q) + u(t,q) * db₁dx₃(t,q)
dϑ₁dx₄(t, q) = b₁(t,q)
dϑ₂dx₁(t, q) = dA₂dx₁(t,q) + u(t,q) * db₂dx₁(t,q)
dϑ₂dx₂(t, q) = dA₂dx₂(t,q) + u(t,q) * db₂dx₂(t,q)
dϑ₂dx₃(t, q) = dA₂dx₃(t,q) + u(t,q) * db₂dx₃(t,q)
dϑ₂dx₄(t, q) = b₂(t,q)
dϑ₃dx₁(t, q) = dA₃dx₁(t,q) + u(t,q) * db₃dx₁(t,q)
dϑ₃dx₂(t, q) = dA₃dx₂(t,q) + u(t,q) * db₃dx₂(t,q)
dϑ₃dx₃(t, q) = dA₃dx₃(t,q) + u(t,q) * db₃dx₃(t,q)
dϑ₃dx₄(t, q) = b₃(t,q)
dϑ₄dx₁(t, q) = zero(eltype(q))
dϑ₄dx₂(t, q) = zero(eltype(q))
dϑ₄dx₃(t, q) = zero(eltype(q))
dϑ₄dx₄(t, q) = zero(eltype(q))

dβ₁dx₁(t, q) = u(t,q) * dϑ₁dx₁(t,q)
dβ₁dx₂(t, q) = u(t,q) * dϑ₁dx₂(t,q)
dβ₁dx₃(t, q) = u(t,q) * dϑ₁dx₃(t,q)
dβ₁dx₄(t, q) = ϑ₁(t,q)
dβ₂dx₁(t, q) = u(t,q) * dϑ₂dx₁(t,q)
dβ₂dx₂(t, q) = u(t,q) * dϑ₂dx₂(t,q)
dβ₂dx₃(t, q) = u(t,q) * dϑ₂dx₃(t,q)
dβ₂dx₄(t, q) = ϑ₂(t,q)
dβ₃dx₁(t, q) = u(t,q) * dϑ₃dx₁(t,q)
dβ₃dx₂(t, q) = u(t,q) * dϑ₃dx₂(t,q)
dβ₃dx₃(t, q) = u(t,q) * dϑ₃dx₃(t,q)
dβ₃dx₄(t, q) = ϑ₃(t,q)

dγ₁dx₁(t, q, μ) = dϑ₂dx₁(t,q) * dHdx₃(t,q,μ) + ϑ₂(t,q) * d²Hdx₃dx₁(t,q,μ) - d²Hdx₂dx₁(t,q,μ) * ϑ₃(t,q) - dHdx₂(t,q,μ) * dϑ₃dx₁(t,q)
dγ₁dx₂(t, q, μ) = dϑ₂dx₂(t,q) * dHdx₃(t,q,μ) + ϑ₂(t,q) * d²Hdx₃dx₂(t,q,μ) - d²Hdx₂dx₂(t,q,μ) * ϑ₃(t,q) - dHdx₂(t,q,μ) * dϑ₃dx₂(t,q)
dγ₁dx₃(t, q, μ) = dϑ₂dx₃(t,q) * dHdx₃(t,q,μ) + ϑ₂(t,q) * d²Hdx₃dx₃(t,q,μ) - d²Hdx₂dx₃(t,q,μ) * ϑ₃(t,q) - dHdx₂(t,q,μ) * dϑ₃dx₃(t,q)
dγ₁dx₄(t, q, μ) = b₂(t,q) * dHdx₃(t,q,μ) - dHdx₂(t,q,μ) * b₃(t,q)
dγ₂dx₁(t, q, μ) = dϑ₃dx₁(t,q) * dHdx₁(t,q,μ) + ϑ₃(t,q) * d²Hdx₁dx₁(t,q,μ) - d²Hdx₃dx₁(t,q,μ) * ϑ₁(t,q) - dHdx₃(t,q,μ) * dϑ₁dx₁(t,q)
dγ₂dx₂(t, q, μ) = dϑ₃dx₂(t,q) * dHdx₁(t,q,μ) + ϑ₃(t,q) * d²Hdx₁dx₂(t,q,μ) - d²Hdx₃dx₂(t,q,μ) * ϑ₁(t,q) - dHdx₃(t,q,μ) * dϑ₁dx₂(t,q)
dγ₂dx₃(t, q, μ) = dϑ₃dx₃(t,q) * dHdx₁(t,q,μ) + ϑ₃(t,q) * d²Hdx₁dx₃(t,q,μ) - d²Hdx₃dx₃(t,q,μ) * ϑ₁(t,q) - dHdx₃(t,q,μ) * dϑ₁dx₃(t,q)
dγ₂dx₄(t, q, μ) = b₃(t,q) * dHdx₁(t,q,μ) - dHdx₃(t,q,μ) * b₁(t,q)
dγ₃dx₁(t, q, μ) = dϑ₁dx₁(t,q) * dHdx₂(t,q,μ) + ϑ₁(t,q) * d²Hdx₂dx₁(t,q,μ) - d²Hdx₁dx₁(t,q,μ) * ϑ₂(t,q) - dHdx₁(t,q,μ) * dϑ₂dx₁(t,q)
dγ₃dx₂(t, q, μ) = dϑ₁dx₂(t,q) * dHdx₂(t,q,μ) + ϑ₁(t,q) * d²Hdx₂dx₂(t,q,μ) - d²Hdx₁dx₂(t,q,μ) * ϑ₂(t,q) - dHdx₁(t,q,μ) * dϑ₂dx₂(t,q)
dγ₃dx₃(t, q, μ) = dϑ₁dx₃(t,q) * dHdx₂(t,q,μ) + ϑ₁(t,q) * d²Hdx₂dx₃(t,q,μ) - d²Hdx₁dx₃(t,q,μ) * ϑ₂(t,q) - dHdx₁(t,q,μ) * dϑ₂dx₃(t,q)
dγ₃dx₄(t, q, μ) = b₁(t,q) * dHdx₂(t,q,μ) - dHdx₁(t,q,μ) * b₂(t,q)



function v₁(t, q, v, params)
    v[1] = + dβ₃dx₂(t,q)
    v[2] = - dβ₃dx₁(t,q)
    v[3] = 0
    v[4] = 0
end

function v₂(t, q, v, params)
    v[1] = - dβ₂dx₃(t,q)
    v[2] = 0
    v[3] = + dβ₂dx₁(t,q)
    v[4] = 0
end

function v₃(t, q, v, params)
    v[1] = 0
    v[2] = + dβ₁dx₃(t,q)
    v[3] = - dβ₁dx₂(t,q)
    v[4] = 0
end

function v₄(t, q, v, params)
    @unpack μ = params
    v[1] = + dγ₁dx₄(t,q,μ)
    v[2] = 0
    v[3] = 0
    v[4] = - dγ₁dx₁(t,q,μ)
end

function v₅(t, q, v, params)
    @unpack μ = params
    v[1] = 0
    v[2] = + dγ₂dx₄(t,q,μ)
    v[3] = 0
    v[4] = - dγ₂dx₂(t,q,μ)
end

function v₆(t, q, v, params)
    @unpack μ = params
    v[1] = 0
    v[2] = 0
    v[3] = + dγ₃dx₄(t,q,μ)
    v[4] = - dγ₃dx₃(t,q,μ)
end

function v(t, q, v, params)
    @unpack μ = params
    v[1] =   dγ₁dx₄(t,q,μ) + dβ₃dx₂(t,q)   - dβ₂dx₃(t,q)
    v[2] =   dγ₂dx₄(t,q,μ) + dβ₁dx₃(t,q)   - dβ₃dx₁(t,q)
    v[3] =   dγ₃dx₄(t,q,μ) + dβ₂dx₁(t,q)   - dβ₁dx₂(t,q)
    v[4] = - dγ₁dx₁(t,q,μ) - dγ₂dx₂(t,q,μ) - dγ₃dx₃(t,q,μ)
    nothing
end
