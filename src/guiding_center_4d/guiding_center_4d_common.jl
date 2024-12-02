
using Parameters

export hamiltonian, u, ω, ϑ, ϑ₁, ϑ₂, ϑ₃, ϑ₄, dϑ, β₁, β₂, β₃, B, B₁, B₂, B₃, b₁, b₂, b₃, dH


@inline u(t,q) = q[4]

ϑ₁(t,q) = A₁(t,q) + u(t,q) * b₁(t,q)
ϑ₂(t,q) = A₂(t,q) + u(t,q) * b₂(t,q)
ϑ₃(t,q) = A₃(t,q) + u(t,q) * b₃(t,q)
ϑ₄(t,q) = zero(eltype(q))

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


function ϑ(θ::AbstractVector, t::Number, q::AbstractVector)
    θ[1] = ϑ₁(t,q)
    θ[2] = ϑ₂(t,q)
    θ[3] = ϑ₃(t,q)
    θ[4] = ϑ₄(t,q)
    nothing
end

function ϑ(t::Number, q::AbstractVector, k::Int)
    if k == 1
        ϑ₁(t, q)
    elseif k == 2
        ϑ₂(t, q)
    elseif k == 3
        ϑ₃(t, q)
    elseif k == 4
        ϑ₄(t, q)
    else
        nothing
    end
end


function ω(Ω, t, q)
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


function D²ϑ₁(D²ϑ, t, q)
    D²ϑ[1,1] = d²A₁dx₁dx₁(t,q) + u(t,q) * d²b₁dx₁dx₁(t,q)
    D²ϑ[1,2] = d²A₁dx₁dx₂(t,q) + u(t,q) * d²b₁dx₁dx₂(t,q)
    D²ϑ[1,3] = d²A₁dx₁dx₃(t,q) + u(t,q) * d²b₁dx₁dx₃(t,q)
    D²ϑ[1,4] = db₁dx₁(t,q)

    D²ϑ[2,1] = d²A₁dx₂dx₁(t,q) + u(t,q) * d²b₁dx₂dx₁(t,q)
    D²ϑ[2,2] = d²A₁dx₂dx₂(t,q) + u(t,q) * d²b₁dx₂dx₂(t,q)
    D²ϑ[2,3] = d²A₁dx₂dx₃(t,q) + u(t,q) * d²b₁dx₂dx₃(t,q)
    D²ϑ[2,4] = db₁dx₂(t,q)

    D²ϑ[3,1] = d²A₁dx₃dx₁(t,q) + u(t,q) * d²b₁dx₃dx₁(t,q)
    D²ϑ[3,2] = d²A₁dx₃dx₂(t,q) + u(t,q) * d²b₁dx₃dx₂(t,q)
    D²ϑ[3,3] = d²A₁dx₃dx₃(t,q) + u(t,q) * d²b₁dx₃dx₃(t,q)
    D²ϑ[3,4] = db₁dx₃(t,q)

    D²ϑ[4,1] = db₁dx₁(t,q)
    D²ϑ[4,2] = db₁dx₂(t,q)
    D²ϑ[4,3] = db₁dx₃(t,q)
    D²ϑ[4,4] = 0

    nothing
end

function D²ϑ₂(D²ϑ, t, q)
    D²ϑ[1,1] = d²A₂dx₁dx₁(t,q) + u(t,q) * d²b₂dx₁dx₁(t,q)
    D²ϑ[1,2] = d²A₂dx₁dx₂(t,q) + u(t,q) * d²b₂dx₁dx₂(t,q)
    D²ϑ[1,3] = d²A₂dx₁dx₃(t,q) + u(t,q) * d²b₂dx₁dx₃(t,q)
    D²ϑ[1,4] = db₂dx₁(t,q)

    D²ϑ[2,1] = d²A₂dx₂dx₁(t,q) + u(t,q) * d²b₂dx₂dx₁(t,q)
    D²ϑ[2,2] = d²A₂dx₂dx₂(t,q) + u(t,q) * d²b₂dx₂dx₂(t,q)
    D²ϑ[2,3] = d²A₂dx₂dx₃(t,q) + u(t,q) * d²b₂dx₂dx₃(t,q)
    D²ϑ[2,4] = db₂dx₂(t,q)

    D²ϑ[3,1] = d²A₂dx₃dx₁(t,q) + u(t,q) * d²b₂dx₃dx₁(t,q)
    D²ϑ[3,2] = d²A₂dx₃dx₂(t,q) + u(t,q) * d²b₂dx₃dx₂(t,q)
    D²ϑ[3,3] = d²A₂dx₃dx₃(t,q) + u(t,q) * d²b₂dx₃dx₃(t,q)
    D²ϑ[3,4] = db₂dx₃(t,q)

    D²ϑ[4,1] = db₂dx₁(t,q)
    D²ϑ[4,2] = db₂dx₂(t,q)
    D²ϑ[4,3] = db₂dx₃(t,q)
    D²ϑ[4,4] = 0

    nothing
end

function D²ϑ₃(D²ϑ, t, q)
    D²ϑ[1,1] = d²A₃dx₁dx₁(t,q) + u(t,q) * d²b₃dx₁dx₁(t,q)
    D²ϑ[1,2] = d²A₃dx₁dx₂(t,q) + u(t,q) * d²b₃dx₁dx₂(t,q)
    D²ϑ[1,3] = d²A₃dx₁dx₃(t,q) + u(t,q) * d²b₃dx₁dx₃(t,q)
    D²ϑ[1,4] = db₃dx₁(t,q)

    D²ϑ[2,1] = d²A₃dx₂dx₁(t,q) + u(t,q) * d²b₃dx₂dx₁(t,q)
    D²ϑ[2,2] = d²A₃dx₂dx₂(t,q) + u(t,q) * d²b₃dx₂dx₂(t,q)
    D²ϑ[2,3] = d²A₃dx₂dx₃(t,q) + u(t,q) * d²b₃dx₂dx₃(t,q)
    D²ϑ[2,4] = db₃dx₂(t,q)

    D²ϑ[3,1] = d²A₃dx₃dx₁(t,q) + u(t,q) * d²b₃dx₃dx₁(t,q)
    D²ϑ[3,2] = d²A₃dx₃dx₂(t,q) + u(t,q) * d²b₃dx₃dx₂(t,q)
    D²ϑ[3,3] = d²A₃dx₃dx₃(t,q) + u(t,q) * d²b₃dx₃dx₃(t,q)
    D²ϑ[3,4] = db₃dx₃(t,q)

    D²ϑ[4,1] = db₃dx₁(t,q)
    D²ϑ[4,2] = db₃dx₂(t,q)
    D²ϑ[4,3] = db₃dx₃(t,q)
    D²ϑ[4,4] = 0

    nothing
end

function D²ϑ₄(D²ϑ, t, q)
    D²ϑ .= 0
    nothing
end


function D²ϑd₁(D²ϑ, t, q)
    D²ϑ[1,1] = d²A₁dx₁dx₁(t,q) + u(t,q) * d²b₁dx₁dx₁(t,q)
    D²ϑ[1,2] = d²A₁dx₁dx₂(t,q) + u(t,q) * d²b₁dx₁dx₂(t,q)
    D²ϑ[1,3] = d²A₁dx₁dx₃(t,q) + u(t,q) * d²b₁dx₁dx₃(t,q)
    D²ϑ[1,4] = db₁dx₁(t,q)

    D²ϑ[2,1] = d²A₂dx₁dx₁(t,q) + u(t,q) * d²b₂dx₁dx₁(t,q)
    D²ϑ[2,2] = d²A₂dx₁dx₂(t,q) + u(t,q) * d²b₂dx₁dx₂(t,q)
    D²ϑ[2,3] = d²A₂dx₁dx₃(t,q) + u(t,q) * d²b₂dx₁dx₃(t,q)
    D²ϑ[2,4] = db₂dx₁(t,q)

    D²ϑ[3,1] = d²A₃dx₁dx₁(t,q) + u(t,q) * d²b₃dx₁dx₁(t,q)
    D²ϑ[3,2] = d²A₃dx₁dx₂(t,q) + u(t,q) * d²b₃dx₁dx₂(t,q)
    D²ϑ[3,3] = d²A₃dx₁dx₃(t,q) + u(t,q) * d²b₃dx₁dx₃(t,q)
    D²ϑ[3,4] = db₃dx₁(t,q)

    D²ϑ[4,1] = 0
    D²ϑ[4,2] = 0
    D²ϑ[4,3] = 0
    D²ϑ[4,4] = 0

    nothing
end

function D²ϑd₂(D²ϑ, t, q)
    D²ϑ[2,1] = d²A₁dx₂dx₁(t,q) + u(t,q) * d²b₁dx₂dx₁(t,q)
    D²ϑ[2,2] = d²A₁dx₂dx₂(t,q) + u(t,q) * d²b₁dx₂dx₂(t,q)
    D²ϑ[2,3] = d²A₁dx₂dx₃(t,q) + u(t,q) * d²b₁dx₂dx₃(t,q)
    D²ϑ[2,4] = db₁dx₂(t,q)

    D²ϑ[2,1] = d²A₂dx₂dx₁(t,q) + u(t,q) * d²b₂dx₂dx₁(t,q)
    D²ϑ[2,2] = d²A₂dx₂dx₂(t,q) + u(t,q) * d²b₂dx₂dx₂(t,q)
    D²ϑ[2,3] = d²A₂dx₂dx₃(t,q) + u(t,q) * d²b₂dx₂dx₃(t,q)
    D²ϑ[2,4] = db₂dx₂(t,q)

    D²ϑ[3,1] = d²A₂dx₃dx₁(t,q) + u(t,q) * d²b₂dx₃dx₁(t,q)
    D²ϑ[3,2] = d²A₂dx₃dx₂(t,q) + u(t,q) * d²b₂dx₃dx₂(t,q)
    D²ϑ[3,3] = d²A₂dx₃dx₃(t,q) + u(t,q) * d²b₂dx₃dx₃(t,q)
    D²ϑ[3,4] = db₂dx₃(t,q)

    D²ϑ[4,1] = 0
    D²ϑ[4,2] = 0
    D²ϑ[4,3] = 0
    D²ϑ[4,4] = 0

    nothing
end

function D²ϑd₃(D²ϑ, t, q)
    D²ϑ[2,1] = d²A₁dx₃dx₁(t,q) + u(t,q) * d²b₁dx₃dx₁(t,q)
    D²ϑ[2,2] = d²A₁dx₃dx₂(t,q) + u(t,q) * d²b₁dx₃dx₂(t,q)
    D²ϑ[2,3] = d²A₁dx₃dx₃(t,q) + u(t,q) * d²b₁dx₃dx₃(t,q)
    D²ϑ[2,4] = db₁dx₃(t,q)

    D²ϑ[2,1] = d²A₃dx₂dx₁(t,q) + u(t,q) * d²b₃dx₂dx₁(t,q)
    D²ϑ[2,2] = d²A₃dx₂dx₂(t,q) + u(t,q) * d²b₃dx₂dx₂(t,q)
    D²ϑ[2,3] = d²A₃dx₂dx₃(t,q) + u(t,q) * d²b₃dx₂dx₃(t,q)
    D²ϑ[2,4] = db₃dx₂(t,q)

    D²ϑ[3,1] = d²A₃dx₃dx₁(t,q) + u(t,q) * d²b₃dx₃dx₁(t,q)
    D²ϑ[3,2] = d²A₃dx₃dx₂(t,q) + u(t,q) * d²b₃dx₃dx₂(t,q)
    D²ϑ[3,3] = d²A₃dx₃dx₃(t,q) + u(t,q) * d²b₃dx₃dx₃(t,q)
    D²ϑ[3,4] = db₃dx₃(t,q)

    D²ϑ[4,1] = 0
    D²ϑ[4,2] = 0
    D²ϑ[4,3] = 0
    D²ϑ[4,4] = 0

    nothing
end

function D²ϑd₄(D²ϑ, t, q)
    D²ϑ[1,1] = db₁dx₁(t,q)
    D²ϑ[1,2] = db₁dx₂(t,q)
    D²ϑ[1,3] = db₁dx₃(t,q)
    D²ϑ[1,4] = 0

    D²ϑ[2,1] = db₂dx₁(t,q)
    D²ϑ[2,2] = db₂dx₂(t,q)
    D²ϑ[2,3] = db₂dx₃(t,q)
    D²ϑ[2,4] = 0

    D²ϑ[3,1] = db₃dx₁(t,q)
    D²ϑ[3,2] = db₃dx₂(t,q)
    D²ϑ[3,3] = db₃dx₃(t,q)
    D²ϑ[3,4] = 0

    D²ϑ[4,1] = 0
    D²ϑ[4,2] = 0
    D²ϑ[4,3] = 0
    D²ϑ[4,4] = 0

    nothing
end


function dϑ(dϑ, t, q)
    dϑ[1,1] = dϑ₁dx₁(t,q)
    dϑ[1,2] = dϑ₁dx₂(t,q)
    dϑ[1,3] = dϑ₁dx₃(t,q)
    dϑ[1,4] = dϑ₁dx₄(t,q)

    dϑ[2,1] = dϑ₂dx₁(t,q)
    dϑ[2,2] = dϑ₂dx₂(t,q)
    dϑ[2,3] = dϑ₂dx₃(t,q)
    dϑ[2,4] = dϑ₂dx₄(t,q)

    dϑ[3,1] = dϑ₃dx₁(t,q)
    dϑ[3,2] = dϑ₃dx₂(t,q)
    dϑ[3,3] = dϑ₃dx₃(t,q)
    dϑ[3,4] = dϑ₃dx₄(t,q)

    dϑ[4,1] = dϑ₄dx₁(t,q)
    dϑ[4,2] = dϑ₄dx₂(t,q)
    dϑ[4,3] = dϑ₄dx₃(t,q)
    dϑ[4,4] = dϑ₄dx₄(t,q)

    nothing
end


β₁(t,q) = dϑ₃dx₂(t,q) - dϑ₂dx₃(t,q)
β₂(t,q) = dϑ₁dx₃(t,q) - dϑ₃dx₁(t,q)
β₃(t,q) = dϑ₂dx₁(t,q) - dϑ₁dx₂(t,q)

# function β(t,q)
#    return sqrt(β1(t,q)^2 + β2(t,q)^2 + β3(t,q)^2)
# end


function hamiltonian(t, q, params)
    @unpack μ = params
    0.5 * u(t,q)^2 + μ*B(t,q)
end

hamiltonian(t,q,p,params) = hamiltonian(t,q, params)
lagrangian(t,q,v,params) = ϑ₁(t, q) * v[1] + ϑ₂(t, q) * v[2] + ϑ₃(t, q) * v[3] - hamiltonian(t,q,params)


dHdx₁(t, q, μ) = μ * dBdx₁(t,q)
dHdx₂(t, q, μ) = μ * dBdx₂(t,q)
dHdx₃(t, q, μ) = μ * dBdx₃(t,q)
dHdx₄(t, q, μ) = u(t,q)

function dH(dH, t, q, params)
    @unpack μ = params
    dH[1] = dHdx₁(t, q, μ)
    dH[2] = dHdx₂(t, q, μ)
    dH[3] = dHdx₃(t, q, μ)
    dH[4] = dHdx₄(t, q, μ)
    nothing
end


f₁(t,q,v) = dϑ₁dx₁(t,q) * v[1] + dϑ₂dx₁(t,q) * v[2] + dϑ₃dx₁(t,q) * v[3] + dϑ₄dx₁(t,q) * v[4]
f₂(t,q,v) = dϑ₁dx₂(t,q) * v[1] + dϑ₂dx₂(t,q) * v[2] + dϑ₃dx₂(t,q) * v[3] + dϑ₄dx₂(t,q) * v[4]
f₃(t,q,v) = dϑ₁dx₃(t,q) * v[1] + dϑ₂dx₃(t,q) * v[2] + dϑ₃dx₃(t,q) * v[3] + dϑ₄dx₃(t,q) * v[4]
f₄(t,q,v) = dϑ₁dx₄(t,q) * v[1] + dϑ₂dx₄(t,q) * v[2] + dϑ₃dx₄(t,q) * v[3] + dϑ₄dx₄(t,q) * v[4]

g₁(t,q,v) = dϑ₁dx₁(t,q) * v[1] + dϑ₁dx₂(t,q) * v[2] + dϑ₁dx₃(t,q) * v[3] + dϑ₁dx₄(t,q) * v[4]
g₂(t,q,v) = dϑ₂dx₁(t,q) * v[1] + dϑ₂dx₂(t,q) * v[2] + dϑ₂dx₃(t,q) * v[3] + dϑ₂dx₄(t,q) * v[4]
g₃(t,q,v) = dϑ₃dx₁(t,q) * v[1] + dϑ₃dx₂(t,q) * v[2] + dϑ₃dx₃(t,q) * v[3] + dϑ₃dx₄(t,q) * v[4]
g₄(t,q,v) = dϑ₄dx₁(t,q) * v[1] + dϑ₄dx₂(t,q) * v[2] + dϑ₄dx₃(t,q) * v[3] + dϑ₄dx₄(t,q) * v[4]


function g̅₁(t, q, v)
    D²ϑ = zeros(eltype(q), length(q), length(v))
    D²ϑd₁(t, q, D²ϑ)
    return transpose(q) * (D²ϑ * v)
end

function g̅₂(t, q, v)
    D²ϑ = zeros(eltype(q), length(q), length(v))
    D²ϑd₂(t, q, D²ϑ)
    return transpose(q) * (D²ϑ * v)
end

function g̅₃(t, q, v)
    D²ϑ = zeros(eltype(q), length(q), length(v))
    D²ϑd₃(t, q, D²ϑ)
    return transpose(q) * (D²ϑ * v)
end

function g̅₄(t, q, v)
    D²ϑ = zeros(eltype(q), length(q), length(v))
    D²ϑd₄(t, q, D²ϑ)
    return transpose(q) * (D²ϑ * v)
end


function guiding_center_4d_v(v::AbstractVector{DT}, t, q::AbstractVector{DT}, params) where {DT}
    @unpack μ = params

    local lB₁ = B₁(t,q)
    local lB₂ = B₂(t,q)
    local lB₃ = B₃(t,q)

    local lβ₁ = β₁(t,q)
    local lβ₂ = β₂(t,q)
    local lβ₃ = β₃(t,q)
    local lβ  = lβ₁ * dϑ₁dx₄(t,q) + lβ₂ * dϑ₂dx₄(t,q) + lβ₃ * dϑ₃dx₄(t,q)

    local ∇₁B = dBdx₁(t,q)
    local ∇₂B = dBdx₂(t,q)
    local ∇₃B = dBdx₃(t,q)

    v[1] = ( u(t,q) * lβ₁ - μ * ( ∇₂B * dϑ₃dx₄(t,q) - ∇₃B * dϑ₂dx₄(t,q) ) ) / lβ
    v[2] = ( u(t,q) * lβ₂ - μ * ( ∇₃B * dϑ₁dx₄(t,q) - ∇₁B * dϑ₃dx₄(t,q) ) ) / lβ
    v[3] = ( u(t,q) * lβ₃ - μ * ( ∇₁B * dϑ₂dx₄(t,q) - ∇₂B * dϑ₁dx₄(t,q) ) ) / lβ
    v[4] = - μ * ( ∇₁B * lβ₁ + ∇₂B * lβ₂ + ∇₃B * lβ₃ ) / lβ

    nothing
end

function guiding_center_4d_v(v::AbstractVector{DT}, t, q::AbstractVector{DT}, p::AbstractVector{DT}, params) where {DT}
    guiding_center_4d_v(v, t, q, params)
end

# function guiding_center_4d_ϑ(θ::AbstractVector{DT}, t, q::AbstractVector{DT}) where {DT}
#     ϑ(t, q, Θ)
# end

function guiding_center_4d_ϑ(θ::AbstractVector{DT}, t, q::AbstractVector{DT}, v::AbstractVector{DT}, params) where {DT}
    ϑ(θ, t, q)
end

function guiding_center_4d_ϑ(θ::AbstractVector{DT}, t, q::AbstractVector{DT}, v::AbstractVector{DT}, params, κ) where {DT}
    θ[1] = (1-κ) * ϑ₁(t,q) - κ * f₁(t,q,q)
    θ[2] = (1-κ) * ϑ₂(t,q) - κ * f₂(t,q,q)
    θ[3] = (1-κ) * ϑ₃(t,q) - κ * f₃(t,q,q)
    θ[4] = (1-κ) * ϑ₄(t,q) - κ * f₄(t,q,q)
    nothing
end

function guiding_center_4d_f(f::AbstractVector{DT}, t, q::AbstractVector{DT}, v::AbstractVector{DT}, params) where {DT}
    @unpack μ = params
    f[1] = f₁(t,q,v) - dHdx₁(t,q,μ)
    f[2] = f₂(t,q,v) - dHdx₂(t,q,μ)
    f[3] = f₃(t,q,v) - dHdx₃(t,q,μ)
    f[4] = f₄(t,q,v) - dHdx₄(t,q,μ)
    nothing
end

function guiding_center_4d_f(f::AbstractVector{DT}, t, q::AbstractVector{DT}, v::AbstractVector{DT}, params, κ) where {DT}
    @unpack μ = params
    f[1] = (1-κ) * f₁(t,q,v) - κ * (g₁(t,q,v) + g̅₁(t,q,v)) - dHdx₁(t,q,μ)
    f[2] = (1-κ) * f₂(t,q,v) - κ * (g₂(t,q,v) + g̅₂(t,q,v)) - dHdx₂(t,q,μ)
    f[3] = (1-κ) * f₃(t,q,v) - κ * (g₃(t,q,v) + g̅₃(t,q,v)) - dHdx₃(t,q,μ)
    f[4] = (1-κ) * f₄(t,q,v) - κ * (g₄(t,q,v) + g̅₄(t,q,v)) - dHdx₄(t,q,μ)
    nothing
end

function guiding_center_4d_g(g::AbstractVector{DT}, t::Number, q::AbstractVector{DT}, v::AbstractVector{DT}, λ::AbstractVector{DT}, params) where {DT}
    g[1] = f₁(t,q,λ)
    g[2] = f₂(t,q,λ)
    g[3] = f₃(t,q,λ)
    g[4] = f₄(t,q,λ)
    nothing
end

# function guiding_center_4d_g(g::AbstractVector{DT}, t, q::AbstractVector{DT}, λ::AbstractVector{DT}, params, κ::DT=0) where {DT}
#     g[1] = g₁(t,q,λ)
#     g[2] = g₂(t,q,λ)
#     g[3] = g₃(t,q,λ)
#     g[4] = g₄(t,q,λ)
#     nothing
# end

function guiding_center_4d_g(g::AbstractVector{DT}, t, q::AbstractVector{DT}, λ::AbstractVector{DT}, params, κ::DT) where {DT}
    g[1] = (1-κ) * f₁(t,q,λ) - κ * (g₁(t,q,λ) + g̅₁(t,q,λ))
    g[2] = (1-κ) * f₂(t,q,λ) - κ * (g₂(t,q,λ) + g̅₂(t,q,λ))
    g[3] = (1-κ) * f₃(t,q,λ) - κ * (g₃(t,q,λ) + g̅₃(t,q,λ))
    g[4] = (1-κ) * f₄(t,q,λ) - κ * (g₄(t,q,λ) + g̅₄(t,q,λ))
    nothing
end

# function guiding_center_4d_g(g::AbstractVector{DT}, t, q::AbstractVector{DT}, λ::AbstractVector{DT}, params, κ::DT=0) where {DT}
#     g[1] = (1-κ) * g₁(t,q,λ) - κ * (f₁(t,q,λ) + g̅₁(t,q,λ))
#     g[2] = (1-κ) * g₂(t,q,λ) - κ * (f₂(t,q,λ) + g̅₂(t,q,λ))
#     g[3] = (1-κ) * g₃(t,q,λ) - κ * (f₃(t,q,λ) + g̅₃(t,q,λ))
#     g[4] = (1-κ) * g₄(t,q,λ) - κ * (f₄(t,q,λ) + g̅₄(t,q,λ))
#     nothing
# end

function guiding_center_4d_dH(∇H::AbstractVector{DT}, t::Number, q::AbstractVector{DT}, params) where {DT}
    dH(∇H, t, q, params)
end

function guiding_center_4d_ω(Ω::Matrix{DT}, t::Number, q::AbstractVector{DT}, v::AbstractVector{DT}, params) where {DT}
    ω(Ω, t, q)
end
    

function guiding_center_4d_λ(λ::AbstractVector{DT}, t::Number, q::AbstractVector{DT}, μ, Ω::AbstractMatrix{DT}, dh::AbstractVector{DT}) where {DT}
    dH(dh, t, q, μ)
    ω(Ω, t, q)
    λ .= inv(Ω) * dh
    nothing
end

function guiding_center_4d_λ(λ::AbstractVector{DT}, t::Number, q::AbstractVector{DT}, μ) where {DT}
    D = length(q)
    guiding_center_4d_λ(λ, t, q, μ, zeros(DT,D,D), zeros(DT,D))
end


function guiding_center_4d_pᵢ(tᵢ, qᵢ::AbstractArray{T}) where {T <: Number}
    pᵢ = zero(qᵢ)
    ϑ(pᵢ, tᵢ, qᵢ)
    return pᵢ
end

function guiding_center_4d_pᵢ(tᵢ, qᵢ::AbstractVector{<:AbstractArray{T}}) where {T <: Number}
    [guiding_center_4d_pᵢ(tᵢ, q) for q in qᵢ]
end


function guiding_center_4d_λᵢ(qᵢ::AbstractVector{DT}, μ, Δt::Number=1., tᵢ::Number=0.) where {DT}
    D  = length(qᵢ)
    λᵢ = zero(qᵢ)
    Ω  = zeros(DT, D, D)
    dh = zeros(DT, D)
    guiding_center_4d_λ(λᵢ, qᵢ, μ, Ω, dh, Δt, tᵢ)
    return λᵢ
end
