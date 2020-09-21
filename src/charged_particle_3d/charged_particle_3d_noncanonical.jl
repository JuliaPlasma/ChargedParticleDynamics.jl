
ϑ₁(t, q) = g₁₁(t,q) * q[4] + A₁(t,q)
ϑ₂(t, q) = g₂₂(t,q) * q[5] + A₂(t,q)
ϑ₃(t, q) = g₃₃(t,q) * q[6] + A₃(t,q)


function ϑ(t, q, θ)
    θ[1] = ϑ₁(t,q)
    θ[2] = ϑ₂(t,q)
    θ[3] = ϑ₃(t,q)
    θ[4] = zero(eltype(q))
    θ[5] = zero(eltype(q))
    θ[6] = zero(eltype(q))
    nothing
end


function ω(t, q, Β)
    Β[1,1] = 0
    Β[1,2] = dϑ₁dx₂(t,q) - dϑ₂dx₁(t,q)
    Β[1,3] = dϑ₁dx₃(t,q) - dϑ₃dx₁(t,q)
    Β[1,4] = +1
    Β[1,5] = 0
    Β[1,6] = 0

    Β[2,1] = dϑ₂dx₁(t,q) - dϑ₁dx₂(t,q)
    Β[2,2] = 0
    Β[2,3] = dϑ₂dx₃(t,q) - dϑ₃dx₂(t,q)
    Β[2,4] = 0
    Β[2,5] = +1
    Β[2,6] = 0

    Β[3,1] = dϑ₃dx₁(t,q) - dϑ₁dx₃(t,q)
    Β[3,2] = dϑ₃dx₂(t,q) - dϑ₂dx₃(t,q)
    Β[3,3] = 0
    Β[3,4] = 0
    Β[3,5] = 0
    Β[3,6] = +1

    Β[4,1] = -1
    Β[4,2] = 0
    Β[4,3] = 0
    Β[4,4] = 0
    Β[4,5] = 0
    Β[4,6] = 0

    Β[5,1] = 0
    Β[5,2] = -1
    Β[5,3] = 0
    Β[5,4] = 0
    Β[5,5] = 0
    Β[5,6] = 0

    Β[6,1] = 0
    Β[6,2] = 0
    Β[6,3] = -1
    Β[6,4] = 0
    Β[6,5] = 0
    Β[6,6] = 0

    nothing
end


function dϑ(t, q, dϑ)
    dϑ[1,1] = dϑ₁dx₁(t,q)
    dϑ[1,2] = dϑ₁dx₂(t,q)
    dϑ[1,3] = dϑ₁dx₃(t,q)
    dϑ[1,4] = 0
    dϑ[1,5] = 0
    dϑ[1,6] = 0

    dϑ[2,1] = dϑ₂dx₁(t,q)
    dϑ[2,2] = dϑ₂dx₂(t,q)
    dϑ[2,3] = dϑ₂dx₃(t,q)
    dϑ[2,4] = 0
    dϑ[2,5] = 0
    dϑ[2,6] = 0

    dϑ[3,1] = dϑ₃dx₁(t,q)
    dϑ[3,2] = dϑ₃dx₂(t,q)
    dϑ[3,3] = dϑ₃dx₃(t,q)
    dϑ[3,4] = 0
    dϑ[3,5] = 0
    dϑ[3,6] = 0

    dϑ[4,1] = 0
    dϑ[4,2] = 0
    dϑ[4,3] = 0
    dϑ[4,4] = 0
    dϑ[4,5] = 0
    dϑ[4,6] = 0

    dϑ[5,1] = 0
    dϑ[5,2] = 0
    dϑ[5,3] = 0
    dϑ[5,4] = 0
    dϑ[5,5] = 0
    dϑ[5,6] = 0

    dϑ[6,1] = 0
    dϑ[6,2] = 0
    dϑ[6,3] = 0
    dϑ[6,4] = 0
    dϑ[6,5] = 0
    dϑ[6,6] = 0

    nothing
end


β₁(t,q) = dϑ₃dx₂(t,q) - dϑ₂dx₃(t,q)
β₂(t,q) = dϑ₁dx₃(t,q) - dϑ₃dx₁(t,q)
β₃(t,q) = dϑ₂dx₁(t,q) - dϑ₁dx₂(t,q)

β(t,q) = sqrt(β₁(t,q)^2 + β₂(t,q)^2 + β₃(t,q)^2)


hamiltonian(t,q) = 0.5 * (g₁₁(t,q) * q[4]^2 + g₂₂(t,q) * q[5]^2 + g₃₃(t,q) * q[6]^2)
hamiltonian(t,q,p) = hamiltonian(t,q)

toroidal_momentum(t,q) = ϑ₃(t,q)


dHdx₁(t, q) = zero(eltype(q))
dHdx₂(t, q) = zero(eltype(q))
dHdx₃(t, q) = zero(eltype(q))
dHdx₄(t, q) = q[4]
dHdx₅(t, q) = q[5]
dHdx₆(t, q) = q[6]

function dH(t, q, dH)
    dH[1] = dHdx₁(t, q)
    dH[2] = dHdx₂(t, q)
    dH[3] = dHdx₃(t, q)
    dH[4] = dHdx₄(t, q)
    dH[5] = dHdx₅(t, q)
    dH[6] = dHdx₆(t, q)
    nothing
end


f₁(t, q, v) = q[4]
f₂(t, q, v) = q[5]
f₃(t, q, v) = q[6]
f₄(t, q, v) = q[5] * B₃(t,q) - q[6] * B₂(t,q)
f₅(t, q, v) = q[6] * B₁(t,q) - q[4] * B₃(t,q)
f₆(t, q, v) = q[4] * B₂(t,q) - q[5] * B₁(t,q)


function charged_particle_3d_periodicity(qᵢ)
    periodicity = zeros(eltype(qᵢ), size(qᵢ,1))
    periodicity[3] = 2π
    return periodicity
end


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


function charged_particle_3d_v(t, q, v)
    v[1] = f₁(t, q, v)
    v[2] = f₂(t, q, v)
    v[3] = f₃(t, q, v)
    v[4] = f₄(t, q, v)
    v[5] = f₅(t, q, v)
    v[6] = f₅(t, q, v)

    nothing
end


charged_particle_3d_iode_ϑ(t, q, v, p) = ϑ(t, q, p)


function charged_particle_3d_iode_f(t, q, v, f)
    f[1] = dA₁dx₁(t,q) * v[1] + dA₂dx₁(t,q) * v[2] + dA₃dx₁(t,q) * v[3] + E₁(t,q)
    f[2] = dA₁dx₂(t,q) * v[1] + dA₂dx₂(t,q) * v[2] + dA₃dx₂(t,q) * v[3] + E₂(t,q)
    f[3] = dA₁dx₃(t,q) * v[1] + dA₂dx₃(t,q) * v[2] + dA₃dx₃(t,q) * v[3] + E₃(t,q)
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


charged_particle_3d_v(t, q, p, v) = charged_particle_3d_v(t, q, v)


function charged_particle_3d_sode_fx(t, q::Vector{DT}, f::Vector{DT}, Δt) where {DT}
    x = q[1:3]
    v = q[4:6]

    f .= vcat(x .+ Δt .* v, v)

    nothing
end

function charged_particle_3d_sode_fv(t, q::Vector{DT}, f::Vector{DT}, Δt) where {DT}
    x = q[1:3]
    v = q[4:6]

    local lB₁ = B₁(t-Δt, x) * x[1]
    local lB₂ = B₂(t-Δt, x) * x[1]
    local lB₃ = B₃(t-Δt, x) / x[1]

    local Ω = zeros(DT,3,3)

    Ω[1,2] = - lB₃
    Ω[1,3] = + lB₂
    Ω[2,1] = + lB₃
    Ω[2,3] = + lB₁
    Ω[3,1] = - lB₂
    Ω[3,2] = - lB₁

    C1 = eye(DT,3) .+ (0.5*Δt) .* Ω
    C2 = eye(DT,3) .- (0.5*Δt) .* Ω
    R  = inv(C1) * C2

    f .= vcat(x, R * v)

    nothing
end


function charged_particle_3d_ode(qᵢ=qᵢ; periodic=true)
    if periodic
        ODE(charged_particle_3d_v, qᵢ; periodicity=charged_particle_3d_periodicity(qᵢ))
    else
        ODE(charged_particle_3d_v, qᵢ)
    end
end

function charged_particle_3d_sode(qᵢ=qᵢ; periodic=true)
    if periodic
        SODE((charged_particle_3d_sode_fv, charged_particle_3d_sode_fx), qᵢ; periodicity=charged_particle_3d_periodicity(qᵢ))
    else
        SODE((charged_particle_3d_sode_fv, charged_particle_3d_sode_fx), qᵢ)
    end
end

function charged_particle_3d_iode(qᵢ=qᵢ)
    IODE(charged_particle_3d_iode_ϑ, charged_particle_3d_iode_f,
            charged_particle_3d_iode_g, qᵢ, charged_particle_3d_pᵢ(qᵢ);
            v=charged_particle_3d_v, h=hamiltonian)
end

