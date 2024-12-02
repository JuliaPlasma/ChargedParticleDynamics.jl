using GeometricEquations: ODEProblem, IODEProblem, LODEProblem, SODEProblem
using GeometricSolutions: GeometricSolution, DataSeries, TimeSeries

import GeometricProblems.Diagnostics: compute_invariant, compute_invariant_error


const Δt = 0.01
const tspan = (0.0, 10.0)

ϑ₁(t, q) = g₁₁(t,q) * q[4] + A₁(t,q)
ϑ₂(t, q) = g₂₂(t,q) * q[5] + A₂(t,q)
ϑ₃(t, q) = g₃₃(t,q) * q[6] + A₃(t,q)


function ϑ(θ, t, q)
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


hamiltonian(t,q) = 0.5 * ( g₁₁(t,q) * q[4]^2 + g₂₂(t,q) * q[5]^2 + g₃₃(t,q) * q[6]^2 ) + φ(t,q)
hamiltonian(t,q,p,params) = hamiltonian(t,q)
lagrangian(t,q,v,params) = ϑ₁(t, q) * v[1] + ϑ₂(t, q) * v[2] + ϑ₃(t, q) * v[3] - hamiltonian(t,q)


dHdx₁(t, q) = dφdx₁(t,q)
dHdx₂(t, q) = dφdx₂(t,q)
dHdx₃(t, q) = dφdx₃(t,q)
dHdx₄(t, q) = q[4]
dHdx₅(t, q) = q[5]
dHdx₆(t, q) = q[6]

function dH(dH, t, q)
    dH[1] = dHdx₁(t, q)
    dH[2] = dHdx₂(t, q)
    dH[3] = dHdx₃(t, q)
    dH[4] = dHdx₄(t, q)
    dH[5] = dHdx₅(t, q)
    dH[6] = dHdx₆(t, q)
    nothing
end


v₁(t, q, v) = q[4]
v₂(t, q, v) = q[5]
v₃(t, q, v) = q[6]
v₄(t, q, v) = E₁(t,q) + q[5] * B₃(t,q) - q[6] * B₂(t,q)
v₅(t, q, v) = E₂(t,q) + q[6] * B₁(t,q) - q[4] * B₃(t,q)
v₆(t, q, v) = E₃(t,q) + q[4] * B₂(t,q) - q[5] * B₁(t,q)


function charged_particle_3d_periodicity(qᵢ)
    period = zeros(eltype(qᵢ), size(qᵢ,1))
    period[3] = 2π
    return period
end


function charged_particle_3d_pᵢ(tᵢ, qᵢ)
    pᵢ = zero(qᵢ)
    ϑ(pᵢ, tᵢ, qᵢ)
    return pᵢ
end


function charged_particle_3d_v(v, t, q, params)
    v[1] = v₁(t, q, v)
    v[2] = v₂(t, q, v)
    v[3] = v₃(t, q, v)
    v[4] = v₄(t, q, v)
    v[5] = v₅(t, q, v)
    v[6] = v₆(t, q, v)

    nothing
end

charged_particle_3d_v(v, t, q, p, params) = charged_particle_3d_v(v, t, q, params)


charged_particle_3d_iode_ϑ(θ, t, q, v, params) = ϑ(θ, t, q)


function charged_particle_3d_iode_f(f, t, q, v, params)
    f[1] = dA₁dx₁(t,q) * v[1] + dA₂dx₁(t,q) * v[2] + dA₃dx₁(t,q) * v[3] + E₁(t,q)
    f[2] = dA₁dx₂(t,q) * v[1] + dA₂dx₂(t,q) * v[2] + dA₃dx₂(t,q) * v[3] + E₂(t,q)
    f[3] = dA₁dx₃(t,q) * v[1] + dA₂dx₃(t,q) * v[2] + dA₃dx₃(t,q) * v[3] + E₃(t,q)
    f[4] = v[1] - q[4]
    f[5] = v[2] - q[5]
    f[6] = v[3] - q[6]
    nothing
end

function charged_particle_3d_iode_g(g, t, q, v, λ, params)
    g[1] = dA₁dx₁(t,q) * λ[1] + dA₂dx₁(t,q) * λ[2] + dA₃dx₁(t,q) * λ[3]
    g[2] = dA₁dx₂(t,q) * λ[1] + dA₂dx₂(t,q) * λ[2] + dA₃dx₂(t,q) * λ[3]
    g[3] = dA₁dx₃(t,q) * λ[1] + dA₂dx₃(t,q) * λ[2] + dA₃dx₃(t,q) * λ[3]
    g[4] = λ[1]
    g[5] = λ[2]
    g[6] = λ[3]
    nothing
end


function charged_particle_3d_sode_fx(q₁::AbstractArray{DT}, t₁, q₀::AbstractArray{DT}, t₀) where {DT}
    @assert axes(q₁) == axes(q₀)

    x = q₀[1:3]
    v = q₀[4:6]

    q₁[1:3] .= x .+ (t₁ - t₀) .* v
    q₁[4:6] .= v

    nothing
end

function charged_particle_3d_sode_fv(q₁::AbstractArray{DT}, t₁, q₀::AbstractArray{DT}, t₀) where {DT}
    @assert axes(q₁) == axes(q₀)
    
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

    C1 = eye(DT,3) .+ (t₁ - t₀) ./ 2 .* Ω
    C2 = eye(DT,3) .- (t₁ - t₀) ./ 2 .* Ω
    R  = inv(C1) * C2

    q₁[1:3] .= x
    q₁[4:6] .= R * v

    nothing
end


function charged_particle_3d_ode(qᵢ=qᵢ; tspan=tspan, tstep=Δt, periodic=true)
    if periodic
        ODEProblem(charged_particle_3d_v, tspan, tstep, qᵢ; invariants=(h=hamiltonian,), periodicity=charged_particle_3d_periodicity(qᵢ))
    else
        ODEProblem(charged_particle_3d_v, tspan, tstep, qᵢ; invariants=(h=hamiltonian,))
    end
end

function charged_particle_3d_sode(qᵢ=qᵢ; tspan=tspan, tstep=Δt, periodic=true)
    if periodic
        SODEProblem(nothing, (charged_particle_3d_sode_fv, charged_particle_3d_sode_fx),
                    tspan, tstep, x₀, periodicity=charged_particle_3d_periodicity(qᵢ))
    else
        SODEProblem(nothing, (charged_particle_3d_sode_fv, charged_particle_3d_sode_fx),
                    tspan, tstep, x₀)
    end
end

function charged_particle_3d_iode(qᵢ=qᵢ; tspan=tspan, tstep=Δt)
    IODEProblem(
        charged_particle_3d_iode_ϑ,
        charged_particle_3d_iode_f,
        charged_particle_3d_iode_g,
        tspan, tstep, qᵢ, charged_particle_3d_pᵢ(tspan[begin], qᵢ);
        invariants = (h = hamiltonian,), 
        v̄ = charged_particle_3d_v)
end

function charged_particle_3d_lode(qᵢ=qᵢ; tspan=tspan, tstep=Δt)
    LODEProblem(
        charged_particle_3d_iode_ϑ,
        charged_particle_3d_iode_f,
        charged_particle_3d_iode_g,
        ω, lagrangian,
        tspan, tstep, qᵢ, charged_particle_3d_pᵢ(tspan[begin], qᵢ);
        invariants = (h = hamiltonian,), 
        v̄ = charged_particle_3d_v)
end


compute_energy(t::TimeSeries, q::DataSeries) = compute_invariant(t, q, hamiltonian)
compute_energy(sol::GeometricSolution) = compute_energy(sol.t, sol.q)

compute_energy_error(t::TimeSeries, q::DataSeries) = compute_invariant_error(t, q, hamiltonian)
compute_energy_error(sol::GeometricSolution) = compute_energy_error(sol.t, sol.q)
