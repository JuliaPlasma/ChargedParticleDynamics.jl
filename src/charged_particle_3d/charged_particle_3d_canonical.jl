using GeometricEquations: IODEProblem, LODEProblem, PODEProblem
using GeometricSolutions: GeometricSolution, DataSeries, TimeSeries

import GeometricProblems.Diagnostics: compute_invariant, compute_invariant_error


const Δt = 0.01
const tᵢ = 0.0
const tspan = (tᵢ, 10.0)

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

lagrangian(t,q,v) = (g₁₁(t,q) * v[1]^2 + g₂₂(t,q) * v[2]^2 + g₃₃(t,q) * v[3]^2) / 2 - φ(t,q)
hamiltonian(t,q,p) = (g₁₁(t,q) * v¹(t,q,p)^2 + g₂₂(t,q) * v²(t,q,p)^2 + g₃₃(t,q) * v³(t,q,p)^2) / 2 + φ(t,q)
toroidal_momentum(t,q,p) = p[3]


function charged_particle_3d_pᵢ(tᵢ, qᵢ, vᵢ)
    pᵢ = zero(qᵢ)
    ϑ(pᵢ, tᵢ, qᵢ, vᵢ)
    pᵢ
end


function charged_particle_3d_pode_v(v, t, q, p, params)
    v[1] = v¹(t,q,p)
    v[2] = v²(t,q,p)
    v[3] = v³(t,q,p)
    nothing
end

function charged_particle_3d_pode_f(f, t, q, p, params)
    f[1] = dA₁dx₁(t,q) * v¹(t,q,p) + dA₂dx₁(t,q) * v²(t,q,p) + dA₃dx₁(t,q) * v³(t,q,p) + E₁(t,q) +
           (dg₁₁dx₁(t,q) * v¹(t,q,p)^2 + dg₂₂dx₁(t,q) * v²(t,q,p)^2 + dg₃₃dx₁(t,q) * v³(t,q,p)^2) / 2
    f[2] = dA₁dx₂(t,q) * v¹(t,q,p) + dA₂dx₂(t,q) * v²(t,q,p) + dA₃dx₂(t,q) * v³(t,q,p) + E₂(t,q) +
           (dg₁₁dx₂(t,q) * v¹(t,q,p)^2 + dg₂₂dx₂(t,q) * v²(t,q,p)^2 + dg₃₃dx₂(t,q) * v³(t,q,p)^2) / 2
    f[3] = dA₁dx₃(t,q) * v¹(t,q,p) + dA₂dx₃(t,q) * v²(t,q,p) + dA₃dx₃(t,q) * v³(t,q,p) + E₃(t,q) +
           (dg₁₁dx₃(t,q) * v¹(t,q,p)^2 + dg₂₂dx₃(t,q) * v²(t,q,p)^2 + dg₃₃dx₃(t,q) * v³(t,q,p)^2) / 2
    nothing
end

charged_particle_3d_iode_ϑ(θ, t, q, v, params) = ϑ(θ, t, q, v)

function charged_particle_3d_iode_f(f, t, q, v, params)
    f[1] = dA₁dx₁(t,q) * v[1] + dA₂dx₁(t,q) * v[2] + dA₃dx₁(t,q) * v[3] + E₁(t,q) +
           (dg₁₁dx₁(t,q) * v[1]^2 + dg₂₂dx₁(t,q) * v[2]^2 + dg₃₃dx₁(t,q) * v[3]^2) / 2
    f[2] = dA₁dx₂(t,q) * v[1] + dA₂dx₂(t,q) * v[2] + dA₃dx₂(t,q) * v[3] + E₂(t,q) +
           (dg₁₁dx₂(t,q) * v[1]^2 + dg₂₂dx₂(t,q) * v[2]^2 + dg₃₃dx₂(t,q) * v[3]^2) / 2
    f[3] = dA₁dx₃(t,q) * v[1] + dA₂dx₃(t,q) * v[2] + dA₃dx₃(t,q) * v[3] + E₃(t,q) +
           (dg₁₁dx₃(t,q) * v[1]^2 + dg₂₂dx₃(t,q) * v[2]^2 + dg₃₃dx₃(t,q) * v[3]^2) / 2
    nothing
end

function charged_particle_3d_iode_g(g, t, q, v, λ, params)
    g[1] = dA₁dx₁(t,q) * λ[1] + dA₂dx₁(t,q) * λ[2] + dA₃dx₁(t,q) * λ[3] +
           (dg₁₁dx₁(t,q) * λ[1]^2 + dg₂₂dx₁(t,q) * λ[2]^2 + dg₃₃dx₁(t,q) * λ[3]^2) / 2
    g[2] = dA₁dx₂(t,q) * λ[1] + dA₂dx₂(t,q) * λ[2] + dA₃dx₂(t,q) * λ[3] +
           (dg₁₁dx₂(t,q) * λ[1]^2 + dg₂₂dx₂(t,q) * λ[2]^2 + dg₃₃dx₂(t,q) * λ[3]^2) / 2
    g[3] = dA₁dx₃(t,q) * λ[1] + dA₂dx₃(t,q) * λ[2] + dA₃dx₃(t,q) * λ[3] +
           (dg₁₁dx₃(t,q) * λ[1]^2 + dg₂₂dx₃(t,q) * λ[2]^2 + dg₃₃dx₃(t,q) * λ[3]^2) / 2
    nothing
end


# function charged_particle_3d_pode(q₀=qᵢ, v₀=vᵢ)
#     PODE(charged_particle_3d_pode_v, charged_particle_3d_pode_f, q₀, charged_particle_3d_pᵢ(q₀, v₀))
# end
function charged_particle_3d_pode(q₀=qᵢ, p₀=pᵢ; tspan=tspan, tstep=Δt)
    PODEProblem(
        charged_particle_3d_pode_v,
        charged_particle_3d_pode_f,
        tspan, tstep, q₀, p₀;
        invariants = (h = hamiltonian,))
end


function charged_particle_3d_iode(q₀=qᵢ, p₀=pᵢ; tspan=tspan, tstep=Δt)
    IODEProblem(
        charged_particle_3d_iode_ϑ,
        charged_particle_3d_iode_f,
        charged_particle_3d_iode_g,
        tspan, tstep, q₀, p₀;
        invariants = (h = hamiltonian,))
end

function charged_particle_3d_lode(q₀=qᵢ, p₀=pᵢ; tspan=tspan, tstep=Δt)
    LODEProblem(
        charged_particle_3d_iode_ϑ,
        charged_particle_3d_iode_f,
        charged_particle_3d_iode_g,
        ω, lagrangian,
        tspan, tstep, q₀, p₀;
        invariants = (h = hamiltonian,))
end


compute_energy(t::TimeSeries, q::DataSeries, p::DataSeries) = compute_invariant(t, q, p, hamiltonian)
compute_energy(sol::GeometricSolution) = compute_energy(sol.t, sol.q, sol.p)

compute_energy_error(t::TimeSeries, q::DataSeries, p::DataSeries) = compute_invariant_error(t, q, p, hamiltonian)
compute_energy_error(sol::GeometricSolution) = compute_energy_error(sol.t, sol.q, sol.p)
