import GeometricEquations
using GeometricEquations: IODEProblem, LODEProblem, PODEProblem
using GeometricSolutions: GeometricSolution, DataSeries, ScalarDataSeries, TimeSeries
using GeometricSolutions: compute_invariant, compute_invariant_error


const Δt = 0.01
const tᵢ = 0.0
const tspan = (tᵢ, 10.0)

ϑ₁(t, q, v) = g₁₁(t, q) * v[1] + A₁(t, q)
ϑ₂(t, q, v) = g₂₂(t, q) * v[2] + A₂(t, q)
ϑ₃(t, q, v) = g₃₃(t, q) * v[3] + A₃(t, q)

function ϑ(θ, t, q, v)
    θ[1] = ϑ₁(t, q, v)
    θ[2] = ϑ₂(t, q, v)
    θ[3] = ϑ₃(t, q, v)
    nothing
end

v¹(t, q, p) = g¹¹(t, q) * (p[1] - A₁(t, q))
v²(t, q, p) = g²²(t, q) * (p[2] - A₂(t, q))
v³(t, q, p) = g³³(t, q) * (p[3] - A₃(t, q))

# L = ½ gᵢⱼ vⁱ vʲ + A·v - φ. The A·v term is what makes ϑ = ∂L/∂v the one-form above; without it
# the Lagrangian is not the Legendre dual of `hamiltonian` and `charged_particle_3d_lode` describes
# a different system than `charged_particle_3d_pode`.
lagrangian(t, q, v) = (g₁₁(t, q) * v[1]^2 + g₂₂(t, q) * v[2]^2 + g₃₃(t, q) * v[3]^2) / 2 +
                      A₁(t, q) * v[1] + A₂(t, q) * v[2] + A₃(t, q) * v[3] - φ(t, q)
hamiltonian(t, q, p) = (g₁₁(t, q) * v¹(t, q, p)^2 + g₂₂(t, q) * v²(t, q, p)^2 + g₃₃(t, q) * v³(t, q, p)^2) / 2 + φ(t, q)
toroidal_momentum(t, q, p) = p[3]

lagrangian(t, q, v, params) = lagrangian(t, q, v)
hamiltonian(t, q, p, params) = hamiltonian(t, q, p)


function charged_particle_3d_pᵢ(tᵢ, qᵢ, vᵢ)
    pᵢ = zero(qᵢ)
    ϑ(pᵢ, tᵢ, qᵢ, vᵢ)
    pᵢ
end


function charged_particle_3d_pode_v(v, t, q, p, params)
    v[1] = v¹(t, q, p)
    v[2] = v²(t, q, p)
    v[3] = v³(t, q, p)
    nothing
end

function charged_particle_3d_pode_f(f, t, q, p, params)
    f[1] = dA₁dx₁(t, q) * v¹(t, q, p) + dA₂dx₁(t, q) * v²(t, q, p) + dA₃dx₁(t, q) * v³(t, q, p) + E₁(t, q) +
           (dg₁₁dx₁(t, q) * v¹(t, q, p)^2 + dg₂₂dx₁(t, q) * v²(t, q, p)^2 + dg₃₃dx₁(t, q) * v³(t, q, p)^2) / 2
    f[2] = dA₁dx₂(t, q) * v¹(t, q, p) + dA₂dx₂(t, q) * v²(t, q, p) + dA₃dx₂(t, q) * v³(t, q, p) + E₂(t, q) +
           (dg₁₁dx₂(t, q) * v¹(t, q, p)^2 + dg₂₂dx₂(t, q) * v²(t, q, p)^2 + dg₃₃dx₂(t, q) * v³(t, q, p)^2) / 2
    f[3] = dA₁dx₃(t, q) * v¹(t, q, p) + dA₂dx₃(t, q) * v²(t, q, p) + dA₃dx₃(t, q) * v³(t, q, p) + E₃(t, q) +
           (dg₁₁dx₃(t, q) * v¹(t, q, p)^2 + dg₂₂dx₃(t, q) * v²(t, q, p)^2 + dg₃₃dx₃(t, q) * v³(t, q, p)^2) / 2
    nothing
end

charged_particle_3d_iode_ϑ(θ, t, q, v, params) = ϑ(θ, t, q, v)

function charged_particle_3d_iode_f(f, t, q, v, params)
    f[1] = dA₁dx₁(t, q) * v[1] + dA₂dx₁(t, q) * v[2] + dA₃dx₁(t, q) * v[3] + E₁(t, q) +
           (dg₁₁dx₁(t, q) * v[1]^2 + dg₂₂dx₁(t, q) * v[2]^2 + dg₃₃dx₁(t, q) * v[3]^2) / 2
    f[2] = dA₁dx₂(t, q) * v[1] + dA₂dx₂(t, q) * v[2] + dA₃dx₂(t, q) * v[3] + E₂(t, q) +
           (dg₁₁dx₂(t, q) * v[1]^2 + dg₂₂dx₂(t, q) * v[2]^2 + dg₃₃dx₂(t, q) * v[3]^2) / 2
    f[3] = dA₁dx₃(t, q) * v[1] + dA₂dx₃(t, q) * v[2] + dA₃dx₃(t, q) * v[3] + E₃(t, q) +
           (dg₁₁dx₃(t, q) * v[1]^2 + dg₂₂dx₃(t, q) * v[2]^2 + dg₃₃dx₃(t, q) * v[3]^2) / 2
    nothing
end

# The projection force is gᵢ = (∂ϑⱼ/∂qⁱ) λʲ = (∂Aⱼ/∂qⁱ + ∂gⱼⱼ/∂qⁱ vʲ) λʲ, which is *linear* in λ.
# It was previously written with λʲ² in the metric term, which is neither linear in λ nor the
# derivative of the one-form; `v` is available in the IODE `g` signature and belongs in that term.
function charged_particle_3d_iode_g(g, t, q, v, λ, params)
    g[1] = dA₁dx₁(t, q) * λ[1] + dA₂dx₁(t, q) * λ[2] + dA₃dx₁(t, q) * λ[3] +
           dg₁₁dx₁(t, q) * v[1] * λ[1] + dg₂₂dx₁(t, q) * v[2] * λ[2] + dg₃₃dx₁(t, q) * v[3] * λ[3]
    g[2] = dA₁dx₂(t, q) * λ[1] + dA₂dx₂(t, q) * λ[2] + dA₃dx₂(t, q) * λ[3] +
           dg₁₁dx₂(t, q) * v[1] * λ[1] + dg₂₂dx₂(t, q) * v[2] * λ[2] + dg₃₃dx₂(t, q) * v[3] * λ[3]
    g[3] = dA₁dx₃(t, q) * λ[1] + dA₂dx₃(t, q) * λ[2] + dA₃dx₃(t, q) * λ[3] +
           dg₁₁dx₃(t, q) * v[1] * λ[1] + dg₂₂dx₃(t, q) * v[2] * λ[2] + dg₃₃dx₃(t, q) * v[3] * λ[3]
    nothing
end


@doc raw"""
The symplectic two-form of the canonical formulation. The Lagrangian is regular, so the phase
space is ``(q,p) \in \mathbb{R}^{6}`` and the form is the canonical one,
``\Omega = \begin{pmatrix} 0 & \mathbb{1} \\ - \mathbb{1} & 0 \end{pmatrix}``,
in the sign convention ``\Omega_{ij} = \partial \vartheta_{i} / \partial z^{j} -
\partial \vartheta_{j} / \partial z^{i}`` used throughout this package.

Note that no `GeometricIntegrators` integrator currently evaluates the two-form of an `LODE`; it
is carried by the problem for completeness and for use by projection methods.
"""
function ω(Ω, t, q, params)
    Ω .= 0
    for i in 1:3
        Ω[i, 3+i] = +1
        Ω[3+i, i] = -1
    end
    nothing
end

ω(Ω, t, q, v, params) = ω(Ω, t, q, params)


# function charged_particle_3d_pode(q₀=qᵢ, v₀=vᵢ)
#     PODE(charged_particle_3d_pode_v, charged_particle_3d_pode_f, q₀, charged_particle_3d_pᵢ(q₀, v₀))
# end
function charged_particle_3d_pode(q₀=qᵢ, p₀=pᵢ; tspan=tspan, tstep=Δt)
    PODEProblem(
        charged_particle_3d_pode_v,
        charged_particle_3d_pode_f,
        tspan, tstep, q₀, p₀;
        invariants=(h=hamiltonian,))
end


function charged_particle_3d_iode(q₀=qᵢ, p₀=pᵢ; tspan=tspan, tstep=Δt)
    IODEProblem(
        charged_particle_3d_iode_ϑ,
        charged_particle_3d_iode_f,
        charged_particle_3d_iode_g,
        tspan, tstep, q₀, p₀;
        invariants=(h=hamiltonian,),
        v̄=charged_particle_3d_pode_v
    )
end

function charged_particle_3d_lode(q₀=qᵢ, p₀=pᵢ; tspan=tspan, tstep=Δt)
    LODEProblem(
        charged_particle_3d_iode_ϑ,
        charged_particle_3d_iode_f,
        charged_particle_3d_iode_g,
        ω, lagrangian,
        tspan, tstep, q₀, p₀;
        invariants=(h=hamiltonian,),
        v̄=charged_particle_3d_pode_v)
end


# `Union{TimeSeries, ScalarDataSeries}`: since GeometricSolutions 0.6 a solution's `sol.t` is a
# `ScalarDataSeries`, not a `TimeSeries`, so annotating only the latter made the `GeometricSolution`
# methods below fail to dispatch.
compute_energy(t::Union{TimeSeries, ScalarDataSeries}, q::DataSeries, p::DataSeries, params) = compute_invariant(t, q, p, params, hamiltonian)
compute_energy(sol::GeometricSolution) = compute_energy(sol.t, sol.q, sol.p, GeometricEquations.parameters(sol.problem))

compute_energy_error(t::Union{TimeSeries, ScalarDataSeries}, q::DataSeries, p::DataSeries, params) = compute_invariant_error(t, q, p, params, hamiltonian)
compute_energy_error(sol::GeometricSolution) = compute_energy_error(sol.t, sol.q, sol.p, GeometricEquations.parameters(sol.problem))
