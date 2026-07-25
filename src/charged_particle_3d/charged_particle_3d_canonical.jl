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

# Derivatives of the one-form with respect to the position coordinates. The metric depends on
# position, hence the gᵢᵢ derivative terms; they are what makes these differ from ∂Aᵢ/∂xⱼ in every
# curvilinear coordinate system.
dϑ₁dx₁(t, q, v) = dg₁₁dx₁(t, q) * v[1] + dA₁dx₁(t, q)
dϑ₁dx₂(t, q, v) = dg₁₁dx₂(t, q) * v[1] + dA₁dx₂(t, q)
dϑ₁dx₃(t, q, v) = dg₁₁dx₃(t, q) * v[1] + dA₁dx₃(t, q)

dϑ₂dx₁(t, q, v) = dg₂₂dx₁(t, q) * v[2] + dA₂dx₁(t, q)
dϑ₂dx₂(t, q, v) = dg₂₂dx₂(t, q) * v[2] + dA₂dx₂(t, q)
dϑ₂dx₃(t, q, v) = dg₂₂dx₃(t, q) * v[2] + dA₂dx₃(t, q)

dϑ₃dx₁(t, q, v) = dg₃₃dx₁(t, q) * v[3] + dA₃dx₁(t, q)
dϑ₃dx₂(t, q, v) = dg₃₃dx₂(t, q) * v[3] + dA₃dx₂(t, q)
dϑ₃dx₃(t, q, v) = dg₃₃dx₃(t, q) * v[3] + dA₃dx₃(t, q)

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

# The projection force is gᵢ = (∂ϑⱼ/∂qⁱ) λʲ, which is *linear* in λ. It was previously written with
# λʲ² in the metric term, which is neither linear in λ nor the derivative of the one-form; `v` is
# available in the IODE `g` signature and belongs in that term. Written through the `dϑⱼdxᵢ` above
# so that it cannot drift away from the one-form it is meant to be the gradient of.
function charged_particle_3d_iode_g(g, t, q, v, λ, params)
    g[1] = dϑ₁dx₁(t, q, v) * λ[1] + dϑ₂dx₁(t, q, v) * λ[2] + dϑ₃dx₁(t, q, v) * λ[3]
    g[2] = dϑ₁dx₂(t, q, v) * λ[1] + dϑ₂dx₂(t, q, v) * λ[2] + dϑ₃dx₂(t, q, v) * λ[3]
    g[3] = dϑ₁dx₃(t, q, v) * λ[1] + dϑ₂dx₃(t, q, v) * λ[2] + dϑ₃dx₃(t, q, v) * λ[3]
    nothing
end


@doc raw"""
The symplectic two-form of the canonical formulation,
``\Omega_{ij} = \partial \vartheta_{i} / \partial q^{j} - \partial \vartheta_{j} / \partial q^{i}``
with ``\vartheta_{i} = g_{ii} v^{i} + A_{i}``, in the sign convention used throughout this package.

This is the ``D \times D`` form on the *configuration* space that the `LODE` interface asks for,
with ``D = 3`` the length of `q` — not the ``6 \times 6`` canonical form on ``(q,p)``. The two are
different objects, and `GeometricIntegrators` sizes the buffer it hands in from the problem
dimension, so the canonical one does not fit.

It depends on `v` through the metric terms and reduces to the magnetic field, ``\Omega_{ij} =
\partial_{j} A_{i} - \partial_{i} A_{j}``, in cartesian coordinates.

Note that no `GeometricIntegrators` integrator currently evaluates the two-form of an `LODE`; it
is carried by the problem for completeness and for use by projection methods.
"""
function ω(Ω, t, q, v, params)
    Ω .= 0

    Ω[1, 2] = dϑ₁dx₂(t, q, v) - dϑ₂dx₁(t, q, v)
    Ω[1, 3] = dϑ₁dx₃(t, q, v) - dϑ₃dx₁(t, q, v)
    Ω[2, 3] = dϑ₂dx₃(t, q, v) - dϑ₃dx₂(t, q, v)

    Ω[2, 1] = -Ω[1, 2]
    Ω[3, 1] = -Ω[1, 3]
    Ω[3, 2] = -Ω[2, 3]

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
