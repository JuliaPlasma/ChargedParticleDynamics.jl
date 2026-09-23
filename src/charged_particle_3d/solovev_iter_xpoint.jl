@doc raw"""
Analytic ITER-like Solov'ev equilibrium with an X-point, in cylindrical coordinates
``(R, Z, \varphi)``.

The Solov'ev solution of the Grad-Shafranov equation fitted to ITER's shape, with the
separatrix X-point included. See `SolovevIter` for the variant without one.

The hard-coded `μ` and `u` of the `initial_conditions_*` below are inherited from earlier work on
this package and their provenance is not recorded; they are known-good starting points for the
charged particle model rather than values derived here.
"""
module SolovevIterXpoint

using ElectromagneticFields: FieldFunctions, SolovevXpointEquilibriumITER, from_cartesian,
                             b♯, B
using ..Canonical: tᵢ, charged_particle_3d_pᵢ, toroidal_momentum

export podeproblem, iodeproblem,
       hamiltonian, toroidal_momentum

const FIELD = FieldFunctions(SolovevXpointEquilibriumITER())

include("charged_particle_3d_canonical_presets.jl")

export default_parameters

"""
The charged particle has no physical parameters, so the only entry is the field its equations
read. The method exists so that every problem in this package can be constructed the same way.
"""
default_parameters(::Type{T} = Float64) where {T} = (field = FIELD,)

const xᵢ = [7.0 - 1.4, 0.0, 0.0]
const qᵢ = Vector(from_cartesian(FIELD, 0, xᵢ))
const vᵢ = [3.43E-3, 6.75, -3.41E-1]
const pᵢ = charged_particle_3d_pᵢ(tᵢ, qᵢ, vᵢ, default_parameters())

function initial_conditions(x₀, u₀, μ)
    x₀ = Vector(x₀)
    b = b♯(FIELD, 0, x₀)
    vpar = u₀ * b
    vper = sqrt(2 * μ * B(FIELD, 0, x₀))
    v¹ = vper * sqrt(b[3]^2 / (b[1]^2 + b[3]^2))
    v³ = - v¹ * b[1] / b[3]
    v₀ = vpar .+ [v¹, 0, v³]

    (q = x₀, p = charged_particle_3d_pᵢ(tᵢ, x₀, v₀, default_parameters()),
        params = default_parameters())
end

const x₀ = from_cartesian(FIELD, 0, [2.5, 0.0, 0.0])

initial_conditions_barely_passing() = initial_conditions(x₀, 3.425E-1, 1E-2)
initial_conditions_barely_trapped() = initial_conditions(x₀, 3.375E-1, 1E-2)
initial_conditions_deeply_passing() = initial_conditions(x₀, 5E-1, 1E-2)
initial_conditions_deeply_trapped() = initial_conditions(x₀, 1E-1, 1E-2)
initial_conditions_trapped() = initial_conditions(x₀, -2E-3, 1.88E-7)

end
