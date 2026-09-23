"""
Charged Particle in an axisymmetric magnetic field of the form
``B(x,y,z) = (1 + x^2 + y^2) e_z``.
"""
module SymmetricField

using ElectromagneticFields: FieldFunctions, SymmetricQuadraticEquilibrium
using ..Noncanonical: ϑ₁, ϑ₂, fieldpoint, compute_energy, compute_energy_error

export odeproblem, iodeproblem
export hamiltonian, angular_momentum
export compute_energy, compute_energy_error

const FIELD = FieldFunctions(SymmetricQuadraticEquilibrium())

const qᵢ = [1.0, 0.0, 0.0, 0.0, 1.0, 1.0]

function angular_momentum(t, q, params)
    P = fieldpoint(params.field, t, q)
    q[1] * ϑ₂(t, P) - q[2] * ϑ₁(t, P)
end

include("charged_particle_3d_noncanonical_presets.jl")

export default_parameters

"""
The charged particle has no physical parameters, so the only entry is the field its equations
read. The method exists so that every problem in this package can be constructed the same way.
"""
default_parameters(::Type{T} = Float64) where {T} = (field = FIELD,)

end
