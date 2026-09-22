"""
Charged Particle in a singular magnetic field of the form
``B(x,y,z) = (x^2 + y^2)^{-3/2} e_z``.
"""
module SingularFieldCanonical

using ElectromagneticFields: FieldFunctions, SingularEquilibrium
using ..Canonical: tᵢ, charged_particle_3d_pᵢ

export odeproblem, iodeproblem
export hamiltonian#, angular_momentum
export compute_energy, compute_energy_error

const FIELD = FieldFunctions(SingularEquilibrium())

include("charged_particle_3d_canonical_presets.jl")

export default_parameters

"""
The charged particle has no physical parameters, so the only entry is the field its equations
read. The method exists so that every problem in this package can be constructed the same way.
"""
default_parameters(::Type{T} = Float64) where {T} = (field = FIELD,)

const qᵢ = [1.0, 0.0, 0.0]
const vᵢ = [0.0, -1.0, 0.0]
const pᵢ = charged_particle_3d_pᵢ(tᵢ, qᵢ, vᵢ, default_parameters())

# angular_momentum(t,q) = q[1] * ϑ₂(t,q) - q[2] * ϑ₁(t,q)

end
