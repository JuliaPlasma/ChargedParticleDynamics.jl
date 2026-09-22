"""
Charged Particle in an uniform magnetic field of the form
``B(x,y,z) = B_0 e_z``.
"""
module ThetaPinchCanonical

using ElectromagneticFields: FieldFunctions, ThetaPinchEquilibrium
using ..Canonical: tᵢ, charged_particle_3d_pᵢ

export iodeproblem, hamiltonian, angular_momentum

const FIELD = FieldFunctions(ThetaPinchEquilibrium())

include("charged_particle_3d_canonical_presets.jl")

export default_parameters

"""
The charged particle has no physical parameters, so the only entry is the field its equations
read. The method exists so that every problem in this package can be constructed the same way.
"""
default_parameters(::Type{T} = Float64) where {T} = (field = FIELD,)

const qᵢ = [2.5, 0.0, 0.0]
const vᵢ = [0.0, 0.2, 0.1]
const pᵢ = charged_particle_3d_pᵢ(tᵢ, qᵢ, vᵢ, default_parameters())

end
