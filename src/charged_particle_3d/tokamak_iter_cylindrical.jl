@doc raw"""
ITER-size axisymmetric tokamak equilibrium, in cylindrical coordinates ``(R, Z, \varphi)``.

Major radius ``R_{0} = 6.2``, magnetic field on axis ``B_{0} = 5.3`` and safety factor
``q_{0} = \sqrt{2}`` — the ITER parameters of `ElectromagneticFields`.

The hard-coded `μ` and `u` of the `initial_conditions_*` below are inherited from earlier work on
this package and their provenance is not recorded; they are known-good starting points for the
charged particle model rather than values derived here.
"""
module TokamakIterCylindrical

using ElectromagneticFields: FieldFunctions, AxisymmetricTokamakCylindricalITER
using ..Canonical: tᵢ, charged_particle_3d_pᵢ

export podeproblem, iodeproblem,
       hamiltonian, toroidal_momentum

const FIELD = FieldFunctions(AxisymmetricTokamakCylindricalITER())

include("charged_particle_3d_canonical_presets.jl")

export default_parameters

"""
The charged particle has no physical parameters, so the only entry is the field its equations
read. The method exists so that every problem in this package can be constructed the same way.
"""
default_parameters(::Type{T} = Float64) where {T} = (field = FIELD,)

const qᵢ = [7.0, 0.0, 0.0]
const vᵢ = [3.43E-3, 6.75, -3.41E-1]
const pᵢ = charged_particle_3d_pᵢ(tᵢ, qᵢ, vᵢ, default_parameters())

end
