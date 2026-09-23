@doc raw"""
Small axisymmetric tokamak equilibrium, in cartesian coordinates ``(x, y, z)``.

Major radius ``R_{0} = 1``, magnetic field on axis ``B_{0} = 1`` and safety factor
``q_{0} = 2``.

The hard-coded `μ` and `u` of the `initial_conditions_*` below are inherited from earlier work on
this package and their provenance is not recorded; they are known-good starting points for the
charged particle model rather than values derived here.
"""
module TokamakSmallCartesian

using ElectromagneticFields: FieldFunctions, AxisymmetricTokamakCartesianEquilibrium
using ..Canonical: tᵢ, charged_particle_3d_pᵢ, toroidal_momentum

export podeproblem, iodeproblem,
       hamiltonian, toroidal_momentum

const FIELD = FieldFunctions(AxisymmetricTokamakCartesianEquilibrium())

include("charged_particle_3d_canonical_presets.jl")

export default_parameters

"""
The charged particle has no physical parameters, so the only entry is the field its equations
read. The method exists so that every problem in this package can be constructed the same way.
"""
default_parameters(::Type{T} = Float64) where {T} = (field = FIELD,)

const qᵢ = [1.05, 0.0, 0.0]
const vᵢ = [2.1E-3, 4.3E-4, 0.0]
const pᵢ = charged_particle_3d_pᵢ(tᵢ, qᵢ, vᵢ, default_parameters())

end
