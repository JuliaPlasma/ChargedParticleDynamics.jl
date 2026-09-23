@doc raw"""
Charged particle in a uniform magnetic field of the form ``B(x,y,z) = B_{0} e_{z}``, in the
noncanonical formulation on the phasespace ``z = (x, v)``.

The same equilibrium as `ThetaPinchCanonical`, in cartesian coordinates, so the metric is trivial
and the Boris splitting of `sodeproblem` is a valid splitting of the model here.
"""
module ThetaPinchNoncanonical

using ElectromagneticFields: FieldFunctions, ThetaPinchEquilibrium
using ..Noncanonical: ϑ, ϑ₃, fieldpoint
using ...ChargedParticleDynamics: check_chart

export odeproblem, sodeproblem, iodeproblem,
       hamiltonian, toroidal_momentum, ϑ

const qᵢ = [2.5, 0.0, 0.0, 0.0, 0.2, 0.1]

const FIELD = FieldFunctions(ThetaPinchEquilibrium())

function toroidal_momentum(t, q, params)
    check_chart(params.field, FIELD)
    ϑ₃(t, fieldpoint(params.field, t, q))
end

include("charged_particle_3d_noncanonical_presets.jl")

export default_parameters

"""
The charged particle has no physical parameters, so the only entry is the field its equations
read. The method exists so that every problem in this package can be constructed the same way.
"""
default_parameters(::Type{T} = Float64) where {T} = (field = FIELD,)

end
