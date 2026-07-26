"""
Guiding centre dynamics in the field of a magnetic dipole.
"""
module Dipole3d

import ElectromagneticFields.Dipole

export initial_conditions_dipole

export hamiltonian

Dipole.@code() # inject magnetic field code

const DEFAULT_TIMESTEP = 0.03
const DEFAULT_TIMESPAN = (0.0, 90_000.0)

initial_conditions_dipole() = merge(initial_conditions(0.0, [1.0, 2.0, 1.0, 0.01]), (params=(μ=1E-2,),))

export default_parameters

"The magnetic moment μ this equilibrium is set up for."
default_parameters(::Type{T}=Float64) where {T} = (μ = T(1E-2),)

const parameters = default_parameters()

include("guiding_center_3d_equations.jl")
include("guiding_center_3d_canonical.jl")
include("guiding_center_3d_diagnostics.jl")

end
