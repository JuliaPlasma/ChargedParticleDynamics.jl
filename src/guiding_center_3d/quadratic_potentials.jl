"""
Analytic toy problem with quadratic potentials.
"""
module QuadraticPotentials3d

import ElectromagneticFields.QuadraticPotentials

export initial_conditions_quadratic

export hamiltonian

QuadraticPotentials.@code() # inject magnetic field code

const DEFAULT_TIMESTEP = 0.5
const DEFAULT_TIMESPAN = (0.0, 2.5E4)

# const DEFAULT_TIMESTEP = 0.35
# const DEFAULT_TIMESPAN = (0.0, 12250.)

initial_conditions_quadratic() = merge(initial_conditions(0.0, [0.3, 0.2, -1.4, 0.3]), (params=(μ=2.5E-3,),))

export default_parameters, default_constraints

"The magnetic moment μ this equilibrium is set up for."
default_parameters(::Type{T}=Float64) where {T} = (μ = T(2.5E-3),)

"""
The constraint pair [`hodeproblem`](@ref) and its siblings use here by default.

All three pairs are well conditioned in this equilibrium, so this keeps the
historical choice.
"""
default_constraints() = :g31


include("guiding_center_3d_equations.jl")
include("guiding_center_3d_canonical.jl")
include("guiding_center_3d_diagnostics.jl")

end
