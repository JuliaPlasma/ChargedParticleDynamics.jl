"""
Analytic toy problem with quadratic potentials.
"""
module QuadraticPotentials3d

import ElectromagneticFields.QuadraticPotentials

export initial_conditions_quadratic

export hamiltonian

QuadraticPotentials.@code() # inject magnetic field code

const Δt = 0.5
const tspan = (0.0, 2.5E4)

# const Δt = 0.35
# const tspan = (0.0, 12250.)

initial_conditions_quadratic() = merge(initial_conditions(0.0, [0.3, 0.2, -1.4, 0.3]), (params=(μ=2.5E-3,),))

include("guiding_center_3d_equations.jl")
include("guiding_center_3d_canonical.jl")
include("guiding_center_3d_diagnostics.jl")

end
