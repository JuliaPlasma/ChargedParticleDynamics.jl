"""
Analytic toy problem with quadratic potentials.
"""
module Dipole3d

import ElectromagneticFields.Dipole

export initial_conditions_dipole

export hamiltonian

Dipole.@code() # inject magnetic field code

const Δt = 0.035
const tspan = (0.0, 35.0)
# const tspan = (0.0, 999.99)

initial_conditions_dipole() = (initial_conditions(0.0, [1.0, 2.0, 1.0, 0.01])..., (μ=1E-2,))

include("guiding_center_3d_equations.jl")
include("guiding_center_3d_diagnostics.jl")

end
