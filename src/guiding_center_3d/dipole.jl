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

export default_parameters, default_constraints

"The magnetic moment μ this equilibrium is set up for."
default_parameters(::Type{T}=Float64) where {T} = (μ = T(1E-2),)

"""
The constraint pair [`hodeproblem`](@ref) and its siblings use here by default.

All three pairs are regular at the initial condition, but only this one stays that way. The orbit
takes `b₁` and `b₂` through zero, so `(g³, g¹)` and `(g², g³)` — which divide by them — hold the
constraints to 1E-3 and the energy to 4E-5 over a hundred steps, and lose the orbit entirely past a
few hundred: `|ΔH/H|` reaches 2E3 by `t = 30`, and `hodeproblem_canonical` meets a NaN there.
`(g¹, g²)` divides by `b₃`, which does not vanish along this orbit, and holds both to 2E-8 out to
`t = 300`. This is the one equilibrium in the package where the pair has to be chosen on its
conditioning along the orbit rather than at the initial condition.
"""
default_constraints() = :g12


include("guiding_center_3d_equations.jl")
include("guiding_center_3d_canonical.jl")
include("guiding_center_3d_diagnostics.jl")

end
