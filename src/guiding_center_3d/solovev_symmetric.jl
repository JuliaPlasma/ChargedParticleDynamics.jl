"""
Analytic, quadratic Solov'ev equilibrium.
"""
module SolovevSymmetricField

import ElectromagneticFields.SolovevSymmetric

export initial_conditions_barely_passing, initial_conditions_barely_trapped,
    initial_conditions_deeply_passing, initial_conditions_deeply_trapped

export hamiltonian, toroidal_momentum

SolovevSymmetric.@code(2.0, 5.0, 1.0, 1.0) # inject magnetic field code

const DEFAULT_TIMESTEP = 0.1
const DEFAULT_TIMESPAN = (0.0, 1E3)

initial_conditions_barely_passing() = merge(initial_conditions(0, [2.5, 0.0, 0.0, 3.425E-1]), (params=(μ=1E-2,),)) # Δt=2.5, nt=50
initial_conditions_barely_trapped() = merge(initial_conditions(0, [2.5, 0.0, 0.0, 3.375E-1]), (params=(μ=1E-2,),)) # Δt=3.0, nt=100
initial_conditions_deeply_passing() = merge(initial_conditions(0, [2.5, 0.0, 0.0, 5E-1]), (params=(μ=1E-2,),))     # Δt=2.5, nt=25
initial_conditions_deeply_trapped() = merge(initial_conditions(0, [2.5, 0.0, 0.0, 1E-1]), (params=(μ=1E-2,),))     # Δt=5.0, nt=50

export default_parameters

"The magnetic moment μ this equilibrium is set up for."
default_parameters(::Type{T}=Float64) where {T} = (μ = T(1E-2),)

const parameters = default_parameters()

include("guiding_center_3d_equations.jl")
include("guiding_center_3d_canonical.jl")

# The canonical toroidal momentum is the covariant φ-component of the momentum. Since p = ϑ on the
# constraint manifold this is simply p₃ — the previous R(t,q) * ϑ₃(t,q) both indexed q[4] of a
# three-component state and carried a spurious factor of R that destroys the conservation.
toroidal_momentum(t, q, p) = p[3]

include("guiding_center_3d_diagnostics.jl")

end
