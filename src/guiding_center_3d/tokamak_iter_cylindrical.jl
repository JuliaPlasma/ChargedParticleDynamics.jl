"""
Analytic ITER-like Solov'ev equilibrium with X-point.
"""
module TokamakIterCylindrical

import ElectromagneticFields.AxisymmetricTokamakCylindrical

export initial_conditions_barely_passing, initial_conditions_barely_trapped,
    initial_conditions_deeply_passing, initial_conditions_deeply_trapped,
    initial_conditions_trapped

export hamiltonian, toroidal_momentum

AxisymmetricTokamakCylindrical.@code_iter() # inject magnetic field code

const DEFAULT_TIMESTEP = 1.0
const DEFAULT_TIMESPAN = (0.0, 1000.0)

const qᵢ = [7.0 - 1.4, 0.0, 0.0, 2.8166280889939737]

export default_parameters, default_constraints

"The magnetic moment μ this equilibrium is set up for."
default_parameters(::Type{T}=Float64) where {T} = (μ = T(4.607782183567846),)

"""
The constraint pair [`hodeproblem`](@ref) and its siblings use here by default.

`b₁ = b_R` vanishes on the midplane, where every initial condition of
this equilibrium sits, so `(g³, g¹)` is singular there. `(g¹, g²)` divides by `b₃ = b_φ` instead.
"""
default_constraints() = :g12


include("guiding_center_3d_equations.jl")
include("guiding_center_3d_canonical.jl")

# The canonical toroidal momentum is the covariant φ-component of the momentum. Since p = ϑ on the
# constraint manifold this is simply p₃ — the previous R(t,q) * ϑ₃(t,q) both indexed q[4] of a
# three-component state and carried a spurious factor of R that destroys the conservation.
toroidal_momentum(t, q, p) = p[3]

include("guiding_center_3d_diagnostics.jl")

initial_conditions_barely_passing() = merge(initial_conditions(0, [from_cartesian(0, [2.5, 0.0, 0.0])..., 3.425E-1]), (params=(μ=1E-2,),))
initial_conditions_barely_trapped() = merge(initial_conditions(0, [from_cartesian(0, [2.5, 0.0, 0.0])..., 3.375E-1]), (params=(μ=1E-2,),))
initial_conditions_deeply_passing() = merge(initial_conditions(0, [from_cartesian(0, [2.5, 0.0, 0.0])..., 5E-1]), (params=(μ=1E-2,),))
initial_conditions_deeply_trapped() = merge(initial_conditions(0, [from_cartesian(0, [2.5, 0.0, 0.0])..., 1E-1]), (params=(μ=1E-2,),))
initial_conditions_trapped() = merge(initial_conditions(0, [from_cartesian(0, [7.0, 0.0, 0.0])..., -2E-3]), (params=(μ=1.88E-7,),))

end
