"""
Analytic, quadratic Solov'ev equilibrium.
"""
module SolovevSymmetricField

import ElectromagneticFields.SolovevSymmetric

export initial_conditions_barely_passing, initial_conditions_barely_trapped,
    initial_conditions_deeply_passing, initial_conditions_deeply_trapped

export hamiltonian, toroidal_momentum

SolovevSymmetric.@code(2.0, 5.0, 1.0, 1.0) # inject magnetic field code

# `hodeproblem_compact` holds this span, and beyond it: over `t ∈ [0, 1E3]` it stays at 1.1E-7.
# `hodeproblem` and `hodeproblem_canonical` do not — they lose the orbit between `t = 10` and
# `t = 50` at *any* step size, reaching 2.4E+6 and 1.5E+228 over that longer span, and a finer step
# does not rescue them (`hodeproblem` is 2.2E+9 at `Δt = 0.05` over `t ∈ [0, 50]`). Over
# `t ∈ [0, 1E1]` both sit at 4.5E-12, so they carry a shorter span where they are integrated. That
# only the compact form survives a long run here is the point worth knowing about this equilibrium.
const DEFAULT_TIMESTEP = 0.1
const DEFAULT_TIMESPAN = (0.0, 1E2)

initial_conditions_barely_passing() = merge(initial_conditions(0, [2.5, 0.0, 0.0, 3.425E-1]), (params=(μ=1E-2,),)) # Δt=2.5, nt=50
initial_conditions_barely_trapped() = merge(initial_conditions(0, [2.5, 0.0, 0.0, 3.375E-1]), (params=(μ=1E-2,),)) # Δt=3.0, nt=100
initial_conditions_deeply_passing() = merge(initial_conditions(0, [2.5, 0.0, 0.0, 5E-1]), (params=(μ=1E-2,),))     # Δt=2.5, nt=25
initial_conditions_deeply_trapped() = merge(initial_conditions(0, [2.5, 0.0, 0.0, 1E-1]), (params=(μ=1E-2,),))     # Δt=5.0, nt=50

export default_parameters, default_constraints

"The magnetic moment μ this equilibrium is set up for."
default_parameters(::Type{T}=Float64) where {T} = (μ = T(1E-2),)

"""
The constraint pair [`hodeproblem`](@ref) and its siblings use here by default.

Both `b₁` and `b₃` vanish at the initial conditions of this equilibrium,
where `b = e₂`, which leaves `(g², g³)` as the only regular pair.
"""
default_constraints() = :g23


include("guiding_center_3d_equations.jl")
include("guiding_center_3d_canonical.jl")

# The canonical toroidal momentum is the covariant φ-component of the momentum. Since p = ϑ on the
# constraint manifold this is simply p₃ — the previous R(t,q) * ϑ₃(t,q) both indexed q[4] of a
# three-component state and carried a spurious factor of R that destroys the conservation.
toroidal_momentum(t, q, p) = p[3]

include("guiding_center_3d_diagnostics.jl")

end
