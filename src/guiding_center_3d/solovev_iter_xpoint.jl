"""
Analytic ITER-like Solov'ev equilibrium with X-point.
"""
module SolovevIterXpoint

import ElectromagneticFields.Solovev

export initial_conditions_barely_passing, initial_conditions_barely_trapped,
    initial_conditions_deeply_passing, initial_conditions_deeply_trapped,
    initial_conditions_trapped

export hamiltonian, toroidal_momentum

Solovev.@code_iter_xpoint() # inject magnetic field code

# A thousand steps of the shipped initial conditions. The relative energy error over the run is
# 1.0E-2, against 1.6E-2 at `Δt = 0.2` and 2.7 at `Δt = 0.5`; `Δt = 1` is not integrable here at all
# (`hodeproblem` reaches `NaN` with 89 of its steps at the solver's iteration cap, and
# `hodeproblem_canonical` throws).
#
# The span is what a thousand steps buys, not a limit of one formulation: over `t ∈ [0, 1E3]` all
# three lose `barely_passing`, `barely_trapped` and `deeply_passing` alike (2.7, 2.7, and `NaN` or
# 5.8), while `deeply_trapped` and `trapped` survive. The orbit runs out, not a discretisation of
# it, so the example is short rather than one formulation carrying an exception.
const DEFAULT_TIMESTEP = 0.1
const DEFAULT_TIMESPAN = (0.0, 1E2)

const xᵢ = [7.0 - 1.4, 0.0, 0.0]
const qᵢ = [from_cartesian(0, xᵢ)..., 2.8166280889939737]

export default_parameters, default_constraints

"The magnetic moment μ this equilibrium is set up for."
default_parameters(::Type{T}=Float64) where {T} = (μ = T(4.607782183567846),)

"""
The constraint pair [`hodeproblem`](@ref) and its siblings use here by default.

`b₁` is small but non-zero at every initial condition of this
equilibrium, so the historical pair is still usable.
"""
default_constraints() = :g31


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
