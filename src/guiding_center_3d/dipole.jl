"""
Guiding centre dynamics in the field of a magnetic dipole.
"""
module Dipole3d

import ElectromagneticFields.Dipole

export initial_conditions_dipole

export hamiltonian

Dipole.@code() # inject magnetic field code

# A thousand steps, at the same rate as the other 3D guiding centre equilibria.
#
# `Δt = 0.1` is chosen to *exhibit* the constraint drift rather than to minimise it. The constraints
# are this model's characteristic diagnostic — they vanish along the exact flow, so their magnitude
# is the drift off the manifold on which `p = ϑ(q)` — and this equilibrium is the one where the
# three formulations separate most clearly. Over the run below:
#
#   formulation             max|g¹|     max|g²|     max|g³|   max|ΔH/H|
#   hodeproblem            4.8E-06     3.7E-06     3.7E-06     3.3E-07
#   hodeproblem_canonical  7.1E-04     6.9E-04     9.5E-05     2.3E-06
#   hodeproblem_compact    8.1E-06     8.5E-06     5.0E-06     2.1E-07
#
# `Δt = 0.03` puts all of these about two orders lower (`hodeproblem` at 2.2E-08) and the difference
# between the formulations stops being legible, which is the wrong choice for an example. `Δt = 0.3`
# is too far the other way: `hodeproblem_canonical` diverges to 9.1E+3 with 988 of its steps at the
# solver's iteration cap.
const DEFAULT_TIMESTEP = 0.1
const DEFAULT_TIMESPAN = (0.0, 1E2)

initial_conditions_dipole() = merge(initial_conditions(0.0, [1.0, 2.0, 1.0, 0.01]), (params=(μ=1E-2,),))

export default_parameters, default_constraints

"The magnetic moment μ this equilibrium is set up for."
default_parameters(::Type{T}=Float64) where {T} = (μ = T(1E-2),)

"""
The constraint pair [`hodeproblem`](@ref) and its siblings use here by default.

All three pairs are regular at the initial condition — `λₒ` is -34, 34 and 68 — but only this one stays
that way, because the orbit takes `b₁` and `b₂` through zero and the other two pairs divide by them.
Over the module's declared example, a thousand steps at `Δt = 0.1`:

    pair              max|gᵏ|     max|ΔH/H|
    (g¹, g²)  :g12   4.8E-06      3.3E-07
    (g³, g¹)  :g31   7.4E-03      1.3E-06
    (g², g³)  :g23   8.8E-03      9.8E-07

Three orders of magnitude in the constraints, which are this model's characteristic diagnostic.
`(g¹, g²)` divides by `b₃`, which does not vanish along this orbit. The separation is not an artifact
of the coarse step either: at `Δt = 0.03` the same ordering holds at 1.8E-08 against 1.3E-05 and
6.8E-06. What it is a function of is the *span*, since the pairs are all regular at the initial
condition and it is the orbit that spoils two of them — at `Δt = 0.03` and `t = 300`, `:g23` has
reached 8.3E-01 where `:g12` is still at 2.1E-07.

This is the one equilibrium in the package where the pair has to be chosen on its conditioning *along
the orbit* rather than at the initial condition; `TokamakMediumCartesian` is also chosen on
conditioning, but at the initial condition. The full comparison across spans, including
`hodeproblem_canonical` — which loses `:g31` by `t = 30`, well before the Hamilton-Dirac form does —
is tabulated in `docs/src/findings.md` under "Regular is not the same as well conditioned".
"""
default_constraints() = :g12


include("guiding_center_3d_equations.jl")
include("guiding_center_3d_canonical.jl")
include("guiding_center_3d_diagnostics.jl")

end
