"""
Analytic axisymmetric small tokamak equilibrium in circular coordinates.
"""
module TokamakSmallToroidal

import ElectromagneticFields.AxisymmetricTokamakToroidal

export initial_conditions_barely_passing, initial_conditions_barely_trapped,
    initial_conditions_deeply_passing, initial_conditions_deeply_trapped,
    initial_conditions_pauli

export hamiltonian, toroidal_momentum

AxisymmetricTokamakToroidal.@code() # inject magnetic field code

# Δt = 500, matching this equilibrium's `GuidingCenter4d` and `PauliParticle3d` counterparts.
# `hodeproblem` and `hodeproblem_compact` hold it at 4.8E-3 relative energy with constraints of
# 1.5E-8 and 2.2E-8. `hodeproblem_canonical` does not — it throws a `SingularException` at 500 on
# either regular pair — and carries a smaller step where it is integrated; 250 is the largest it
# tolerates. See `test/guiding_center_3d_tests.jl`.
const DEFAULT_TIMESTEP = 500.0
const DEFAULT_TIMESPAN = (0.0, 5E5)
# const DEFAULT_TIMESTEP = 500.0
# const DEFAULT_TIMESPAN = (0.0, 5E4)

const xᵢ = [1.05, 0.0, 0.0]
# As for `TokamakSmallCylindrical`: `(r, θ, ϕ)` is left-handed too, so this carried a compensating
# minus for the reversed `b` that `ElectromagneticFields` 0.7.0 removed the need for. All three charts
# of this equilibrium now start the same physical particle with the same `u` and `μ`.
const uᵢ = 0.00045135897235326736
const qᵢ = [from_cartesian(0, xᵢ)..., uᵢ]


# These four read `[1.05, 0.1, 0.0]` until the three families were aligned. The offset bought nothing
# — this module's pair is `:g12`, which divides by `b₃ = 1.054` and is regular on the midplane, so
# unlike the cartesian charts there is no vanishing component of `b` here to move away from. It was
# also inconsistent three ways over: with `initial_conditions_pauli` just below, with this
# equilibrium's `Cylindrical` chart, and with its `GuidingCenter4d` and `PauliParticle3d` modules. In
# the toroidal chart the offset is not even a toroidal rotation, since it moves `r` from 0.05 to
# 0.0548 and so changes the flux surface. At `y = 0` every formulation holds the same step or better:
# `hodeproblem` on `barely_passing` improves from 4.7E-7 to 6.4E-8 relative energy.
initial_conditions_barely_passing() = merge(initial_conditions(0, [from_cartesian(0, [1.05, 0.0, 0.0])..., 8.117E-4]), (params=(μ=2.448E-6,),))
initial_conditions_barely_trapped() = merge(initial_conditions(0, [from_cartesian(0, [1.05, 0.0, 0.0])..., 7.610E-4]), (params=(μ=2.250E-6,),))
initial_conditions_deeply_passing() = merge(initial_conditions(0, [from_cartesian(0, [1.05, 0.0, 0.0])..., 1.623E-3]), (params=(μ=2.448E-6,),))
initial_conditions_deeply_trapped() = merge(initial_conditions(0, [from_cartesian(0, [1.05, 0.0, 0.0])..., 4.306E-4]), (params=(μ=2.250E-6,),))
initial_conditions_pauli() = merge(initial_conditions(0, [from_cartesian(0, [1.05, 0.0, 0.0])..., 4.3E-4]), (params=(μ=2.310E-6,),))

u_loop() = 4.0E-4
μ_loop() = 2.5E-6
u_surface() = 4.0E-4
μ_surface() = 2.5E-6

function f_loop(t)
    R0 = 1.0
    Z0 = 0.0
    φ0 = 0.0
    u0 = u_loop()
    r0 = 0.05

    Rt = R0 + r0 * cos(2π * t)
    Zt = Z0 + r0 * sin(2π * t)

    qt = [Rt, Zt, φ0, u0]

    return qt
end

function f_surface(s, t)
    R0 = 1.0
    Z0 = 0.0
    φ0 = 0.0
    u0 = u_surface()
    r0 = 0.1

    Rt = R0 + 2r0 * (s - 0.5)
    Zt = Z0 + 2r0 * (t - 0.5)

    qt = [Rt, Zt, φ0, u0]

    return qt
end

export default_parameters, default_constraints

"The magnetic moment μ this equilibrium is set up for."
default_parameters(::Type{T}=Float64) where {T} = (μ = T(2.314593645825811e-6),)

"""
The constraint pair [`hodeproblem`](@ref) and its siblings use here by default.

`b₁ = b_r` vanishes on the midplane, where every initial condition of
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

end
