"""
Analytic axisymmetric medium-size tokamak equilibrium in cartesian coordinates.
"""
module TokamakMediumCylindrical

import ElectromagneticFields.AxisymmetricTokamakCylindrical

export initial_conditions_barely_passing, initial_conditions_barely_trapped,
    initial_conditions_deeply_passing, initial_conditions_deeply_trapped

export hamiltonian, toroidal_momentum

AxisymmetricTokamakCylindrical.@code(2.0, 5.0, 2.0) # inject magnetic field code

const DEFAULT_TIMESTEP = 0.1
const DEFAULT_TIMESPAN = (0.0, 1E2)

# These read `[2.5, 0.1, 0.0]` until the families were aligned, matching neither this equilibrium's
# `GuidingCenter4d` module nor its own `Cartesian` chart's reason for the offset. In a cylindrical
# chart a cartesian `y` displacement is very nearly a toroidal rotation — it leaves `Z = 0` and so
# leaves `b₁ = b_R = 0` either way — so it cannot make a singular pair regular, which is what the
# offset is for in the cartesian chart. All three formulations agree to three digits at both values.
initial_conditions_barely_passing() = merge(initial_conditions(0, [from_cartesian(0, [2.5, 0.0, 0.0])..., 3.425E-1]), (params=(μ=1E-2,),))
initial_conditions_barely_trapped() = merge(initial_conditions(0, [from_cartesian(0, [2.5, 0.0, 0.0])..., 3.375E-1]), (params=(μ=1E-2,),))
initial_conditions_deeply_passing() = merge(initial_conditions(0, [from_cartesian(0, [2.5, 0.0, 0.0])..., 5E-1]), (params=(μ=1E-2,),))
initial_conditions_deeply_trapped() = merge(initial_conditions(0, [from_cartesian(0, [2.5, 0.0, 0.0])..., 1E-1]), (params=(μ=1E-2,),))

μ_loop() = 1E-3
μ_surface() = 1E-3

function f_loop(s)
    R0 = 1.75
    Z0 = 0.0
    φ0 = 0.0
    u0 = 0.5
    rx = 0.1
    ry = 0.1

    Rs = R0 + rx * cos(2π * s)
    Zs = Z0 + ry * sin(2π * s)

    qs = [Rs, Zs, φ0, u0]

    return qs
end

function f_surface(s, t)
    R0 = 1.75
    Z0 = 0.0
    φ0 = 0.0
    u0 = 0.5
    r0 = 0.1

    Rt = R0 + r0 * (s - 0.5)
    Zt = Z0 + r0 * (t - 0.5)

    qt = [Rt, Zt, φ0, u0]

    return qt
end

export default_parameters, default_constraints

"The magnetic moment μ this equilibrium is set up for."
default_parameters(::Type{T}=Float64) where {T} = (μ = T(1E-2),)

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

end
