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

initial_conditions_barely_passing() = merge(initial_conditions(0, [from_cartesian(0, [2.5, 0.1, 0.0])..., 3.425E-1]), (params=(μ=1E-2,),)) # Δt=2.5, nt=50
initial_conditions_barely_trapped() = merge(initial_conditions(0, [from_cartesian(0, [2.5, 0.1, 0.0])..., 3.375E-1]), (params=(μ=1E-2,),)) # Δt=3.0, nt=100
initial_conditions_deeply_passing() = merge(initial_conditions(0, [from_cartesian(0, [2.5, 0.1, 0.0])..., 5E-1]), (params=(μ=1E-2,),))     # Δt=2.5, nt=25
initial_conditions_deeply_trapped() = merge(initial_conditions(0, [from_cartesian(0, [2.5, 0.1, 0.0])..., 1E-1]), (params=(μ=1E-2,),))     # Δt=5.0, nt=50

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

export default_parameters

"The magnetic moment μ this equilibrium is set up for."
default_parameters(::Type{T}=Float64) where {T} = (μ = T(1E-2),)


include("guiding_center_3d_equations.jl")
include("guiding_center_3d_canonical.jl")

# The canonical toroidal momentum is the covariant φ-component of the momentum. Since p = ϑ on the
# constraint manifold this is simply p₃ — the previous R(t,q) * ϑ₃(t,q) both indexed q[4] of a
# three-component state and carried a spurious factor of R that destroys the conservation.
toroidal_momentum(t, q, p) = p[3]

include("guiding_center_3d_diagnostics.jl")

end
