"""
Analytic axisymmetric small tokamak equilibrium in cylindrical coordinates.
"""
module TokamakSmallCylindrical

import ElectromagneticFields.AxisymmetricTokamakCylindrical

export initial_conditions_barely_passing, initial_conditions_barely_trapped,
    initial_conditions_deeply_passing, initial_conditions_deeply_trapped,
    initial_conditions_pauli

export hamiltonian, toroidal_momentum

AxisymmetricTokamakCylindrical.@code() # inject magnetic field code

const DEFAULT_TIMESTEP = 0.1
const DEFAULT_TIMESPAN = (0.0, 1E2)
# const DEFAULT_TIMESTEP = 500.0
# const DEFAULT_TIMESPAN = (0.0, 5E4)

const xᵢ = [1.05, 0.0, 0.0]
const uᵢ = -0.00045135897235326736
const qᵢ = [from_cartesian(0, xᵢ)..., uᵢ]


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

export default_parameters

"The magnetic moment μ this equilibrium is set up for."
default_parameters(::Type{T}=Float64) where {T} = (μ = T(2.314593645825811e-6),)

const parameters = default_parameters()

include("guiding_center_3d_equations.jl")
include("guiding_center_3d_canonical.jl")

# The canonical toroidal momentum is the covariant φ-component of the momentum. Since p = ϑ on the
# constraint manifold this is simply p₃ — the previous R(t,q) * ϑ₃(t,q) both indexed q[4] of a
# three-component state and carried a spurious factor of R that destroys the conservation.
toroidal_momentum(t, q, p) = p[3]

include("guiding_center_3d_diagnostics.jl")

end
