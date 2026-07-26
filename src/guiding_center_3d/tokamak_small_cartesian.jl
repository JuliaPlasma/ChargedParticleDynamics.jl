"""
Analytic axisymmetric small tokamak equilibrium in cartesian coordinates.
"""
module TokamakSmallCartesian

import ElectromagneticFields.AxisymmetricTokamakCartesian

export initial_conditions_barely_passing, initial_conditions_barely_trapped,
    initial_conditions_deeply_passing, initial_conditions_deeply_trapped,
    initial_conditions_pauli, initial_conditions_default

export hamiltonian, toroidal_momentum

AxisymmetricTokamakCartesian.@code() # inject magnetic field code

# const DEFAULT_TIMESTEP = 0.1
# const DEFAULT_TIMESPAN = (0.0, 1E2)

const DEFAULT_TIMESTEP = 200.0
const DEFAULT_TIMESPAN = (0.0, 2E6)

# const DEFAULT_TIMESTEP = 500.0
# const DEFAULT_TIMESPAN = (0.0, 5E6)

initial_conditions_default() = merge(initial_conditions(0, [1.05, 0.0, -0.005, 0.000045]), (params=(μ=3.2e-7,),))
initial_conditions_barely_passing() = merge(initial_conditions(0, [1.05, 0.0, 0.0, 8.117E-4]), (params=(μ=2.448E-6,),))
initial_conditions_barely_trapped() = merge(initial_conditions(0, [1.05, 0.0, 0.0, 7.610E-4]), (params=(μ=2.250E-6,),))
initial_conditions_deeply_passing() = merge(initial_conditions(0, [1.05, 0.0, 0.0, 1.623E-3]), (params=(μ=2.448E-6,),))
initial_conditions_deeply_trapped() = merge(initial_conditions(0, [1.05, 0.0, 0.0, 4.306E-4]), (params=(μ=2.250E-6,),))
initial_conditions_pauli() = merge(initial_conditions(0, [1.05, 0.0, 0.0, 4.3E-4]), (params=(μ=2.310E-6,),))

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
default_parameters(::Type{T}=Float64) where {T} = (μ = T(3.2e-7),)

const parameters = default_parameters()

include("guiding_center_3d_equations.jl")
include("guiding_center_3d_canonical.jl")

# In cartesian coordinates the third coordinate is z, not an angle, so the toroidal momentum is the
# generator of rotation about the z-axis, x p₂ - y p₁, rather than p₃.
toroidal_momentum(t, q, p) = q[1] * p[2] - q[2] * p[1]

include("guiding_center_3d_diagnostics.jl")

end
