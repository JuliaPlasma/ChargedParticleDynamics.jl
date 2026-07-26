"""
Analytic axisymmetric medium-size tokamak equilibrium in cartesian coordinates.
"""
module TokamakMediumCartesian

import ElectromagneticFields.AxisymmetricTokamakCartesian

export initial_conditions_barely_passing, initial_conditions_barely_trapped,
    initial_conditions_deeply_passing, initial_conditions_deeply_trapped

export hamiltonian

AxisymmetricTokamakCartesian.@code(2.0, 5.0, 2.0) # inject magnetic field code

const DEFAULT_TIMESTEP = 0.1
const DEFAULT_TIMESPAN = (0.0, 1E2)

initial_conditions_barely_passing() = merge(initial_conditions(0, [2.5, 0.1, 0.0, 3.425E-1]), (params=(μ=1E-2,),)) # Δt=2.5, nt=50
initial_conditions_barely_trapped() = merge(initial_conditions(0, [2.5, 0.1, 0.0, 3.375E-1]), (params=(μ=1E-2,),)) # Δt=3.0, nt=100
initial_conditions_deeply_passing() = merge(initial_conditions(0, [2.5, 0.1, 0.0, 5E-1]), (params=(μ=1E-2,),))     # Δt=2.5, nt=25
initial_conditions_deeply_trapped() = merge(initial_conditions(0, [2.5, 0.1, 0.0, 1E-1]), (params=(μ=1E-2,),))     # Δt=5.0, nt=50

μ_loop() = 1E-3
μ_surface() = 1E-3

function f_loop(s)
    X0 = 1.75
    Y0 = 0.0
    Z0 = 0.0
    u0 = 0.5
    rx = 0.1
    ry = 0.1

    Xs = X0 + rx * cos(2π * s)
    Zs = Z0 + ry * sin(2π * s)

    qs = [Xs, Y0, Zs, u0]

    return qs
end

function f_surface(s, t)
    X0 = 1.75
    Y0 = 0.0
    Z0 = 0.0
    u0 = 0.5
    r0 = 0.1

    Xt = X0 + r0 * (s - 0.5)
    Zt = Z0 + r0 * (t - 0.5)

    qt = [Xt, Y0, Zt, u0]

    return qt
end

export default_parameters, default_constraints

"The magnetic moment μ this equilibrium is set up for."
default_parameters(::Type{T}=Float64) where {T} = (μ = T(1E-2),)

"""
The constraint pair [`hodeproblem`](@ref) and its siblings use here by default.

`b₁ ≠ 0` at every initial condition of this equilibrium, so the
historical pair is still usable.
"""
default_constraints() = :g31


include("guiding_center_3d_equations.jl")
include("guiding_center_3d_canonical.jl")
include("guiding_center_3d_diagnostics.jl")

end
