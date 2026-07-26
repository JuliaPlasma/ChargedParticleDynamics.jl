@doc raw"""
Medium-size axisymmetric tokamak equilibrium, in cylindrical coordinates ``(R, Z, \varphi)``, with
major radius ``R_{0} = 2``, magnetic field ``B_{0} = 5`` and safety factor ``q_{0} = 2``.

The gyrokinetic counterpart of [`GuidingCenter4d.TokamakMediumCylindrical`](@ref), sharing its
initial conditions. The independent variable is the rescaled time ``s``, with
``dt = B^{\star}_{\parallel} \, ds`` and ``B^{\star}_{\parallel} \approx 11`` at the deeply passing
initial condition.
"""
module GuidingCenter4dTokamakMediumCylindrical

    import ElectromagneticFields.AxisymmetricTokamakCylindrical

    export initial_conditions_barely_passing, initial_conditions_barely_trapped,
           initial_conditions_deeply_passing, initial_conditions_deeply_trapped

    export hamiltonian, toroidal_momentum

    AxisymmetricTokamakCylindrical.@code(2., 5., 2.) # inject magnetic field code

    include("coordinate_transformations.jl")
    include("gc_common.jl")
    include("gc_equations.jl")

    # Rescaled time: the 4D guiding centre uses Δt = 1.0 over (0, 10³) with B*∥ ≈ 11 here.
    const DEFAULT_TIMESTEP = 1E-1
    const DEFAULT_TIMESPAN = (0.0, 1E2)

    const x₀ = from_cartesian(0, [2.5, 0.0, 0.0])

    export default_parameters

    """
    The magnetic moment `μ` of the shipped initial conditions.
    """
    default_parameters(::Type{T}=Float64) where {T} = (μ = T(1E-2),)

    const parameters = default_parameters()

    const qᵢ = [x₀..., 5E-1]

    toroidal_momentum(t, q) = ϑ₃(t, q)

    initial_conditions_barely_passing() = (q = [x₀..., 3.425E-1], params = (μ = 1E-2,))
    initial_conditions_barely_trapped() = (q = [x₀..., 3.375E-1], params = (μ = 1E-2,))
    initial_conditions_deeply_passing() = (q = [x₀..., 5E-1    ], params = (μ = 1E-2,))
    initial_conditions_deeply_trapped() = (q = [x₀..., 1E-1    ], params = (μ = 1E-2,))

    include("../guiding_center_4d/guiding_center_4d_diagnostics.jl")

end
