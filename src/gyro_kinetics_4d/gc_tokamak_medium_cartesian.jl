@doc raw"""
Medium-size axisymmetric tokamak equilibrium, in cartesian coordinates, with major radius
``R_{0} = 2``, magnetic field ``B_{0} = 5`` and safety factor ``q_{0} = 2``.

The gyrokinetic counterpart of `GuidingCenter4d.TokamakMediumCartesian`, sharing its initial
conditions. The independent variable is the rescaled time ``s``, with
``dt = B^{\star}_{\parallel} \, ds`` and ``B^{\star}_{\parallel} \approx 3.8`` at the deeply passing
initial condition.

Like its 4D guiding centre counterpart this module defines no `toroidal_momentum`: the third
coordinate is ``z`` rather than an angle, so ``\vartheta_{3}`` is not a conserved momentum here.
"""
module GuidingCenter4dTokamakMediumCartesian

    import ElectromagneticFields.AxisymmetricTokamakCartesian

    export initial_conditions_barely_passing, initial_conditions_barely_trapped,
           initial_conditions_deeply_passing, initial_conditions_deeply_trapped

    export hamiltonian

    AxisymmetricTokamakCartesian.@code(2., 5., 2.) # inject magnetic field code

    include("coordinate_transformations.jl")
    include("gc_common.jl")
    include("gc_equations.jl")

    # Rescaled time: the 4D guiding centre uses Δt = 1.0 over (0, 10³) with B*∥ ≈ 3.8 here.
    const DEFAULT_TIMESTEP = 2.5E-1
    const DEFAULT_TIMESPAN = (0.0, 2.5E2)

    export default_parameters

    """
    The magnetic moment `μ` of the shipped initial conditions.
    """
    default_parameters(::Type{T}=Float64) where {T} = (μ = T(1E-2),)


    const qᵢ = [2.5, 0.0, 0.0, 5E-1]

    initial_conditions_barely_passing() = (q = [2.5, 0.0, 0.0, 3.425E-1], params = (μ = 1E-2,))
    initial_conditions_barely_trapped() = (q = [2.5, 0.0, 0.0, 3.375E-1], params = (μ = 1E-2,))
    initial_conditions_deeply_passing() = (q = [2.5, 0.0, 0.0, 5E-1    ], params = (μ = 1E-2,))
    initial_conditions_deeply_trapped() = (q = [2.5, 0.0, 0.0, 1E-1    ], params = (μ = 1E-2,))

    include("../guiding_center_4d/guiding_center_4d_diagnostics.jl")

end
