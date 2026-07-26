@doc raw"""
Small axisymmetric tokamak equilibrium, in cartesian coordinates.

The gyrokinetic counterpart of `GuidingCenter4d.TokamakSmallCartesian`, sharing its initial
conditions. The independent variable is the rescaled time ``s``, with
``dt = B^{\star}_{\parallel} \, ds``; here ``B^{\star}_{\parallel} \approx 0.95`` at the deeply
passing initial condition, so the rescaling is close to a no-op and the defaults nearly match the
4D guiding centre's.
"""
module GuidingCenter4dTokamakSmallCartesian

    import ElectromagneticFields.AxisymmetricTokamakCartesian

    export initial_conditions_barely_passing, initial_conditions_barely_trapped,
           initial_conditions_deeply_passing, initial_conditions_deeply_trapped,
           initial_conditions_pauli

    export hamiltonian, toroidal_momentum

    AxisymmetricTokamakCartesian.@code() # inject magnetic field code

    include("coordinate_transformations.jl")
    include("gc_common.jl")
    include("gc_equations.jl")

    # Rescaled time: the 4D guiding centre uses Δt = 500 over (0, 5E4) with B*∥ ≈ 0.95 here.
    const DEFAULT_TIMESTEP = 5E2
    const DEFAULT_TIMESPAN = (0.0, 5E4)

    export default_parameters

    """
    The magnetic moment `μ` of the default initial condition, shared with the 4D guiding centre
    module of the same equilibrium.
    """
    default_parameters(::Type{T}=Float64) where {T} = (μ = T(2.314593645825811e-6),)


    const qᵢ = [1.05, 0.0, 0.0, 0.00045135897235326736]

    # In cartesian coordinates the third coordinate is `z` rather than an angle, so the conserved
    # quantity is the generator of rotation about the z-axis, not `ϑ₃`.
    toroidal_momentum(t, q) = q[1] * ϑ₂(t, q) - q[2] * ϑ₁(t, q)

    initial_conditions_barely_passing() = (q = [1.05, 0.0, 0.0, 8.117E-4], params = (μ = 2.448E-6,))
    initial_conditions_barely_trapped() = (q = [1.05, 0.0, 0.0, 7.610E-4], params = (μ = 2.250E-6,))
    initial_conditions_deeply_passing() = (q = [1.05, 0.0, 0.0, 1.623E-3], params = (μ = 2.448E-6,))
    initial_conditions_deeply_trapped() = (q = [1.05, 0.0, 0.0, 4.306E-4], params = (μ = 2.250E-6,))
    initial_conditions_pauli()          = (q = [1.05, 0.0, 0.0, 4.3E-4  ], params = (μ = 2.310E-6,))

    include("../guiding_center_4d/guiding_center_4d_diagnostics.jl")

end
