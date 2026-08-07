@doc raw"""
Small axisymmetric tokamak equilibrium, in cylindrical coordinates ``(R, Z, \varphi)``.

The gyrokinetic counterpart of `GuidingCenter4d.TokamakSmallCylindrical`, sharing its
initial conditions. The independent variable is the rescaled time ``s``, with
``dt = B^{\star}_{\parallel} \, ds``; here ``B^{\star}_{\parallel} \approx 1.0`` at the deeply
passing initial condition, so the rescaling is close to a no-op.
"""
module GuidingCenter4dTokamakSmallCylindrical

    import ElectromagneticFields.AxisymmetricTokamakCylindrical

    export initial_conditions_barely_passing, initial_conditions_barely_trapped,
           initial_conditions_deeply_passing, initial_conditions_deeply_trapped,
           initial_conditions_pauli

    export hamiltonian, toroidal_momentum

    AxisymmetricTokamakCylindrical.@code() # inject magnetic field code

    include("coordinate_transformations.jl")
    include("gc_common.jl")
    include("gc_equations.jl")

    # Rescaled time: the 4D guiding centre uses Δt = 500 over (0, 5E4) with B*∥ ≈ 1.0 here.
    const DEFAULT_TIMESTEP = 5E2
    const DEFAULT_TIMESPAN = (0.0, 5E5)

    const x₀ = from_cartesian(0, [1.05, 0.0, 0.0])

    export default_parameters

    """
    The magnetic moment `μ` of the default initial condition, shared with the 4D guiding centre
    module of the same equilibrium.
    """
    default_parameters(::Type{T}=Float64) where {T} = (μ = T(2.314593645825811e-6),)


    # Shared with the 4D guiding centre module of the same equilibrium, whose `uᵢ` this tracks. The
    # minus this carried until `ElectromagneticFields` 0.7.0 compensated for the reversed `b` of the
    # left-handed `(R, Z, ϕ)` chart; see the note on `uᵢ` in `guiding_center_4d/tokamak_small_cylindrical.jl`.
    const qᵢ = [x₀..., 0.00045135897235326736]

    toroidal_momentum(t, q) = ϑ₃(t, q)

    initial_conditions_barely_passing() = (q = [x₀..., 8.117E-4], params = (μ = 2.448E-6,))
    initial_conditions_barely_trapped() = (q = [x₀..., 7.610E-4], params = (μ = 2.250E-6,))
    initial_conditions_deeply_passing() = (q = [x₀..., 1.623E-3], params = (μ = 2.448E-6,))
    initial_conditions_deeply_trapped() = (q = [x₀..., 4.306E-4], params = (μ = 2.250E-6,))
    initial_conditions_pauli()          = (q = [x₀..., 4.3E-4  ], params = (μ = 2.310E-6,))

    include("../guiding_center_4d/guiding_center_4d_diagnostics.jl")

end
