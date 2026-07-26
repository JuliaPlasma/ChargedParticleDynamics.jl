@doc raw"""
Small axisymmetric tokamak equilibrium, in toroidal coordinates ``(r, \theta, \varphi)``.

The gyrokinetic counterpart of `GuidingCenter4d.TokamakSmallToroidal`, sharing its initial
conditions. The independent variable is the rescaled time ``s``, with
``dt = B^{\star}_{\parallel} \, ds`` and ``B^{\star}_{\parallel} \approx 5 \times 10^{-2}`` at the
deeply passing initial condition — the smallest of the shipped equilibria, so the rescaled step is
correspondingly the largest.
"""
module GuidingCenter4dTokamakSmallToroidal

    import ElectromagneticFields.AxisymmetricTokamakToroidal

    export initial_conditions_barely_passing, initial_conditions_barely_trapped,
           initial_conditions_deeply_passing, initial_conditions_deeply_trapped,
           initial_conditions_pauli

    export hamiltonian, toroidal_momentum

    AxisymmetricTokamakToroidal.@code() # inject magnetic field code

    include("coordinate_transformations.jl")
    include("gc_common.jl")
    include("gc_equations.jl")

    # Rescaled time: the 4D guiding centre uses Δt = 500 over (0, 5E4) with B*∥ ≈ 0.05 here.
    const DEFAULT_TIMESTEP = 1E4
    const DEFAULT_TIMESPAN = (0.0, 1E6)

    const x₀ = from_cartesian(0, [1.05, 0.0, 0.0])

    export default_parameters

    """
    The magnetic moment `μ` of the default initial condition, shared with the 4D guiding centre
    module of the same equilibrium.
    """
    default_parameters(::Type{T}=Float64) where {T} = (μ = T(2.314593645825811e-6),)


    const qᵢ = [x₀..., -0.00045135897235326736]

    toroidal_momentum(t, q) = ϑ₃(t, q)

    initial_conditions_barely_passing() = (q = [x₀..., 8.117E-4], params = (μ = 2.448E-6,))
    initial_conditions_barely_trapped() = (q = [x₀..., 7.610E-4], params = (μ = 2.250E-6,))
    initial_conditions_deeply_passing() = (q = [x₀..., 1.623E-3], params = (μ = 2.448E-6,))
    initial_conditions_deeply_trapped() = (q = [x₀..., 4.306E-4], params = (μ = 2.250E-6,))
    initial_conditions_pauli()          = (q = [x₀..., 4.3E-4  ], params = (μ = 2.310E-6,))

    include("../guiding_center_4d/guiding_center_4d_diagnostics.jl")

end
