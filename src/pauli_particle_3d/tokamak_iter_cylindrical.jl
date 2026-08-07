@doc raw"""
ITER-size axisymmetric tokamak equilibrium, in cylindrical coordinates ``(R, Z, \varphi)``.

Major radius ``R_{0} = 6.2``, magnetic field on axis ``B_{0} = 5.3`` and safety factor
``q_{0} = \sqrt{2}`` — the ITER parameters of `ElectromagneticFields`.

The hard-coded `μ` and `u` of the `initial_conditions_*` below are inherited from earlier work on
this package and their provenance is not recorded; they are known-good starting points for the
Pauli particle model rather than values derived here.
"""
module TokamakIterCylindrical

    import ElectromagneticFields.AxisymmetricTokamakCylindrical

    export podeproblem, hamiltonian, toroidal_momentum
    export initial_conditions_barely_passing, initial_conditions_barely_trapped
    export initial_conditions_deeply_passing, initial_conditions_deeply_trapped
    export initial_conditions_trapped

    AxisymmetricTokamakCylindrical.@code_iter() # inject magnetic field code

    const qᵢ = [7.0-1.4, 0.0, 0.0]
    const vᵢ = [3.43E-3, 6.75, -3.41E-1]

    include("pauli_particle_3d.jl")

    export default_parameters

    """
    The magnetic moment μ of the default initial condition `(qᵢ, vᵢ)`, obtained by splitting
    `vᵢ` into its parallel and perpendicular parts at `qᵢ`.
    """
    default_parameters(::Type{T}=Float64) where {T} = (μ = T(initial_conditions(qᵢ, vᵢ).params.μ),)


    const DEFAULT_TIMESTEP = 1.0
    const DEFAULT_TIMESPAN = (0.0, 1E3)

    # The same `(x, u, μ)` as this equilibrium's `GuidingCenter3d` and `GuidingCenter4d` modules.
    # `barely_passing` was missing, which is why nothing exercised the condition the other two
    # families lead with.
    initial_conditions_barely_passing() = initial_conditions(from_cartesian(0, [2.5, 0., 0.]), 3.425E-1, 1E-2)
    initial_conditions_barely_trapped() = initial_conditions(from_cartesian(0, [2.5, 0., 0.]), 3.375E-1, 1E-2)
    initial_conditions_deeply_passing() = initial_conditions(from_cartesian(0, [2.5, 0., 0.]),  5E-1,    1E-2)
    initial_conditions_deeply_trapped() = initial_conditions(from_cartesian(0, [2.5, 0., 0.]),  1E-1,    1E-2)
    initial_conditions_trapped()        = initial_conditions(from_cartesian(0, [7.0, 0., 0.]), -2E-3,    1.88E-7)

end
