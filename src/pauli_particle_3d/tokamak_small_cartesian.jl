@doc raw"""
Small axisymmetric tokamak equilibrium, in cartesian coordinates ``(x, y, z)``.

Major radius ``R_{0} = 1``, magnetic field on axis ``B_{0} = 1`` and safety factor
``q_{0} = 2``.

The hard-coded `μ` and `u` of the `initial_conditions_*` below are inherited from earlier work on
this package and their provenance is not recorded; they are known-good starting points for the
Pauli particle model rather than values derived here.
"""
module TokamakSmallCartesian

    import ElectromagneticFields.AxisymmetricTokamakCartesian

    export podeproblem, hamiltonian, toroidal_momentum

    const equ = AxisymmetricTokamakCartesian.@code() # inject magnetic field code

    const qᵢ = [1.05,   0.0,    0.0]
    const vᵢ = [2.1E-3, 4.3E-4, 0.0]

    const DEFAULT_TIMESTEP = 500.0
    const DEFAULT_TIMESPAN = (0.0, 5E4)

    include("pauli_particle_3d.jl")

    export default_parameters

    """
    The magnetic moment μ of the default initial condition `(qᵢ, vᵢ)`, obtained by splitting
    `vᵢ` into its parallel and perpendicular parts at `qᵢ`.
    """
    default_parameters(::Type{T}=Float64) where {T} = (μ = T(initial_conditions(qᵢ, vᵢ).params.μ),)

    const parameters = default_parameters()

end
