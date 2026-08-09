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
    export initial_conditions_barely_passing, initial_conditions_barely_trapped
    export initial_conditions_deeply_passing, initial_conditions_deeply_trapped
    export initial_conditions_pauli

    const equ = AxisymmetricTokamakCartesian.@code() # inject magnetic field code

    const qᵢ = [1.05,   0.0,    0.0]
    const vᵢ = [2.1E-3, 4.3E-4, 0.0]

    # A thousand steps, at which all three formulations agree to a relative energy error of 4.5E-8.
    # `Δt = 500` is not usable here: the `pode` and `hode` forms reach 1.9E+37 with 441 of their
    # steps at the solver's iteration cap, and the `iode` form throws.
    #
    # The `Cylindrical` and `Toroidal` charts of this same tokamak do run at `Δt = 500`, at 2.5E-6
    # and 1.2E-6. The asymmetry is real: it is the cartesian chart that is stiff, not the
    # equilibrium — the toroidal transit appears there as a rotation of (x, y) that has to be
    # resolved, where in the other two charts it is a coordinate that simply winds.
    const DEFAULT_TIMESTEP = 10.0
    const DEFAULT_TIMESPAN = (0.0, 1E4)

    include("pauli_particle_3d.jl")

    export default_parameters

    """
    The magnetic moment μ of the default initial condition `(qᵢ, vᵢ)`, obtained by splitting
    `vᵢ` into its parallel and perpendicular parts at `qᵢ`.
    """
    default_parameters(::Type{T}=Float64) where {T} = (μ = T(initial_conditions(qᵢ, vᵢ).params.μ),)


    # The same `(x, u, μ)` as this equilibrium's `GuidingCenter3d` and `GuidingCenter4d` modules, so
    # that the three families can be started from one condition and compared. This module had none,
    # which is why the Pauli family could only ever be run at its own `(qᵢ, vᵢ)` default.
    initial_conditions_barely_passing() = initial_conditions([1.05, 0.0, 0.0], 8.117E-4, 2.448E-6)
    initial_conditions_barely_trapped() = initial_conditions([1.05, 0.0, 0.0], 7.610E-4, 2.250E-6)
    initial_conditions_deeply_passing() = initial_conditions([1.05, 0.0, 0.0], 1.623E-3, 2.448E-6)
    initial_conditions_deeply_trapped() = initial_conditions([1.05, 0.0, 0.0], 4.306E-4, 2.250E-6)

    # The guiding centre modules call this one `initial_conditions_pauli` because it is the condition
    # that matches this family's own `(qᵢ, vᵢ)` default, to the three digits they carry it to.
    initial_conditions_pauli() = initial_conditions([1.05, 0.0, 0.0], 4.3E-4, 2.310E-6)

end
