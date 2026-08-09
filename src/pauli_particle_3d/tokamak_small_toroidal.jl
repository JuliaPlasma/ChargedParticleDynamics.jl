@doc raw"""
Small axisymmetric tokamak equilibrium, in toroidal coordinates ``(r, \theta, \varphi)``.

Major radius ``R_{0} = 1``, magnetic field on axis ``B_{0} = 1`` and safety factor
``q_{0} = 2``.

The hard-coded `μ` and `u` of the `initial_conditions_*` below are inherited from earlier work on
this package and their provenance is not recorded; they are known-good starting points for the
Pauli particle model rather than values derived here.
"""
module TokamakSmallToroidal

    import ElectromagneticFields.AxisymmetricTokamakToroidal

    export podeproblem, hamiltonian, toroidal_momentum
    export initial_conditions_barely_passing, initial_conditions_barely_trapped
    export initial_conditions_deeply_passing, initial_conditions_deeply_trapped
    export initial_conditions_pauli

    AxisymmetricTokamakToroidal.@code() # inject magnetic field code

    const xᵢ = [1.05, 0.0, 0.0]
    const qᵢ = from_cartesian(0, xᵢ)
    # `DF̄` is the inverse Jacobian of the chart map, so it takes the chart coordinates `qᵢ = (r, θ, ϕ)`
    # and not the cartesian `xᵢ`. Here the two genuinely differ — `r = 0.05` against `x = 1.05` — so
    # this is the chart in which the question the comment used to ask has an answer that matters.
    const vᵢ = DF̄(0, qᵢ) * [2.1E-3, 4.3E-4, 0.0]

    const DEFAULT_TIMESTEP = 500.0
    const DEFAULT_TIMESPAN = (0.0, 5E5)

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
    #
    # The `u` here needed no adjustment for `ElectromagneticFields` 0.7.0. It is the same value the
    # cartesian chart of this equilibrium carries, and now denotes the same physical particle in both
    # — which it did not while `b` was reversed here; see `uᵢ` in
    # `guiding_center_4d/tokamak_small_toroidal.jl`.
    initial_conditions_barely_passing() = initial_conditions(from_cartesian(0, [1.05, 0.0, 0.0]), 8.117E-4, 2.448E-6)
    initial_conditions_barely_trapped() = initial_conditions(from_cartesian(0, [1.05, 0.0, 0.0]), 7.610E-4, 2.250E-6)
    initial_conditions_deeply_passing() = initial_conditions(from_cartesian(0, [1.05, 0.0, 0.0]), 1.623E-3, 2.448E-6)
    initial_conditions_deeply_trapped() = initial_conditions(from_cartesian(0, [1.05, 0.0, 0.0]), 4.306E-4, 2.250E-6)

    # The guiding centre modules call this one `initial_conditions_pauli` because it is the condition
    # that matches this family's own `(qᵢ, vᵢ)` default, to the three digits they carry it to.
    initial_conditions_pauli() = initial_conditions(from_cartesian(0, [1.05, 0.0, 0.0]), 4.3E-4, 2.310E-6)

end
