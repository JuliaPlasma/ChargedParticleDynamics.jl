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

    # Left-handed chart: `det(DF) < 0`. See `ORIENTATION` in `gc_common.jl`.
    const ORIENTATION = -1

    include("coordinate_transformations.jl")
    include("gc_common.jl")
    include("gc_equations.jl")

    # This model integrates in rescaled time, `Δs = Δt / ωabs`. The factor is the phasespace Jacobian
    # `ωabs = J B*∥` and *not* the physical `B*∥` — see the note on `ωabs` in `gc_common.jl`. The
    # distinction is at its widest in this chart: at `initial_conditions_deeply_passing`
    # `ωabs = 0.04993` while `B*∥ = 0.9511`, the same `B*∥` the cartesian and cylindrical charts of
    # this equilibrium give. The factor of twenty between them is the volume element `J = 0.0525` of a
    # toroidal chart at `r = 0.05`. The 4D guiding centre runs this equilibrium at `Δt = 500` over
    # `(0, 5E5)`, so `Δs = 500 / 0.04993 ≈ 1E4`.
    #
    # It is deliberately not the 500 the `Cartesian` and `Cylindrical` charts carry: there `ωabs ≈ 1`
    # and the rescaled step nearly coincides with the physical one. Levelling this to 500 would leave
    # the integration looking *better* rather than broken — a thousand Strang steps give 5.8E-4 in
    # relative energy at `Δs = 1E4` against 1.4E-6 at `Δs = 500` — because the second run silently
    # covers a twentieth of the physical span, `2.5E4` against `4.99E5`, not because it resolves the
    # orbit better. That is why the step is pinned by the rescaling assertion in
    # `test/gyro_kinetics_4d_tests.jl` rather than by an energy threshold, which would prefer the
    # wrong value.
    const DEFAULT_TIMESTEP = 1E4
    const DEFAULT_TIMESPAN = (0.0, 1E7)

    const x₀ = from_cartesian(0, [1.05, 0.0, 0.0])

    export default_parameters

    """
    The magnetic moment `μ` of the default initial condition, shared with the 4D guiding centre
    module of the same equilibrium.
    """
    default_parameters(::Type{T}=Float64) where {T} = (μ = T(2.314593645825811e-6),)


    # Shared with the 4D guiding centre module of the same equilibrium, whose `uᵢ` this tracks. The
    # minus this carried until `ElectromagneticFields` 0.7.0 compensated for the reversed `b` of the
    # left-handed `(r, θ, ϕ)` chart; see the note on `uᵢ` in `guiding_center_4d/tokamak_small_toroidal.jl`.
    const qᵢ = [x₀..., 0.00045135897235326736]

    toroidal_momentum(t, q) = ϑ₃(t, q)

    initial_conditions_barely_passing() = (q = [x₀..., 8.117E-4], params = (μ = 2.448E-6,))
    initial_conditions_barely_trapped() = (q = [x₀..., 7.610E-4], params = (μ = 2.250E-6,))
    initial_conditions_deeply_passing() = (q = [x₀..., 1.623E-3], params = (μ = 2.448E-6,))
    initial_conditions_deeply_trapped() = (q = [x₀..., 4.306E-4], params = (μ = 2.250E-6,))
    initial_conditions_pauli()          = (q = [x₀..., 4.3E-4  ], params = (μ = 2.310E-6,))

    include("../guiding_center_4d/guiding_center_4d_diagnostics.jl")

end
