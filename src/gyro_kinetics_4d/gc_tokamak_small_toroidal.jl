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

    # This model integrates in rescaled time, `Δs = Δt / ωabs`. The factor is `ωabs` and *not* the
    # physical `B*∥` — see the note on `ωabs` in `gc_common.jl`. The distinction is at its widest in
    # this chart: at `initial_conditions_deeply_passing` `ωabs = -0.04993` while `B*∥ = +0.9511`, the
    # same `B*∥` the cartesian and cylindrical charts of this equilibrium give. The factor of twenty
    # is the volume element `J = 0.0525` of a toroidal chart at `r = 0.05`, and the sign is its
    # handedness. The 4D guiding centre runs this equilibrium at `Δt = 500` over `(0, 5E5)`, so
    # `|Δs| = 500 / 0.04993 ≈ 1E4`.
    #
    # It is deliberately not the 500 the `Cartesian` and `Cylindrical` charts carry: there `|ωabs| ≈ 1`
    # and the rescaled step nearly coincides with the physical one. Levelling this to 500 leaves the
    # integration looking healthy — a thousand steps at either value gives the same 4.277E-3 relative
    # energy error, since each covers a different physical span — and is caught only by the rescaling
    # assertion in `test/gyro_kinetics_4d_tests.jl`.
    #
    # The step below is positive, so `s` advances while physical time runs *backwards*, as in the
    # cylindrical chart. `TODO.md` carries the question of whether to rescale by `B*∥` instead.
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
