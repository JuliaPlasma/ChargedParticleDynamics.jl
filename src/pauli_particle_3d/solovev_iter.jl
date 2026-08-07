@doc raw"""
Analytic ITER-like Solov'ev equilibrium, in cylindrical coordinates ``(R, Z, \varphi)``.

The Solov'ev solution of the Grad-Shafranov equation fitted to ITER's shape, without an
X-point. See `SolovevIterXpoint` for the diverted variant.

The hard-coded `μ` and `u` of the `initial_conditions_*` below are inherited from earlier work on
this package and their provenance is not recorded; they are known-good starting points for the
Pauli particle model rather than values derived here.
"""
module SolovevIter

    import ElectromagneticFields.Solovev

    export podeproblem, hamiltonian, toroidal_momentum
    export initial_conditions_barely_passing, initial_conditions_barely_trapped
    export initial_conditions_deeply_passing, initial_conditions_deeply_trapped

    # const R₀ = 6.2
    # const B₀ = 5.3
    # const E₀ = 0.
    # const q  = 2.

    Solovev.@code_iter() # inject magnetic field code

    const xᵢ = [7.0-1.4, 0.0, 0.0]
    const qᵢ = from_cartesian(0, xᵢ)
    const vᵢ = [3.43E-3, 6.75, -3.41E-1]

    # A thousand steps, at a relative energy error of 9.7E-3 (1.1E-2 at `Δt = 0.2`, 9.4E-2 at
    # `Δt = 0.5`). `Δt = 1` is not usable: all three formulations throw `NaN detected in direction
    # vector`.
    #
    # The span is what a thousand steps buys, not a limit of one formulation: over `t ∈ [0, 1E3]`
    # `pode`, `hode` and `iode` alike reach 2.7 for `barely_passing`, 2.7 for `barely_trapped` and
    # 5.6 for `deeply_passing`, while `deeply_trapped` and `trapped` survive. The orbit runs out,
    # not a discretisation of it.
    const DEFAULT_TIMESTEP = 0.1
    const DEFAULT_TIMESPAN = (0.0, 1E2)

    include("pauli_particle_3d.jl")

    export default_parameters

    """
    The magnetic moment μ of the default initial condition `(qᵢ, vᵢ)`, obtained by splitting
    `vᵢ` into its parallel and perpendicular parts at `qᵢ`.
    """
    default_parameters(::Type{T}=Float64) where {T} = (μ = T(initial_conditions(qᵢ, vᵢ).params.μ),)


    initial_conditions_barely_passing() = initial_conditions(from_cartesian(0, [2.5, 0., 0.]), 3.425E-1, 1E-2)
    initial_conditions_barely_trapped() = initial_conditions(from_cartesian(0, [2.5, 0., 0.]), 3.375E-1, 1E-2)
    initial_conditions_deeply_passing() = initial_conditions(from_cartesian(0, [2.5, 0., 0.]),  5E-1,    1E-2)
    initial_conditions_deeply_trapped() = initial_conditions(from_cartesian(0, [2.5, 0., 0.]),  1E-1,    1E-2)
    initial_conditions_trapped()        = initial_conditions(from_cartesian(0, [7.0, 0., 0.]), -2E-3,    1.88E-7)

end
