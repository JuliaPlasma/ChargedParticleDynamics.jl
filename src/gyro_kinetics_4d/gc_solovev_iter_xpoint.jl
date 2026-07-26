"""
Analytic ITER-like Solov'ev equilibrium with X-point.
"""
module GuidingCenter4dSolovevIterXpoint

    using ElectromagneticFields.Solovev

    export initial_conditions_barely_passing, initial_conditions_barely_trapped,
           initial_conditions_deeply_passing, initial_conditions_deeply_trapped,
           initial_conditions_trapped

    # export toroidal_momentum

    Solovev.@code_iter_xpoint() # inject magnetic field code

    include("coordinate_transformations.jl")
    include("gc_common.jl")
    include("gc_equations.jl")

    # The state is (R, u) directly. It used to be handed to the problem constructors scaled by
    # ω₀ = B*∥(t₀, q₀), a literal reading of the "R̃ = B*∥ R" of the notes — but that shorthand
    # denotes the time reparametrisation ds = B*∥ dt, which the vector field of `gc_common.jl`
    # already carries, not a change of variables. Scaling the state on top of it, with the
    # Jacobian left off the vector field as well, integrated a different system entirely.
    function solovev_xpoint_iter_initial_conditions(t₀, q₀, μ)
        local params = (μ=μ, R₀=R₀, ω₀=ωabs(t₀, q₀))
        (q = q₀, params = params)
    end

    # In the rescaled time s, with dt = B*∥ ds and B*∥ ≈ 2·10² here, this span covers a few
    # hundred units of physical time — comparable to the defaults of the 4D guiding centre model,
    # whose `Δt = 1.0` would be a factor of B*∥ too large for this system.
    const DEFAULT_TIMESTEP = 1E-3
    const DEFAULT_TIMESPAN = (0.0, 1.0)

    const x₀ = from_cartesian(0, [2.5, 0., 0.])
    const μ₀ = 1E-2

    export default_parameters

    """
    The parameters of the default initial condition: the magnetic moment `μ`, the major radius
    `R₀` of the equilibrium, and `ω₀ = B*∥` evaluated at the deeply passing initial condition.

    `ω₀` is not used by the equations of motion — it is the scale factor of the coordinate
    transformation in `coordinate_transformations.jl`, which is retained as a utility but not
    applied; see the module documentation.
    """
    default_parameters(::Type{T}=Float64) where {T} =
        (μ = T(μ₀), R₀ = T(R₀), ω₀ = T(ωabs(0, [x₀..., 5.0E-1])))

    const parameters = default_parameters()

    initial_conditions_barely_passing() = solovev_xpoint_iter_initial_conditions(0, [x₀..., 3.425E-1], μ₀)
    initial_conditions_barely_trapped() = solovev_xpoint_iter_initial_conditions(0, [x₀..., 3.375E-1], μ₀)
    initial_conditions_deeply_passing() = solovev_xpoint_iter_initial_conditions(0, [x₀..., 5.0E-1  ], μ₀)
    initial_conditions_deeply_trapped() = solovev_xpoint_iter_initial_conditions(0, [x₀..., 1.0E-1  ], μ₀)

    initial_conditions_trapped() = solovev_xpoint_iter_initial_conditions(0, [from_cartesian(0, [7.0, 0., 0.])..., -2E-3], 1.88E-7)

    include("../guiding_center_4d/guiding_center_4d_diagnostics.jl")

end
