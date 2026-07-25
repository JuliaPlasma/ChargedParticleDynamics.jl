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

    function solovev_xpoint_iter_initial_conditions(t₀, q₀, μ)
        local ω₀ = ωabs(t₀, q₀)
        local params = (μ=μ, R₀=R₀, ω₀=ω₀)
        (transform_q_to_q̃(t₀, q₀, params), params)
    end

    const Δt = 1.0
    const tspan = (0.0, 1000.0)

    const x₀ = from_cartesian(0, [2.5, 0., 0.])
    const μ₀ = 1E-2

    export default_parameters

    """
    The parameters of the default initial condition: the magnetic moment `μ`, the major radius
    `R₀` of the equilibrium, and the gyro frequency `ω₀` that sets the scaling of the transformed
    coordinates `q̃`. `ω₀` is evaluated at the deeply passing initial condition.
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
