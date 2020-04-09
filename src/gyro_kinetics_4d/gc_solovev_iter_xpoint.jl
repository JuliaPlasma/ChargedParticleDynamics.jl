"""
Analytic ITER-like Solov'ev equilibrium with X-point.
"""
module GuidingCenter4dSolovevIterXpoint

    using ElectromagneticFields: load_equilibrium, SolovevXpointITER

    export initial_conditions_barely_passing, initial_conditions_barely_trapped,
           initial_conditions_deeply_passing, initial_conditions_deeply_trapped,
           initial_conditions_trapped


    equ = SolovevXpointITER()
    load_equilibrium(equ; target_module=GuidingCenter4dSolovevIterXpoint)

    include("coordinate_transformations.jl")
    include("gc_common.jl")
    include("gc_equations.jl")

    function solovev_xpoint_iter_initial_conditions(t₀, q₀, μ)
        local ω₀ = ωabs(t₀, q₀)
        local params = (μ=μ, x₀=[1., 0., 0.], R₀=equ.R₀, scaling_factor=[ω₀ * equ.R₀, ω₀ * equ.R₀, ω₀ , ω₀])
        (transform_q_to_q̃(t₀, q₀, params), params)
    end

    const x₀ = [2.5/equ.R₀, 0., 0.]
    const μ₀ = 1E-2

    initial_conditions_barely_passing() = solovev_xpoint_iter_initial_conditions(0., [x₀..., 3.425E-1], μ₀)
    initial_conditions_barely_trapped() = solovev_xpoint_iter_initial_conditions(0., [x₀..., 3.375E-1], μ₀)
    initial_conditions_deeply_passing() = solovev_xpoint_iter_initial_conditions(0., [x₀..., 5.0E-1  ], μ₀)
    initial_conditions_deeply_trapped() = solovev_xpoint_iter_initial_conditions(0., [x₀..., 1.0E-1  ], μ₀)

    initial_conditions_trapped() = solovev_xpoint_iter_initial_conditions(0., [7.0 / equ.R₀, 0., 0., -2E-3], 1.88E-7)

    include("../guiding_center_4d/guiding_center_4d_diagnostics.jl")

end
