"""
Analytic ITER-like Solov'ev equilibrium with X-point.
"""
module SolovevIterXpoint

    import ElectromagneticFields.Solovev

    export initial_conditions_barely_passing, initial_conditions_barely_trapped,
           initial_conditions_deeply_passing, initial_conditions_deeply_trapped,
           initial_conditions_trapped

    export hamiltonian, toroidal_momentum

    Solovev.@code_iter_xpoint() # inject magnetic field code

    const Δt = 1.0
    const tspan = (0.0, 1000.0)

    const xᵢ = [7.0-1.4, 0.0, 0.0]
    const qᵢ = [from_cartesian(0, xᵢ)..., 2.8166280889939737]
    const parameters = (μ = 4.607782183567846,)

    include("guiding_center_3d_equations.jl")

    function toroidal_momentum(t,q)
        R(t,q) * ϑ₃(t,q)
    end

    include("guiding_center_3d_diagnostics.jl")

    initial_conditions_barely_passing() = ([from_cartesian(0, 2.5, 0., 0.)..., 3.425E-1], (μ = 1E-2,))
    initial_conditions_barely_trapped() = ([from_cartesian(0, 2.5, 0., 0.)..., 3.375E-1], (μ = 1E-2,))
    initial_conditions_deeply_passing() = ([from_cartesian(0, 2.5, 0., 0.)...,  5E-1],    (μ = 1E-2,))
    initial_conditions_deeply_trapped() = ([from_cartesian(0, 2.5, 0., 0.)...,  1E-1],    (μ = 1E-2,))
    initial_conditions_trapped()        = ([from_cartesian(0, 7.0, 0., 0.)..., -2E-3],    (μ = 1.88E-7,))

end
