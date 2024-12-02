"""
Analytic ITER-like Solov'ev equilibrium with X-point.
"""
module TokamakIterCylindrical

    import ElectromagneticFields.AxisymmetricTokamakCylindrical

    export initial_conditions_barely_passing, initial_conditions_barely_trapped,
           initial_conditions_deeply_passing, initial_conditions_deeply_trapped,
           initial_conditions_trapped

    export hamiltonian, toroidal_momentum

    AxisymmetricTokamakCylindrical.@code_iter() # inject magnetic field code

    const Δt = 1.0
    const tspan = (0.0, 1000.0)

    const qᵢ = [7.0-1.4, 0.0, 0.0, 2.8166280889939737]
    const parameters = (μ = 4.607782183567846,)

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")

    function toroidal_momentum(t,q)
        R(t,q) * ϑ₃(t,q)
    end

    include("guiding_center_4d_diagnostics.jl")

    initial_conditions_barely_passing() = ([2.5, 0., 0., 3.425E-1], (μ = 1E-2,))
    initial_conditions_barely_trapped() = ([2.5, 0., 0., 3.375E-1], (μ = 1E-2,))
    initial_conditions_deeply_passing() = ([2.5, 0., 0.,  5E-1],    (μ = 1E-2,))
    initial_conditions_deeply_trapped() = ([2.5, 0., 0.,  1E-1],    (μ = 1E-2,))
    initial_conditions_trapped()        = ([7.0, 0., 0., -2E-3],    (μ = 1.88E-7,))

end
