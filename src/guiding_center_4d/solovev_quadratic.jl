"""
Analytic, quadratic Solov'ev equilibrium.
"""
module SolovevQuadraticField

    using ElectromagneticFields.SolovevQuadratic

    export initial_conditions_barely_passing, initial_conditions_barely_trapped,
           initial_conditions_deeply_passing, initial_conditions_deeply_trapped

    export hamiltonian, toroidal_momentum

    const equ = SolovevQuadratic.init(2., 5., 1., 1.)

    initial_conditions_barely_passing() = ([2.5, 0., 0., 3.425E-1], (μ = 1E-2,)) # Δt=2.5, nt=50
    initial_conditions_barely_trapped() = ([2.5, 0., 0., 3.375E-1], (μ = 1E-2,)) # Δt=3.0, nt=100
    initial_conditions_deeply_passing() = ([2.5, 0., 0., 5E-1], (μ = 1E-2,))     # Δt=2.5, nt=25
    initial_conditions_deeply_trapped() = ([2.5, 0., 0., 1E-1], (μ = 1E-2,))     # Δt=5.0, nt=50

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")

    function toroidal_momentum(t,q)
        R(t,q) * ϑ3(t,q)
    end

    include("guiding_center_4d_diagnostics.jl")

end
