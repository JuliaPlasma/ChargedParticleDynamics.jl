"""
Analytic, quadratic Solov'ev equilibrium.
"""
module SolovevSymmetricField

    import ElectromagneticFields.SolovevSymmetric

    export initial_conditions_barely_passing, initial_conditions_barely_trapped,
           initial_conditions_deeply_passing, initial_conditions_deeply_trapped

    export hamiltonian, toroidal_momentum

    SolovevSymmetric.@code(2., 5., 1., 1.) # inject magnetic field code

    const Δt = 5E3
    const tspan = (0.0, 5E5)

    initial_conditions_barely_passing() = ([2.5, 0., 0., 3.425E-1], (μ = 1E-2,)) # Δt=2.5, nt=50
    initial_conditions_barely_trapped() = ([2.5, 0., 0., 3.375E-1], (μ = 1E-2,)) # Δt=3.0, nt=100
    initial_conditions_deeply_passing() = ([2.5, 0., 0., 5E-1], (μ = 1E-2,))     # Δt=2.5, nt=25
    initial_conditions_deeply_trapped() = ([2.5, 0., 0., 1E-1], (μ = 1E-2,))     # Δt=5.0, nt=50

    export default_parameters

    "The magnetic moment μ this equilibrium is set up for."
    default_parameters(::Type{T}=Float64) where {T} = (μ = T(1E-2),)

    const parameters = default_parameters()

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")

    # The canonical toroidal momentum is the covariant φ-component of the one-form, ϑ₃. It was
    # previously multiplied by R, which destroys the conservation: on the small tokamak the
    # relative variation over 10³ time units is 2e-13 for ϑ₃ and 3e-3 for R ϑ₃.
    function toroidal_momentum(t,q)
        ϑ₃(t,q)
    end

    include("guiding_center_4d_diagnostics.jl")

end
