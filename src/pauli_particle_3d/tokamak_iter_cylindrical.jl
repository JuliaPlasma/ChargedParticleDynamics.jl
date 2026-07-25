module TokamakIterCylindrical

    import ElectromagneticFields.AxisymmetricTokamakCylindrical

    export pauli_particle_3d_pode, hamiltonian, toroidal_momentum

    AxisymmetricTokamakCylindrical.@code_iter() # inject magnetic field code

    const qᵢ = [7.0-1.4, 0.0, 0.0]
    const vᵢ = [3.43E-3, 6.75, -3.41E-1]

    include("pauli_particle_3d.jl")

    export default_parameters

    """
    The magnetic moment μ of the default initial condition `(qᵢ, vᵢ)`, obtained by splitting
    `vᵢ` into its parallel and perpendicular parts at `qᵢ`.
    """
    default_parameters(::Type{T}=Float64) where {T} = (μ = T(initial_conditions(qᵢ, vᵢ)[3].μ),)

    const parameters = default_parameters()

    const Δt = 1.0
    const tspan = (0.0, 100.0)

    function initial_conditions(x₀, u₀, μ)
        v₀ = u₀ * [b¹(0, x₀), b²(0, x₀), b³(0, x₀)]
        (x₀, v₀, (μ = μ,))
    end

    initial_conditions_barely_trapped() = initial_conditions(from_cartesian(0, [2.5, 0., 0.]), 3.375E-1, 1E-2)
    initial_conditions_deeply_passing() = initial_conditions(from_cartesian(0, [2.5, 0., 0.]),  5E-1,    1E-2)
    initial_conditions_deeply_trapped() = initial_conditions(from_cartesian(0, [2.5, 0., 0.]),  1E-1,    1E-2)
    initial_conditions_trapped()        = initial_conditions(from_cartesian(0, [7.0, 0., 0.]), -2E-3,    1.88E-7)

end
