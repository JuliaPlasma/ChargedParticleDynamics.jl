module TokamakSmallCylindrical

    import ElectromagneticFields.AxisymmetricTokamakCylindrical

    export pauli_particle_3d_pode, hamiltonian, toroidal_momentum

    AxisymmetricTokamakCylindrical.@code() # inject magnetic field code

    const xᵢ = [1.05, 0.0, 0.0]
    const qᵢ = from_cartesian(0, xᵢ)
    const vᵢ = DF̄(0, qᵢ) * [2.1E-3, 4.3E-4, 0.0] # should DF be evaluated at qᵢ or xᵢ ?

    const Δt = 500.0
    const tspan = (0.0, 5E4)

    include("pauli_particle_3d.jl")

    export default_parameters

    """
    The magnetic moment μ of the default initial condition `(qᵢ, vᵢ)`, obtained by splitting
    `vᵢ` into its parallel and perpendicular parts at `qᵢ`.
    """
    default_parameters(::Type{T}=Float64) where {T} = (μ = T(initial_conditions(qᵢ, vᵢ)[3].μ),)

    const parameters = default_parameters()

end
