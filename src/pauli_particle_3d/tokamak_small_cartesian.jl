module TokamakSmallCartesian

    import ElectromagneticFields.AxisymmetricTokamakCartesian

    export pauli_particle_3d_pode, hamiltonian, toroidal_momentum

    const equ = AxisymmetricTokamakCartesian.@code() # inject magnetic field code

    const qᵢ = [1.05,   0.0,    0.0]
    const vᵢ = [2.1E-3, 4.3E-4, 0.0]

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
