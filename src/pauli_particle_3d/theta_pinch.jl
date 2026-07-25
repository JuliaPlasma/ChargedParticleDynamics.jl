"""
Charged Particle in an uniform magnetic field of the form
``B(x,y,z) = B_0 e_z``.
"""
module ThetaPinchField

    import ElectromagneticFields.ThetaPinch

    export pauli_particle_3d_pode, hamiltonian, angular_momentum

    ThetaPinch.@code() # inject magnetic field code

    const qᵢ = [1., 0., 0.]
    const vᵢ = [0., 1., 1.]

    const Δt = 10.0
    const tspan = (0.0, 1000.0)

    include("pauli_particle_3d.jl")

    export default_parameters

    """
    The magnetic moment μ of the default initial condition `(qᵢ, vᵢ)`, obtained by splitting
    `vᵢ` into its parallel and perpendicular parts at `qᵢ`.
    """
    default_parameters(::Type{T}=Float64) where {T} = (μ = T(initial_conditions(qᵢ, vᵢ)[3].μ),)

    const parameters = default_parameters()

end
