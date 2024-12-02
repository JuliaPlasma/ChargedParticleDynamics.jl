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
    const parameters = (μ = 1E-2,)

    const Δt = 1.0
    const tspan = (0.0, 1000.0)

    include("pauli_particle_3d.jl")

end
