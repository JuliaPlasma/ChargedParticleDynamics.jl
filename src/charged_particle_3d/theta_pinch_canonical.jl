"""
Charged Particle in an uniform magnetic field of the form
``B(x,y,z) = B_0 e_z``.
"""
module ThetaPinchCanonical

    import ElectromagneticFields.ThetaPinch

    export charged_particle_3d_iode, hamiltonian, angular_momentum

    ThetaPinch.@code() # inject magnetic field code

    include("charged_particle_3d_canonical.jl")

    const qᵢ = [2.5, 0.0, 0.0]
    const vᵢ = [0.0, 0.2, 0.1]
    const pᵢ = charged_particle_3d_pᵢ(tᵢ, qᵢ, vᵢ)

end
