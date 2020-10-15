"""
Charged Particle in an uniform magnetic field of the form
``B(x,y,z) = B_0 e_z``.
"""
module ThetaPinchCanonical

    import ElectromagneticFields.ThetaPinch

    export charged_particle_3d_iode, hamiltonian, angular_momentum

    ThetaPinch.@code() # inject magnetic field code

    const q₀ = [2.5, 0.0, 0.0, 0.0, 0.2, 0.1]

    include("charged_particle_3d.jl")

end
