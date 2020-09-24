"""
Charged Particle in an uniform magnetic field of the form
``B(x,y,z) = B_0 e_z``.
"""
module ChargedParticle3dUniform

    using ElectromagneticFields.ThetaPinch

    export charged_particle_3d_iode, hamiltonian, angular_momentum

    const B₀ = 2.
    const E₀ = 2.

    const equ = ThetaPinch.init(B₀)

    const q₀ = [2.5, 0.0, 0.0, 0.0, 0.2, 0.1]

    include("charged_particle_3d.jl")

end
