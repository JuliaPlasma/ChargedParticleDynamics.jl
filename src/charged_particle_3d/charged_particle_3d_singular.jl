"""
Charged Particle in a singular magnetic field of the form
``B(x,y,z) = (x^2 + y^2)^{-3/2} e_z``.
"""
module ChargedParticle3dSingular

    using ElectromagneticFields.Singular

    export charged_particle_3d_iode, hamiltonian, angular_momentum

    const B₀ = 1.
    const E₀ = 2.

    const equ = Singular.init(B₀)

    const q₀ = [1., 0., 0., 0., -1., 0.]

    include("charged_particle_3d.jl")

end
