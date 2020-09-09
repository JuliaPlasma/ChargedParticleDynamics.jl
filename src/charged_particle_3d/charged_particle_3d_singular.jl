"""
Charged Particle in a singular magnetic field of the form
``B(x,y,z) = (x^2 + y^2)^{-3/2} e_z``.
"""
module ChargedParticle3dSingular

    export charged_particle_3d_iode, hamiltonian, angular_momentum

    const E₀ = 0
    const q₀ = [1., 0., 0., 0., -1., 0.]

    include("../electromagneticfields/singular.jl")
    include("charged_particle_3d.jl")

end
