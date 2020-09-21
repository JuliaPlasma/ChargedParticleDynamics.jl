"""
Charged Particle in an axisymmetric magnetic field of the form
``B(x,y,z) = (1 + x^2 + y^2) e_z``.
"""
module ChargedParticle3dSymmetric

    export charged_particle_3d_iode, hamiltonian, angular_momentum

    const E₀ = 0#2
    const q₀ = [1., 0., 0., 0., 1., 1.]

    include("../electromagneticfields/symmetric.jl")
    include("charged_particle_3d.jl")

end
