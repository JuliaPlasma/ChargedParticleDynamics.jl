"""
Charged Particle in an axisymmetric magnetic field of the form
``B(x,y,z) = (1 + x^2 + y^2) e_z``.
"""
module PauliParticle3dSymmetric

    export pauli_particle_3d_pode, hamiltonian, angular_momentum

    const E₀ = 2
    const qᵢ = [1., 0., 0.]
    const vᵢ = [0., 1., 1.]
    const parameters = (μ = 1E-2,)

    include("../electromagneticfields/symmetric.jl")
    include("pauli_particle_3d.jl")

end
