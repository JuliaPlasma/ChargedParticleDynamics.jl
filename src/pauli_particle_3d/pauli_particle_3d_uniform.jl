"""
Charged Particle in an uniform magnetic field of the form
``B(x,y,z) = B_0 e_z``.
"""
module PauliParticle3dUniform

    export pauli_particle_3d_pode, hamiltonian, angular_momentum

    const E₀ = 2
    const qᵢ = [1., 0., 0.]
    const vᵢ = [0., 1., 1.]
    const parameters = (μ = 1E-2,)

    include("../electromagneticfields/uniform.jl")
    include("pauli_particle_3d.jl")

end
