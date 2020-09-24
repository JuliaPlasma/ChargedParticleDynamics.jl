"""
Charged Particle in an axisymmetric magnetic field of the form
``B(x,y,z) = (1 + x^2 + y^2) e_z``.
"""
module SymmetricField

    using ElectromagneticFields.SymmetricQuadratic

    export pauli_particle_3d_pode, hamiltonian, toroidal_momentum

    const B₀ = 1.
    const E₀ = 2.

    const equ = SymmetricQuadratic.init(B₀)

    const qᵢ = [1., 0., 0.]
    const vᵢ = [0., 1., 1.]
    const parameters = (μ = 1E-2,)

    include("pauli_particle_3d.jl")

end
