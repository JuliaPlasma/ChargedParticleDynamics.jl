"""
Charged Particle in an axisymmetric magnetic field of the form
``B(x,y,z) = (1 + x^2 + y^2) e_z``.
"""
module SymmetricField

    import ElectromagneticFields.SymmetricQuadratic

    export charged_particle_3d_iode, hamiltonian, angular_momentum

    SymmetricQuadratic.@code() # inject magnetic field code

    const q₀ = [1., 0., 0., 0., 1., 1.]

    include("charged_particle_3d.jl")

end
