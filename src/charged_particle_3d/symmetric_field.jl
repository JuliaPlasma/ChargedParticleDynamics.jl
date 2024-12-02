"""
Charged Particle in an axisymmetric magnetic field of the form
``B(x,y,z) = (1 + x^2 + y^2) e_z``.
"""
module SymmetricField

    import ElectromagneticFields.SymmetricQuadratic

    export charged_particle_3d_ode, charged_particle_3d_iode
    export hamiltonian, angular_momentum
    export compute_energy, compute_energy_error

    SymmetricQuadratic.@code() # inject magnetic field code

    const qᵢ = [1., 0., 0., 0., 1., 1.]

    angular_momentum(t,q) = q[1] * ϑ₂(t,q) - q[2] * ϑ₁(t,q)

    include("charged_particle_3d_noncanonical.jl")

end
