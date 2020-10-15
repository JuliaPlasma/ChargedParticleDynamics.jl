"""
Charged Particle in a singular magnetic field of the form
``B(x,y,z) = (x^2 + y^2)^{-3/2} e_z``.
"""
module SingularField

    import ElectromagneticFields.Singular

    export charged_particle_3d_iode, hamiltonian, angular_momentum

    Singular.@code() # inject magnetic field code

    const q₀ = [1., 0., 0., 0., -1., 0.]

    include("charged_particle_3d.jl")

end
