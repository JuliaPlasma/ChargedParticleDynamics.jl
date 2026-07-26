@doc raw"""
The Pauli particle in three dimensions: the charged particle plus the ``\mu \vert B \vert`` term of
the guiding centre, on the full six-dimensional phasespace.

It sits between the charged particle and the guiding centre — it keeps a *regular* Lagrangian, so it
needs no projection onto a constraint manifold. See the [Pauli Particles in 3D](@ref) page.
"""
module PauliParticle3d

    include("pauli_particle_3d/symmetric_field.jl")
    include("pauli_particle_3d/theta_pinch.jl")
    include("pauli_particle_3d/solovev_iter.jl")
    include("pauli_particle_3d/solovev_iter_xpoint.jl")
    include("pauli_particle_3d/tokamak_iter_cylindrical.jl")
    include("pauli_particle_3d/tokamak_small_cartesian.jl")
    include("pauli_particle_3d/tokamak_small_cylindrical.jl")
    include("pauli_particle_3d/tokamak_small_toroidal.jl")

end
