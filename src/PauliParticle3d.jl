module PauliParticle3d

    include("pauli_particle_3d/pauli_particle_3d_symmetric.jl")
    include("pauli_particle_3d/pauli_particle_3d_uniform.jl")
    include("pauli_particle_3d/solovev_iter_xpoint.jl")
    include("pauli_particle_3d/tokamak_iter_cylindrical.jl")
    include("pauli_particle_3d/tokamak_small_cartesian.jl")
    include("pauli_particle_3d/tokamak_small_cylindrical.jl")
    include("pauli_particle_3d/tokamak_small_toroidal.jl")

end
