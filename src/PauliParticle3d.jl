module PauliParticle3d

    include("pauli_particle_3d/symmetric_field.jl")
    include("pauli_particle_3d/theta_pinch.jl")
    include("pauli_particle_3d/solovev_iter_xpoint.jl")
    include("pauli_particle_3d/tokamak_iter_cylindrical.jl")
    include("pauli_particle_3d/tokamak_small_cartesian.jl")
    include("pauli_particle_3d/tokamak_small_circular.jl")
    include("pauli_particle_3d/tokamak_small_cylindrical.jl")

end
