module ChargedParticle3d

    include("charged_particle_3d/charged_particle_3d_singular.jl")
    include("charged_particle_3d/charged_particle_3d_symmetric.jl")
    include("charged_particle_3d/charged_particle_3d_uniform.jl")

    include("charged_particle_3d/solovev_iter_xpoint.jl")
    include("charged_particle_3d/theta_pinch_noncanonical.jl")
    include("charged_particle_3d/tokamak_small_noncanonical.jl")
    include("charged_particle_3d/tokamak_small_cartesian.jl")
    include("charged_particle_3d/tokamak_small_cylindrical.jl")
    include("charged_particle_3d/tokamak_small_toroidal.jl")

end
