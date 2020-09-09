module ChargedParticle3d

    include("charged_particle_3d/charged_particle_3d_singular.jl")
    include("charged_particle_3d/charged_particle_3d_symmetric.jl")
    include("charged_particle_3d/charged_particle_3d_uniform.jl")

    include("charged_particle_3d/charged_particle_3d_theta_pinch_noncanonical.jl")
    include("charged_particle_3d/charged_particle_3d_tokamak_canonical.jl")
    include("charged_particle_3d/charged_particle_3d_tokamak_noncanonical.jl")

end
