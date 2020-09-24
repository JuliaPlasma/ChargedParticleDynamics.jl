module ChargedParticle3d

    include("charged_particle_3d/singular_field.jl")
    include("charged_particle_3d/symmetric_field.jl")
    include("charged_particle_3d/solovev_iter_xpoint.jl")
    include("charged_particle_3d/theta_pinch_canonical.jl")
    include("charged_particle_3d/theta_pinch_noncanonical.jl")
    include("charged_particle_3d/tokamak_small_noncanonical.jl")
    include("charged_particle_3d/tokamak_small_cartesian.jl")
    include("charged_particle_3d/tokamak_small_circular.jl")
    include("charged_particle_3d/tokamak_small_cylindrical.jl")

end
