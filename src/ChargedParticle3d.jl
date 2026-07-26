@doc raw"""
The charged particle in a static electromagnetic field: the full Lorentz-force dynamics in three
dimensions, with no averaging over the gyration.

Available in a canonical formulation on the phasespace ``(x, p)`` and a noncanonical one on
``(x, v)``, one module per equilibrium. See the [Charged Particles in 3D](@ref) page.
"""
module ChargedParticle3d

    include("charged_particle_3d/singular_field_canonical.jl")
    include("charged_particle_3d/singular_field.jl")
    include("charged_particle_3d/symmetric_field.jl")
    include("charged_particle_3d/solovev_iter.jl")
    include("charged_particle_3d/solovev_iter_xpoint.jl")
    include("charged_particle_3d/theta_pinch_canonical.jl")
    include("charged_particle_3d/theta_pinch_noncanonical.jl")
    include("charged_particle_3d/tokamak_iter_cylindrical.jl")
    include("charged_particle_3d/tokamak_small_noncanonical.jl")
    include("charged_particle_3d/tokamak_small_cartesian.jl")
    include("charged_particle_3d/tokamak_small_cylindrical.jl")
    include("charged_particle_3d/tokamak_small_toroidal.jl")

end
