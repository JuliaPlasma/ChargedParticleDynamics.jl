@doc raw"""
The charged particle in a static electromagnetic field: the full Lorentz-force dynamics in three
dimensions, with no averaging over the gyration.

Available in a canonical formulation on the phasespace ``(x, p)`` and a noncanonical one on
``(x, v)``, one module per equilibrium. See the [Charged Particles in 3D](@ref) page.
"""
module ChargedParticle3d

# The equations of the two formulations, written once and reading the field from `params.field`.
# The equilibrium modules below each hold a field, an initial condition and the problem
# constructors that default to them.
"""
The canonical formulation of the charged particle, on the phasespace ``(x, p)``, in any field whose
chart is orthogonal.
"""
module Canonical
include("charged_particle_3d/charged_particle_3d_canonical.jl")
end

"""
The noncanonical formulation of the charged particle, on ``(x, v)``, in any field whose chart is
orthogonal.
"""
module Noncanonical
include("charged_particle_3d/charged_particle_3d_noncanonical.jl")
end

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
