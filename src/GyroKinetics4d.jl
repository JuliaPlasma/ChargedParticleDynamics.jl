@doc raw"""
The characteristics of the gyrokinetic Vlasov equation, in the rescaled time of the
volume-preserving formulation with ``dt = B^{\star}_{\parallel} \, ds``.

The same physics as `GuidingCenter4d`, written so that the right-hand side is divergence-free and
admits an exactly volume-preserving splitting. See the
[Gyrokinetic Guiding Centre Dynamics in 4D](@ref) page.
"""
module GyroKinetics4d

# The gyrokinetic guiding centre model, one module per equilibrium, mirroring `GuidingCenter4d`.
#
# Three of that module's eleven equilibria are absent, for two reasons:
#
#   * `SymmetricField` and `ThetaPinchField` carry only a Poincaré loop and surface
#     parameterisation and no point initial condition, and this model has no loop/surface
#     machinery, so there is nothing for them to integrate here.
#   * `SolovevSymmetricField` has no module here. The field is a value in `params` and no
#     longer defines the equilibrium parameters `α` and `β` as constants of the module, so its
#     `β` no longer collides with the vector potential `β` of `gc_common.jl`; the module has
#     simply not been written.
#
# The equations, written once and reading the field from `params.field`. The equilibrium modules
# below each hold a field, an initial condition and the problem constructors that default to them.
include("gyro_kinetics_4d/coordinate_transformations.jl")
include("gyro_kinetics_4d/gc_common.jl")
include("gyro_kinetics_4d/gc_equations.jl")

include("gyro_kinetics_4d/gc_solovev_iter.jl")
include("gyro_kinetics_4d/gc_solovev_iter_xpoint.jl")
include("gyro_kinetics_4d/gc_tokamak_iter_cylindrical.jl")
include("gyro_kinetics_4d/gc_tokamak_medium_cartesian.jl")
include("gyro_kinetics_4d/gc_tokamak_medium_cylindrical.jl")
include("gyro_kinetics_4d/gc_tokamak_small_cartesian.jl")
include("gyro_kinetics_4d/gc_tokamak_small_cylindrical.jl")
include("gyro_kinetics_4d/gc_tokamak_small_toroidal.jl")

end
