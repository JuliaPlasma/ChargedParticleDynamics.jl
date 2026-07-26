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
    #   * `SolovevSymmetricField` cannot host this model at all: `SolovevSymmetric.@code` injects
    #     the equilibrium parameters `α` and `β` into the module as constants, and `β` collides
    #     with the vector potential `β` of `gc_common.jl` — "cannot define function β; it already
    #     has a value". Renaming the model's `β` would diverge from the notes it follows and change
    #     its exported interface, so the equilibrium is left out instead.
    include("gyro_kinetics_4d/gc_solovev_iter.jl")
    include("gyro_kinetics_4d/gc_solovev_iter_xpoint.jl")
    include("gyro_kinetics_4d/gc_tokamak_iter_cylindrical.jl")
    include("gyro_kinetics_4d/gc_tokamak_medium_cartesian.jl")
    include("gyro_kinetics_4d/gc_tokamak_medium_cylindrical.jl")
    include("gyro_kinetics_4d/gc_tokamak_small_cartesian.jl")
    include("gyro_kinetics_4d/gc_tokamak_small_cylindrical.jl")
    include("gyro_kinetics_4d/gc_tokamak_small_toroidal.jl")

end
