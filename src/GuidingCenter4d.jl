@doc raw"""
The guiding centre dynamics in four dimensions: the gyro-averaged motion in the guiding centre
position and the parallel velocity, ``q = (X, u)``, with the magnetic moment ``\mu`` a parameter.

A noncanonical, degenerate system built from the one-form ``\vartheta = A + u b``. See the
[Guiding Center Dynamics in 4D](@ref) page.
"""
module GuidingCenter4d

# The equations, written once and reading the field from `params.field`. The equilibrium modules
# below each hold a field, an initial condition and the problem constructors that default to them,
# together with the diagnostics and Poincaré invariants that read their own `toroidal_momentum`,
# `f_loop` and `f_surface`.
include("guiding_center_4d/guiding_center_4d_common.jl")
include("guiding_center_4d/guiding_center_4d_equations.jl")

include("guiding_center_4d/solovev_iter.jl")
include("guiding_center_4d/solovev_iter_xpoint.jl")
include("guiding_center_4d/solovev_symmetric.jl")
include("guiding_center_4d/symmetric_quadratic.jl")
include("guiding_center_4d/theta_pinch.jl")

include("guiding_center_4d/tokamak_iter_cylindrical.jl")
include("guiding_center_4d/tokamak_medium_cartesian.jl")
include("guiding_center_4d/tokamak_medium_cylindrical.jl")
include("guiding_center_4d/tokamak_small_cartesian.jl")
include("guiding_center_4d/tokamak_small_cylindrical.jl")
include("guiding_center_4d/tokamak_small_toroidal.jl")

end
