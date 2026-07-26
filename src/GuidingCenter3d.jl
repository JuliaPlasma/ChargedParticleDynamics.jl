@doc raw"""
The guiding centre dynamics in three dimensions: the gyro-averaged motion in the position and its
conjugate momentum, with the parallel velocity eliminated by the two constraints ``v \times b = 0``.

A constrained canonical system, available as `hodeproblem` and `hodeproblem_canonical`. See the
[Guiding Center Dynamics in 3D](@ref) page.
"""
module GuidingCenter3d

    # Based on
    #   Xinjie Li, Ruili Zhang, Jian Liu,
    #   Approximately symplectic Runge-Kutta methods for the guiding center dynamics,
    #   Journal of Computational Physics 563, 115069, 2026, doi:10.1016/j.jcp.2026.115069,
    # and its companion
    #   Ruili Zhang, Jian Liu, Tong Liu, Wenxiang Li, Xiaogang Wang,
    #   Canonical Hamiltonian guiding center theory and classical intrinsic magnetic moment,
    #   Frontiers of Physics 21(2), 026200, 2026, doi:10.15302/frontphys.2026.026200.

    include("guiding_center_3d/dipole.jl")
    include("guiding_center_3d/quadratic_potentials.jl")
    include("guiding_center_3d/solovev_iter.jl")
    include("guiding_center_3d/solovev_iter_xpoint.jl")
    include("guiding_center_3d/solovev_symmetric.jl")
    include("guiding_center_3d/symmetric_quadratic.jl")
    include("guiding_center_3d/theta_pinch.jl")

    include("guiding_center_3d/tokamak_iter_cylindrical.jl")
    include("guiding_center_3d/tokamak_medium_cartesian.jl")
    include("guiding_center_3d/tokamak_medium_cylindrical.jl")
    include("guiding_center_3d/tokamak_small_cartesian.jl")
    include("guiding_center_3d/tokamak_small_cylindrical.jl")
    include("guiding_center_3d/tokamak_small_toroidal.jl")

end
