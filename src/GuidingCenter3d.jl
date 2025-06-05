module GuidingCenter3d

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
