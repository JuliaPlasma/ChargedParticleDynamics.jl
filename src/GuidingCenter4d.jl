module GuidingCenter4d

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
