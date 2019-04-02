__precompile__()

module GuidingCenter4d

    using GeometricIntegrators
    using ElectromagneticFields
    using LinearAlgebra


    export GuidingCenter4dSolovevQuadratic,
           GuidingCenter4dSymmetricQuadratic,
           GuidingCenter4dThetaPinch

    include("guiding_center_4d/solovev_quadratic.jl")
    include("guiding_center_4d/symmetric_quadratic.jl")
    include("guiding_center_4d/theta_pinch.jl")


    export TokamakMediumCartesian,
           TokamakMediumCylindrical,
           TokamakSmallCylindrical

    include("guiding_center_4d/tokamak_medium_cartesian.jl")
    include("guiding_center_4d/tokamak_medium_cylindrical.jl")
    include("guiding_center_4d/tokamak_small_cylindrical.jl")


    export convert_coordinates_RZphi_to_xyz,
           compute_energy, compute_toroidal_momentum, compute_one_form,
           compute_momentum_error, compute_error_drift

    include("guiding_center_4d/guiding_center_4d_diagnostics.jl")

end
