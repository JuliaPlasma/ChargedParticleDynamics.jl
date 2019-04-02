__precompile__()

module GuidingCenter4d

    using GeometricIntegrators
    using ElectromagneticFields
    using LinearAlgebra


    export SymmetricLoop,
           SymmetricSurface,
           UniformLoop

    include("guiding_center_4d/guiding_center_4d_symmetric.jl")
    include("guiding_center_4d/guiding_center_4d_uniform.jl")





    export TokamakCartesianFastBarelyPassing,
           TokamakCartesianFastBarelyTrapped,
           TokamakCartesianFastDeeplyPassing,
           TokamakCartesianFastDeeplyTrapped,
           TokamakCartesianFastLoop,
           TokamakCartesianFastSurface

    include("guiding_center_4d/guiding_center_4d_tokamak_cartesian_fast.jl")


    export TokamakCylindricalSlowBarelyPassing,
           TokamakCylindricalSlowBarelyTrapped,
           TokamakCylindricalSlowDeeplyPassing,
           TokamakCylindricalSlowDeeplyTrapped,
           TokamakCylindricalSlowLoop,
           TokamakCylindricalSlowSurface

    include("guiding_center_4d/guiding_center_4d_tokamak_cylindrical_slow.jl")


    export TokamakCylindricalFastBarelyPassing,
           TokamakCylindricalFastBarelyTrapped,
           TokamakCylindricalFastDeeplyPassing,
           TokamakCylindricalFastDeeplyTrapped,
           TokamakCylindricalFastLoop,
           TokamakCylindricalFastSurface

    include("guiding_center_4d/guiding_center_4d_tokamak_cylindrical_fast.jl")


    export convert_coordinates_RZphi_to_xyz,
           compute_energy, compute_toroidal_momentum, compute_one_form,
           compute_momentum_error, compute_error_drift

    include("guiding_center_4d/guiding_center_4d_diagnostics.jl")

end
