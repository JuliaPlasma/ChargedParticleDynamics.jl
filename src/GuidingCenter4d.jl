__precompile__()

module GuidingCenter4d

    using GeometricIntegrators
    using ElectromagneticFields
    using LinearAlgebra


    export SymmetricLoop,
           SymmetricSurface,
           UniformLoop

    include("guiding_center_4d/guiding_center_4d_symmetric_loop.jl")
    include("guiding_center_4d/guiding_center_4d_symmetric_surface.jl")

    include("guiding_center_4d/guiding_center_4d_uniform_loop.jl")


    export TokamakSlowBarelyPassing,
           TokamakSlowBarelyTrapped,
           TokamakSlowDeeplyPassing,
           TokamakSlowDeeplyTrapped,
           TokamakSlowLoop,
           TokamakSlowSurface

    include("guiding_center_4d/guiding_center_4d_tokamak_cylindrical_slow.jl")


    export TokamakFastBarelyPassing,
           TokamakFastBarelyTrapped,
           TokamakFastDeeplyPassing,
           TokamakFastDeeplyTrapped,
           TokamakFastLoop,
           TokamakFastSurface

    include("guiding_center_4d/guiding_center_4d_tokamak_cylindrical_fast.jl")


    export convert_coordinates_RZphi_to_xyz,
           compute_energy, compute_toroidal_momentum, compute_one_form,
           compute_momentum_error, compute_error_drift

    include("guiding_center_4d/guiding_center_4d_diagnostics.jl")

end
