__precompile__()

module GuidingCenter4d

    using GeometricIntegrators

    export TokamakSlowBarelyPassing,
           TokamakSlowBarelyTrapped,
           TokamakSlowDeeplyPassing,
           TokamakSlowDeeplyTrapped,
           TokamakSlowLoop,
           TokamakSlowSurface,
           TokamakFastBarelyPassing,
           TokamakFastBarelyTrapped,
           TokamakFastDeeplyPassing,
           TokamakFastDeeplyTrapped,
           TokamakFastLoop,
           TokamakFastSurface,
           SymmetricLoop,
           SymmetricSurface,
           UniformLoop

    export convert_coordinates_RZphi_to_xyz,
           compute_energy, compute_toroidal_momentum, compute_one_form,
           compute_momentum_error, compute_error_drift


    include("guiding_center_4d/guiding_center_4d_tokamak_slow_barely_passing.jl")
    include("guiding_center_4d/guiding_center_4d_tokamak_slow_barely_trapped.jl")
    include("guiding_center_4d/guiding_center_4d_tokamak_slow_deeply_passing.jl")
    include("guiding_center_4d/guiding_center_4d_tokamak_slow_deeply_trapped.jl")
    include("guiding_center_4d/guiding_center_4d_tokamak_slow_loop.jl")
    include("guiding_center_4d/guiding_center_4d_tokamak_slow_surface.jl")

    include("guiding_center_4d/guiding_center_4d_tokamak_fast_barely_passing.jl")
    include("guiding_center_4d/guiding_center_4d_tokamak_fast_barely_trapped.jl")
    include("guiding_center_4d/guiding_center_4d_tokamak_fast_deeply_passing.jl")
    include("guiding_center_4d/guiding_center_4d_tokamak_fast_deeply_trapped.jl")
    include("guiding_center_4d/guiding_center_4d_tokamak_fast_loop.jl")
    include("guiding_center_4d/guiding_center_4d_tokamak_fast_surface.jl")

    include("guiding_center_4d/guiding_center_4d_symmetric_loop.jl")
    include("guiding_center_4d/guiding_center_4d_symmetric_surface.jl")

    include("guiding_center_4d/guiding_center_4d_uniform_loop.jl")

    include("guiding_center_4d/guiding_center_4d_diagnostics.jl")

end
