__precompile__()

module ChargedParticleDynamics

    using GeometricIntegrators

    export ChargedParticle3d, GuidingCenter4d

    include("ChargedParticle3d.jl")
    include("GuidingCenter4d.jl")

end
