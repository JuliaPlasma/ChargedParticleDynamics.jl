__precompile__()

module ChargedParticleDynamics

    export ChargedParticle3d, GuidingCenter4d, ChargedParticlePlots

    include("ChargedParticle3d.jl")
    include("GuidingCenter4d.jl")
    include("ChargedParticlePlots.jl")

end
