module ChargedParticleDynamics

    include("initialization/initial_conditions.jl")

    include("ChargedParticle3d.jl")
    include("GuidingCenter4d.jl")
    include("GyroKinetics4d.jl")
    include("PauliParticle3d.jl")
    include("ChargedParticlePlots.jl")

end
