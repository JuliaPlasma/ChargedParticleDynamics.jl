module ChargedParticleDynamics

    export InitialConditions, charged_particle, guiding_center, pauli_particle

    include("initialization/initial_conditions.jl")
    
    include("ChargedParticle3d.jl")
    include("GuidingCenter4d.jl")
    include("GyroKinetics4d.jl")
    include("PauliParticle3d.jl")
    include("ChargedParticlePlots.jl")

end
