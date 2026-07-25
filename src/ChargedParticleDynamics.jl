module ChargedParticleDynamics

    include("utils/coordinates.jl")
    include("utils/initial_conditions.jl")

    include("ChargedParticle3d.jl")
    include("GuidingCenter3d.jl")
    include("GuidingCenter4d.jl")
    include("GyroKinetics4d.jl")
    include("PauliParticle3d.jl")

    include("Plots.jl")

end
