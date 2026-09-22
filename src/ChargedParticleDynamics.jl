module ChargedParticleDynamics

include("utils/field_points.jl")
include("utils/periodicity.jl")
include("utils/coordinates.jl")
include("utils/initial_conditions.jl")

include("ChargedParticle3d.jl")
# [DEBUG-cpdport] families not yet ported
# include("GuidingCenter3d.jl")
include("GuidingCenter4d.jl")
# include("GyroKinetics4d.jl")
include("PauliParticle3d.jl")

include("Plots.jl")

end
