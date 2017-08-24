__precompile__()

module ChargedParticle3d

    using GeometricIntegrators

    export ChargedParticle3dSingular, ChargedParticle3dSymmetric, ChargedParticle3dUniform

    include("charged_particle_3d/charged_particle_3d_singular.jl")
    include("charged_particle_3d/charged_particle_3d_symmetric.jl")
    include("charged_particle_3d/charged_particle_3d_uniform.jl")


end
