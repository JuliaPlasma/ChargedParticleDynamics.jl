module ChargedParticle3d

    using GeometricIntegrators
    using ElectromagneticFields
    using LinearAlgebra

    
    export ChargedParticle3dSingular,
           ChargedParticle3dSymmetric,
           ChargedParticle3dUniform
    
    include("charged_particle_3d/charged_particle_3d_singular.jl")
    include("charged_particle_3d/charged_particle_3d_symmetric.jl")
    include("charged_particle_3d/charged_particle_3d_uniform.jl")


    export ChargedParticle3dThetaPinchNoncanonical,
           ChargedParticle3dTokamakCanonical,
           ChargedParticle3dTokamakNoncanonical

    include("charged_particle_3d/charged_particle_3d_theta_pinch_noncanonical.jl")
    include("charged_particle_3d/charged_particle_3d_tokamak_canonical.jl")
    include("charged_particle_3d/charged_particle_3d_tokamak_noncanonical.jl")

end
