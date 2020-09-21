
module ChargedParticle3dThetaPinchNoncanonical

    using GeometricIntegrators.Equations

    import ElectromagneticFields

    export charged_particle_3d_ode, charged_particle_3d_sode, charged_particle_3d_iode,
           hamiltonian, toroidal_momentum, ϑ

    const R0 = 2.
    
    const qᵢ = [2.5, 0., 0., 0.0, 0.2, 0.1]

    ElectromagneticFields.load_equilibrium(ElectromagneticFields.ThetaPinch(R0); target_module=ChargedParticle3dThetaPinchNoncanonical)

    include("charged_particle_3d_noncanonical.jl")
    
end
