
module ChargedParticle3dTokamakCanonical

    using GeometricIntegrators.Equations

    import ElectromagneticFields

    export charged_particle_3d_pode, hamiltonian, toroidal_momentum

    const R0 = 2.
    const B0 = 5.
    const q  = 2.

    const qᵢ = [1.05, 0., 0.]
    const vᵢ = [0.0, 0.2, 0.1]

    ElectromagneticFields.load_equilibrium(ElectromagneticFields.AxisymmetricTokamakCylindrical(R0, B0, q); target_module=ChargedParticle3dTokamakCanonical)

    include("charged_particle_3d_canonical.jl")

end
