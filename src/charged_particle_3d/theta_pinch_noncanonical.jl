module ThetaPinchNoncanonical

    using GeometricIntegrators.Equations

    using ElectromagneticFields.ThetaPinch

    export charged_particle_3d_ode, charged_particle_3d_sode, charged_particle_3d_iode,
           hamiltonian, toroidal_momentum, ϑ

    const B₀ = 2.

    const qᵢ = [2.5, 0.0, 0.0, 0.0, 0.2, 0.1]

    const equ = ThetaPinch.init(B₀)

    include("charged_particle_3d_noncanonical.jl")
    
end
