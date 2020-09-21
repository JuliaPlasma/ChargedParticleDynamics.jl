module TokamakSmallNoncanonical

    import ElectromagneticFields: load_equilibrium, AxisymmetricTokamakToroidal

    export charged_particle_3d_ode, charged_particle_3d_sode, charged_particle_3d_iode,
           hamiltonian, toroidal_momentum, ϑ

    const R₀ = 1.
    const B₀ = 1.
    const q  = 2.

    load_equilibrium(AxisymmetricTokamakToroidal(R₀, B₀, q); target_module=TokamakSmallNoncanonical)
       
    const qᵢ = [1.05, 0., 0., 2.1E-3, 0.0, -4.3E-4]

    include("charged_particle_3d_noncanonical.jl")
    
end
