module TokamakSmallCylindrical

    import ElectromagneticFields: load_equilibrium, AxisymmetricTokamakCylindrical

    export charged_particle_3d_pode, charged_particle_3d_iode,
           hamiltonian, toroidal_momentum

    const R₀ = 1.
    const B₀ = 1.
    const q  = 2.

    load_equilibrium(AxisymmetricTokamakCylindrical(R₀, B₀, q); target_module=TokamakSmallCylindrical)

    include("charged_particle_3d_canonical.jl")

    const qᵢ = from_cartesian(0, [1.05,   0.0,    0.0])
    const vᵢ = from_cartesian(0, [2.1E-3, 4.3E-4, 0.0])
    const pᵢ = charged_particle_3d_pᵢ(qᵢ, vᵢ)

end
