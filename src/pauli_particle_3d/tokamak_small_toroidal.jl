module TokamakSmallToroidal

    import ElectromagneticFields: load_equilibrium, AxisymmetricTokamakToroidal

    export pauli_particle_3d_pode, hamiltonian, toroidal_momentum

    const R0 = 1.
    const B0 = 1.
    const q  = 2.

    const qᵢ = [0.05, 0., 0.]
    const vᵢ = [2.1E-3, 0., -4.3E-4]

    load_equilibrium(AxisymmetricTokamakToroidal(R0, B0, q); target_module=TokamakSmallToroidal)

    include("pauli_particle_3d.jl")

end
