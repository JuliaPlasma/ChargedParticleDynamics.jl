
module TokamakSmallCylindrical

    import ElectromagneticFields: load_equilibrium, TokamakSmallCylindrical

    export pauli_particle_3d_pode, hamiltonian, toroidal_momentum

    const R0 = 1.
    const B0 = 1.
    const q  = 2.

    const qᵢ = [1.05, 0., 0.]
    const vᵢ = [2.1E-3, 0., -4.3E-4]

    load_equilibrium(AxisymmetricTokamakCylindrical(R0, B0, q); target_module=TokamakSmallCylindrical)

    include("pauli_particle_3d.jl")

end
