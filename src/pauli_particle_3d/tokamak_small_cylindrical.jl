module TokamakSmallCylindrical

    using ElectromagneticFields.AxisymmetricTokamakCylindrical

    export pauli_particle_3d_pode, hamiltonian, toroidal_momentum

    const R₀ = 1.
    const B₀ = 1.
    const E₀ = 0.
    const q  = 2.

    const equ = AxisymmetricTokamakCylindrical.init(R₀, B₀, q)

    const qᵢ = [1.05, 0., 0.]
    const vᵢ = [2.1E-3, 0., -4.3E-4]

    include("pauli_particle_3d.jl")

end
