module TokamakSmallCartesian

    using ElectromagneticFields.AxisymmetricTokamakCartesian

    export pauli_particle_3d_pode, hamiltonian, toroidal_momentum

    const R₀ = 1.
    const B₀ = 1.
    const E₀ = 0.
    const q  = 2.

    const equ = AxisymmetricTokamakCartesian.init(R₀, B₀, q)

    const qᵢ = [1.05,   0.0,    0.0]
    const vᵢ = [2.1E-3, 4.3E-4, 0.0]

    include("pauli_particle_3d.jl")

end
