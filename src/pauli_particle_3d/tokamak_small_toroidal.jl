module TokamakSmallToroidal

    import ElectromagneticFields.AxisymmetricTokamakToroidal

    export pauli_particle_3d_pode, hamiltonian, toroidal_momentum

    AxisymmetricTokamakToroidal.@code() # inject magnetic field code

    const xᵢ = [1.05, 0.0, 0.0]
    const qᵢ = from_cartesian(0, xᵢ)
    const vᵢ = DF̄(0, qᵢ) * [2.1E-3, 4.3E-4, 0.0]

    include("pauli_particle_3d.jl")

end
