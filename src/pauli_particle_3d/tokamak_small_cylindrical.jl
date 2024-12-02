module TokamakSmallCylindrical

    import ElectromagneticFields.AxisymmetricTokamakCylindrical

    export pauli_particle_3d_pode, hamiltonian, toroidal_momentum

    AxisymmetricTokamakCylindrical.@code() # inject magnetic field code

    const xᵢ = [1.05, 0.0, 0.0]
    const qᵢ = from_cartesian(0, xᵢ)
    const vᵢ = DF̄(0, xᵢ) * [2.1E-3, 4.3E-4, 0.0]

    const Δt = 400.0
    const tspan = (0.0, 2E4)

    include("pauli_particle_3d.jl")

end
