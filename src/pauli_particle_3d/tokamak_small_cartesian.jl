module TokamakSmallCartesian

    import ElectromagneticFields.AxisymmetricTokamakCartesian

    export pauli_particle_3d_pode, hamiltonian, toroidal_momentum

    const equ = AxisymmetricTokamakCartesian.@code() # inject magnetic field code

    const qᵢ = [1.05,   0.0,    0.0]
    const vᵢ = [2.1E-3, 4.3E-4, 0.0]

    const Δt = 400.0
    const tspan = (0.0, 2E4)

    include("pauli_particle_3d.jl")

end
