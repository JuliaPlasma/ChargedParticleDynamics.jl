module PauliParticle3dSolovevIterXpoint

    using ElectromagneticFields: load_equilibrium, SolovevXpointITER, AxisymmetricTokamakCylindrical

    export pauli_particle_3d_pode, hamiltonian, toroidal_momentum

    const R₀ = 6.2
    const B₀ = 5.3
    const q  = 2.

    const qᵢ = [7.0-1.4, 0.0, 0.0]
    const vᵢ = [3.43E-3, 6.75, -3.41E-1]
    
    equ = AxisymmetricTokamakCylindrical(R₀, B₀, q)
    # equ = SolovevXpointITER()
    load_equilibrium(equ; target_module=PauliParticle3dSolovevIterXpoint)

    include("pauli_particle_3d.jl")

    function initial_conditions(x₀, u₀, μ)
        v₀ = u₀ * [b¹(0, x₀), b²(0, x₀), b³(0, x₀)]
        (x₀, v₀, (μ = μ,))
    end

    initial_conditions_barely_passing() = initial_conditions([2.5 / equ.R₀, 0., 0.], 3.425E-1, 1E-2)
    initial_conditions_barely_trapped() = initial_conditions([2.5 / equ.R₀, 0., 0.], 3.375E-1, 1E-2)
    initial_conditions_deeply_passing() = initial_conditions([2.5 / equ.R₀, 0., 0.],  5E-1,    1E-2)
    initial_conditions_deeply_trapped() = initial_conditions([2.5 / equ.R₀, 0., 0.],  1E-1,    1E-2)
    initial_conditions_trapped()        = initial_conditions([7.0 / equ.R₀, 0., 0.], -2E-3,    1.88E-7)

end
