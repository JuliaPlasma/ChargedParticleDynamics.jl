module ChargedParticle3dSolovevIterXpoint

    using ElectromagneticFields.SolovevXpoint

    export charged_particle_3d_pode, charged_particle_3d_iode,
           hamiltonian, toroidal_momentum

    const equ = SolovevXpoint.ITER()

    include("charged_particle_3d_canonical.jl")

    const qᵢ = [7.0, 0.0, 0.0]
    const vᵢ = [3.43E-3, 6.75, -3.41E-1]
    const pᵢ = charged_particle_3d_pᵢ(qᵢ, vᵢ)

    function initial_conditions(x₀, u₀, μ)
        vpar = u₀ * [b̂¹(0, x₀), b̂²(0, x₀), b̂³(0, x₀)]
        vper = sqrt(2 * μ * B(0, x₀))
        v¹ = vper * sqrt( b̂³(0, x₀)^2 / (b̂¹(0, x₀)^2 + b̂³(0, x₀)^2) )
        v³ = - v¹ * b̂¹(0, x₀) / b̂³(0, x₀)
        v₀ = vpar .+ [v¹, 0, v³]

        (x₀, charged_particle_3d_pᵢ(x₀, v₀))
    end

    initial_conditions_barely_passing() = initial_conditions([2.5 / equ.R₀, 0., 0.], 3.425E-1, 1E-2)
    initial_conditions_barely_trapped() = initial_conditions([2.5 / equ.R₀, 0., 0.], 3.375E-1, 1E-2)
    initial_conditions_deeply_passing() = initial_conditions([2.5 / equ.R₀, 0., 0.],  5E-1,    1E-2)
    initial_conditions_deeply_trapped() = initial_conditions([2.5 / equ.R₀, 0., 0.],  1E-1,    1E-2)
    initial_conditions_trapped()        = initial_conditions([7.0 / equ.R₀, 0., 0.], -2E-3,    1.88E-7)

end
