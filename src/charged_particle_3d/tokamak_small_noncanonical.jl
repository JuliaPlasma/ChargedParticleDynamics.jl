module TokamakSmallNoncanonical

    using ElectromagneticFields.AxisymmetricTokamakToroidal

    export charged_particle_3d_ode, charged_particle_3d_sode, charged_particle_3d_iode,
           hamiltonian, toroidal_momentum, ϑ

    const R₀ = 1.
    const B₀ = 1.
    const q  = 2.

    const equ = AxisymmetricTokamakToroidal.init(R₀, B₀, q)
       
    const xᵢ = from_cartesian(0, [1.05,   0.0,    0.0])
    const vᵢ = from_cartesian(0, [2.1E-3, 4.3E-4, 0.0])
    const qᵢ = vcat(xᵢ, vᵢ)

    include("charged_particle_3d_noncanonical.jl")
    
end
