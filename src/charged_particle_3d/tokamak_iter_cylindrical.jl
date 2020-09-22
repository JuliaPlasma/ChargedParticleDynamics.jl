module TokamakSmallCylindrical

    using ElectromagneticFields.AxisymmetricTokamakCylindrical

    export charged_particle_3d_pode, charged_particle_3d_iode,
           hamiltonian, toroidal_momentum

    const R₀ = 6.2
    const B₀ = 5.3
    const q  = 2.

    const equ = AxisymmetricTokamakCylindrical.init(R₀, B₀, q)

    include("charged_particle_3d_canonical.jl")

    const qᵢ = [7.0, 0.0, 0.0]
    const vᵢ = [3.43E-3, 6.75, -3.41E-1]
    const pᵢ = charged_particle_3d_pᵢ(qᵢ, vᵢ)

end
