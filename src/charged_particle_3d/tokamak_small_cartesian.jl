module TokamakSmallCartesian

    using ElectromagneticFields.AxisymmetricTokamakCartesian

    export charged_particle_3d_pode, charged_particle_3d_iode,
           hamiltonian, toroidal_momentum

    const R₀ = 1.
    const B₀ = 1.
    const q  = 2.

    const equ = AxisymmetricTokamakCartesian.init(R₀, B₀, q)
       
    include("charged_particle_3d_canonical.jl")

    const qᵢ = [1.05,   0.0,    0.0]
    const vᵢ = [2.1E-3, 4.3E-4, 0.0]
    const pᵢ = charged_particle_3d_pᵢ(qᵢ, vᵢ)

end
