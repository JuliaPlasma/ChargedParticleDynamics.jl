module TokamakSmallCircular

    using ElectromagneticFields.AxisymmetricTokamakCircular
    using LinearAlgebra

    export charged_particle_3d_pode, charged_particle_3d_iode,
           hamiltonian, toroidal_momentum

    const R₀ = 1.
    const B₀ = 1.
    const q  = 2.

    const equ = AxisymmetricTokamakCircular.init(R₀, B₀, q)

    include("charged_particle_3d_canonical.jl")

    const xᵢ = [1.05, 0.0, 0.0]
    const qᵢ = from_cartesian(0, xᵢ)
    const vᵢ = DFinv(0, xᵢ) * [2.1E-3, 4.3E-4, 0.0]
    const pᵢ = charged_particle_3d_pᵢ(qᵢ, vᵢ)

end
