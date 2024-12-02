module TokamakSmallCartesian

    import ElectromagneticFields.AxisymmetricTokamakCartesian

    export charged_particle_3d_pode, charged_particle_3d_iode,
           hamiltonian, toroidal_momentum

    AxisymmetricTokamakCartesian.@code() # inject magnetic field code
       
    include("charged_particle_3d_canonical.jl")

    const qᵢ = [1.05,   0.0,    0.0]
    const vᵢ = [2.1E-3, 4.3E-4, 0.0]
    const pᵢ = charged_particle_3d_pᵢ(tᵢ, qᵢ, vᵢ)

end
