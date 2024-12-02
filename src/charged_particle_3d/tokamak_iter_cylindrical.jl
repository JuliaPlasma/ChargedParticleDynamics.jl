module TokamakIterCylindrical

    import ElectromagneticFields.AxisymmetricTokamakCylindrical

    export charged_particle_3d_pode, charged_particle_3d_iode,
           hamiltonian, toroidal_momentum

    AxisymmetricTokamakCylindrical.@code_iter() # inject magnetic field code

    include("charged_particle_3d_canonical.jl")

    const qᵢ = [7.0, 0.0, 0.0]
    const vᵢ = [3.43E-3, 6.75, -3.41E-1]
    const pᵢ = charged_particle_3d_pᵢ(tᵢ, qᵢ, vᵢ)

end
