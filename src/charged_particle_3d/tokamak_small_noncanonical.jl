module TokamakSmallNoncanonical

    import ElectromagneticFields.AxisymmetricTokamakToroidal

    export charged_particle_3d_ode, charged_particle_3d_sode, charged_particle_3d_iode,
           hamiltonian, toroidal_momentum, ϑ

    AxisymmetricTokamakToroidal.@code() # inject magnetic field code
       
    const xᵢ = [1.05, 0.0, 0.0]
    const qᵢ = vcat(from_cartesian(0, xᵢ), DF̄(0, xᵢ) * [2.1E-3, 4.3E-4, 0.0])

    toroidal_momentum(t,q) = ϑ₃(t,q)

    include("charged_particle_3d_noncanonical.jl")
    
end
