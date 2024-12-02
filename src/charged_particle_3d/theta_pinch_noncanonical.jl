module ThetaPinchNoncanonical

    import ElectromagneticFields.ThetaPinch

    export charged_particle_3d_ode, charged_particle_3d_sode, charged_particle_3d_iode,
           hamiltonian, toroidal_momentum, ϑ

    const qᵢ = [2.5, 0.0, 0.0, 0.0, 0.2, 0.1]

    ThetaPinch.@code() # inject magnetic field code

    toroidal_momentum(t,q) = ϑ₃(t,q)

    include("charged_particle_3d_noncanonical.jl")
    
end
