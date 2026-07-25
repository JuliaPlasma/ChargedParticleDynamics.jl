"""
Charged Particle in a singular magnetic field of the form
``B(x,y,z) = (x^2 + y^2)^{-3/2} e_z``.
"""
module SingularField

    import ElectromagneticFields.Singular

    export charged_particle_3d_ode, charged_particle_3d_iode
    export hamiltonian, angular_momentum
    export compute_energy, compute_energy_error


    Singular.@code() # inject magnetic field code

    const qᵢ = [1., 0., 0., 0., -1., 0.]

    angular_momentum(t,q) = q[1] * ϑ₂(t,q) - q[2] * ϑ₁(t,q)

    include("charged_particle_3d_noncanonical.jl")

    export default_parameters

    """
    The charged particle models are parameter-free — the electromagnetic field is injected as
    code rather than passed as parameters. The method exists so that every problem in this
    package can be constructed the same way.
    """
    default_parameters(::Type{T}=Float64) where {T} = NamedTuple()

    const parameters = default_parameters()

end
