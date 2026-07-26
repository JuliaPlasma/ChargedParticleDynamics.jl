"""
Charged Particle in an uniform magnetic field of the form
``B(x,y,z) = B_0 e_z``.
"""
module ThetaPinchCanonical

    import ElectromagneticFields.ThetaPinch

    export iodeproblem, hamiltonian, angular_momentum

    ThetaPinch.@code() # inject magnetic field code

    include("charged_particle_3d_canonical.jl")

    export default_parameters

    """
    The charged particle models are parameter-free — the electromagnetic field is injected as
    code rather than passed as parameters. The method exists so that every problem in this
    package can be constructed the same way.
    """
    default_parameters(::Type{T}=Float64) where {T} = NamedTuple()

    const parameters = default_parameters()

    const qᵢ = [2.5, 0.0, 0.0]
    const vᵢ = [0.0, 0.2, 0.1]
    const pᵢ = charged_particle_3d_pᵢ(tᵢ, qᵢ, vᵢ)

end
