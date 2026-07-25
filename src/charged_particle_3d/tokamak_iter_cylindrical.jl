module TokamakIterCylindrical

    import ElectromagneticFields.AxisymmetricTokamakCylindrical

    export charged_particle_3d_pode, charged_particle_3d_iode,
           hamiltonian, toroidal_momentum

    AxisymmetricTokamakCylindrical.@code_iter() # inject magnetic field code

    include("charged_particle_3d_canonical.jl")

    export default_parameters

    """
    The charged particle models are parameter-free — the electromagnetic field is injected as
    code rather than passed as parameters. The method exists so that every problem in this
    package can be constructed the same way.
    """
    default_parameters(::Type{T}=Float64) where {T} = NamedTuple()

    const parameters = default_parameters()

    const qᵢ = [7.0, 0.0, 0.0]
    const vᵢ = [3.43E-3, 6.75, -3.41E-1]
    const pᵢ = charged_particle_3d_pᵢ(tᵢ, qᵢ, vᵢ)

end
