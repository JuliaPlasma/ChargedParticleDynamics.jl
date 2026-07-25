module TokamakSmallToroidal

    import ElectromagneticFields.AxisymmetricTokamakToroidal

    export charged_particle_3d_pode, charged_particle_3d_iode,
           hamiltonian, toroidal_momentum

    AxisymmetricTokamakToroidal.@code() # inject magnetic field code

    include("charged_particle_3d_canonical.jl")

    export default_parameters

    """
    The charged particle models are parameter-free — the electromagnetic field is injected as
    code rather than passed as parameters. The method exists so that every problem in this
    package can be constructed the same way.
    """
    default_parameters(::Type{T}=Float64) where {T} = NamedTuple()

    const parameters = default_parameters()

    const xᵢ = [1.05, 0.0, 0.0]
    const qᵢ = from_cartesian(0, xᵢ)
    const vᵢ = DF̄(0, qᵢ) * [2.1E-3, 4.3E-4, 0.0]
    const pᵢ = charged_particle_3d_pᵢ(tᵢ, qᵢ, vᵢ)

end
