module TokamakSmallNoncanonical

    import ElectromagneticFields.AxisymmetricTokamakToroidal

    # `sodeproblem` is deliberately not exported here: this equilibrium is toroidal,
    # and the Boris splitting is only a valid splitting of the model where the metric is trivial.
    # The constructor throws if called anyway; see its docstring.
    export odeproblem, iodeproblem,
           hamiltonian, toroidal_momentum, ϑ

    AxisymmetricTokamakToroidal.@code() # inject magnetic field code
       
    const xᵢ = [1.05, 0.0, 0.0]
    const qᵢ = vcat(from_cartesian(0, xᵢ), DF̄(0, xᵢ) * [2.1E-3, 4.3E-4, 0.0])

    toroidal_momentum(t,q) = ϑ₃(t,q)

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
