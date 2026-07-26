@doc raw"""
Charged particle in a uniform magnetic field of the form ``B(x,y,z) = B_{0} e_{z}``, in the
noncanonical formulation on the phasespace ``z = (x, v)``.

The same equilibrium as `ThetaPinchCanonical`, in cartesian coordinates, so the metric is trivial
and the Boris splitting of `sodeproblem` is a valid splitting of the model here.
"""
module ThetaPinchNoncanonical

    import ElectromagneticFields.ThetaPinch

    export odeproblem, sodeproblem, iodeproblem,
           hamiltonian, toroidal_momentum, ϑ

    const qᵢ = [2.5, 0.0, 0.0, 0.0, 0.2, 0.1]

    ThetaPinch.@code() # inject magnetic field code

    toroidal_momentum(t,q) = ϑ₃(t,q)

    include("charged_particle_3d_noncanonical.jl")

    export default_parameters

    """
    The charged particle models are parameter-free — the electromagnetic field is injected as
    code rather than passed as parameters. The method exists so that every problem in this
    package can be constructed the same way.
    """
    default_parameters(::Type{T}=Float64) where {T} = NamedTuple()

    
end
