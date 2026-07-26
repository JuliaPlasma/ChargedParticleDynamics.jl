@doc raw"""
Small axisymmetric tokamak equilibrium, in toroidal coordinates ``(r, \theta, \varphi)``, in the
noncanonical formulation on the phasespace ``z = (x, v)``.

Major radius ``R_{0} = 1``, magnetic field on axis ``B_{0} = 1`` and safety factor ``q_{0} = 2`` —
the same equilibrium as `TokamakSmallToroidal`.

The only curvilinear module of the four noncanonical ones, which makes it the one where the metric
terms of the vector field matter. `sodeproblem` is therefore unavailable: the frozen-position kick
is quadratic in ``v`` in a curvilinear chart, so the Boris push is not its exact flow. Note also
that the vector field carries a ``1/g_{ii}`` and ``g_{22} = r^{2}``, so it is singular on the
magnetic axis.

The hard-coded initial velocity is inherited from earlier work on this package and its provenance is
not recorded; it is a known-good starting point rather than a value derived here.
"""
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

    
end
