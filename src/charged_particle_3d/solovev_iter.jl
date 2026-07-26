@doc raw"""
Analytic ITER-like Solov'ev equilibrium, in cylindrical coordinates ``(R, Z, \varphi)``.

The Solov'ev solution of the Grad-Shafranov equation fitted to ITER's shape, without an
X-point. See `SolovevIterXpoint` for the diverted variant.

The hard-coded `μ` and `u` of the `initial_conditions_*` below are inherited from earlier work on
this package and their provenance is not recorded; they are known-good starting points for the
charged particle model rather than values derived here.
"""
module SolovevIter

    import ElectromagneticFields.Solovev

    export podeproblem, iodeproblem,
           hamiltonian, toroidal_momentum

    Solovev.@code_iter() # inject magnetic field code

    include("charged_particle_3d_canonical.jl")

    export default_parameters

    """
    The charged particle models are parameter-free — the electromagnetic field is injected as
    code rather than passed as parameters. The method exists so that every problem in this
    package can be constructed the same way.
    """
    default_parameters(::Type{T}=Float64) where {T} = NamedTuple()


    const xᵢ = [7.0 - 1.4, 0.0, 0.0]
    const qᵢ = from_cartesian(0, xᵢ)
    const vᵢ = [3.43E-3, 6.75, -3.41E-1]
    const pᵢ = charged_particle_3d_pᵢ(tᵢ, qᵢ, vᵢ)

    function initial_conditions(x₀, u₀, μ)
        vpar = u₀ * b⃗(0, x₀)
        vper = sqrt(2 * μ * B(0, x₀))
        v¹ = vper * sqrt( b³(0, x₀)^2 / (b¹(0, x₀)^2 + b³(0, x₀)^2) )
        v³ = - v¹ * b¹(0, x₀) / b³(0, x₀)
        v₀ = vpar .+ [v¹, 0, v³]

        (q = x₀, p = charged_particle_3d_pᵢ(tᵢ, x₀, v₀), params = default_parameters())
    end

    initial_conditions_barely_passing() = initial_conditions(from_cartesian(0, [2.5, 0., 0.]), 3.425E-1, 1E-2)
    initial_conditions_barely_trapped() = initial_conditions(from_cartesian(0, [2.5, 0., 0.]), 3.375E-1, 1E-2)
    initial_conditions_deeply_passing() = initial_conditions(from_cartesian(0, [2.5, 0., 0.]),  5E-1,    1E-2)
    initial_conditions_deeply_trapped() = initial_conditions(from_cartesian(0, [2.5, 0., 0.]),  1E-1,    1E-2)
    initial_conditions_trapped()        = initial_conditions(from_cartesian(0, [2.5, 0., 0.]), -2E-3,    1.88E-7)

end
