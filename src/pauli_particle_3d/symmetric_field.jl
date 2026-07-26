"""
Charged Particle in an axisymmetric magnetic field of the form
``B(x,y,z) = (1 + x^2 + y^2) e_z``.
"""
module SymmetricField

    import ElectromagneticFields.SymmetricQuadratic

    export podeproblem, hamiltonian, toroidal_momentum

    SymmetricQuadratic.@code() # inject magnetic field code

    const qᵢ = [1., 0., 0.]
    const vᵢ = [0., 1., 1.]

    const DEFAULT_TIMESTEP = 1.0
    const DEFAULT_TIMESPAN = (0.0, 1000.0)

    include("pauli_particle_3d.jl")

    export default_parameters

    """
    The magnetic moment μ of the default initial condition `(qᵢ, vᵢ)`, obtained by splitting
    `vᵢ` into its parallel and perpendicular parts at `qᵢ`.
    """
    default_parameters(::Type{T}=Float64) where {T} = (μ = T(initial_conditions(qᵢ, vᵢ).params.μ),)

    const parameters = default_parameters()

end
