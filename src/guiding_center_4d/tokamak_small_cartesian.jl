"""
Analytic axisymmetric small tokamak equilibrium in cartesian coordinates.
"""
module TokamakSmallCartesian

    import ElectromagneticFields.AxisymmetricTokamakCartesian

    export initial_conditions_barely_passing, initial_conditions_barely_trapped,
           initial_conditions_deeply_passing, initial_conditions_deeply_trapped,
           initial_conditions_pauli

    export hamiltonian, toroidal_momentum

    AxisymmetricTokamakCartesian.@code() # inject magnetic field code

    const Δt = 500.0
    const tspan = (0.0, 5E4)

    const qᵢ = [1.05, 0., 0., 0.00045135897235326736]

    initial_conditions_barely_passing() = ([1.05, 0., 0., 8.117E-4], (μ = 2.448E-6,))
    initial_conditions_barely_trapped() = ([1.05, 0., 0., 7.610E-4], (μ = 2.250E-6,))
    initial_conditions_deeply_passing() = ([1.05, 0., 0., 1.623E-3], (μ = 2.448E-6,))
    initial_conditions_deeply_trapped() = ([1.05, 0., 0., 4.306E-4], (μ = 2.250E-6,))
    initial_conditions_pauli() = ([1.05, 0., 0., 4.3E-4], (μ = 2.310E-6,))

    u_loop() = 4.0E-4
    μ_loop() = 2.5E-6
    u_surface() = 4.0E-4
    μ_surface() = 2.5E-6

    function f_loop(t)
        R0 = 1.0
        Z0 = 0.0
        φ0 = 0.0
        u0 = u_loop()
        r0 = 0.05

        Rt = R0 + r0*cos(2π*t)
        Zt = Z0 + r0*sin(2π*t)

        qt = [Rt, Zt, φ0, u0]

        return qt
    end

    function f_surface(s,t)
        R0 = 1.0
        Z0 = 0.0
        φ0 = 0.0
        u0 = u_surface()
        r0 = 0.1

        Rt = R0 + 2r0*(s-0.5)
        Zt = Z0 + 2r0*(t-0.5)

        qt = [Rt, Zt, φ0, u0]

        return qt
    end

    export default_parameters

    "The magnetic moment μ this equilibrium is set up for."
    default_parameters(::Type{T}=Float64) where {T} = (μ = T(2.314593645825811e-6),)

    const parameters = default_parameters()

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")
    include("guiding_center_4d_loop.jl")
    include("guiding_center_4d_surface.jl")

    # In cartesian coordinates the third coordinate is z, not an angle, so the toroidal momentum is
    # the generator of rotation about the z-axis, x ϑ₂ - y ϑ₁, rather than ϑ₃. Neither ϑ₃ nor R ϑ₃
    # is conserved here; this is, to the order of the integrator.
    function toroidal_momentum(t,q)
        q[1] * ϑ₂(t,q) - q[2] * ϑ₁(t,q)
    end

    include("guiding_center_4d_diagnostics.jl")

end
