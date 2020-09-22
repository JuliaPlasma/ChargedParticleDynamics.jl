"""
Analytic axisymmetric small tokamak equilibrium in toroidal coordinates.
"""
module TokamakSmallToroidal

    using ElectromagneticFields.AxisymmetricTokamakToroidal

    export initial_conditions_barely_passing, initial_conditions_barely_trapped,
           initial_conditions_deeply_passing, initial_conditions_deeply_trapped,
           initial_conditions_pauli

    export hamiltonian, toroidal_momentum

    const R₀ = 1.
    const B₀ = 1.
    const q  = 2.

    const equ = AxisymmetricTokamakToroidal.init(R₀, B₀, q)

    const qᵢ = [1.05, 0., 0., 0.00045135897235326736]
    const parameters = (μ = 2.314593645825811e-6,)

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

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")
    include("guiding_center_4d_loop.jl")
    include("guiding_center_4d_surface.jl")

    function toroidal_momentum(t,q)
        R(t,q) * ϑ3(t,q)
    end

    include("guiding_center_4d_diagnostics.jl")
    
end
