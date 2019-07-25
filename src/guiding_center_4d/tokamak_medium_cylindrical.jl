"""
Analytic axisymmetric medium-size tokamak equilibrium in cartesian coordinates.
"""
module TokamakMediumCylindrical

    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCylindrical

    export initial_conditions_barely_passing, initial_conditions_barely_trapped,
           initial_conditions_deeply_passing, initial_conditions_deeply_trapped

    export toroidal_momentum


    equ = AxisymmetricTokamakCylindrical(2., 5., 2.)
    load_equilibrium(equ; target_module=TokamakMediumCylindrical)

    initial_conditions_barely_passing() = ([2.5, 0., 0., 3.425E-1], 1E-2) # Δt=2.5, nt=50
    initial_conditions_barely_trapped() = ([2.5, 0., 0., 3.375E-1], 1E-2) # Δt=3.0, nt=100
    initial_conditions_deeply_passing() = ([2.5, 0., 0., 5E-1], 1E-2)     # Δt=2.5, nt=25
    initial_conditions_deeply_trapped() = ([2.5, 0., 0., 1E-1], 1E-2)     # Δt=5.0, nt=50

    μ_loop() = 1E-3
    μ_surface() = 1E-3

    function f_loop(s)
        R0 = 1.75
        Z0 = 0.0
        φ0 = 0.0
        u0 = 0.5
        rx = 0.1
        ry = 0.1

        Rs = R0 + rx*cos(2π*s)
        Zs = Z0 + ry*sin(2π*s)

        qs = [Rs, Zs, φ0, u0]

        return qs
    end

    function f_surface(s,t)
        R0 = 1.75
        Z0 = 0.0
        φ0 = 0.0
        u0 = 0.5
        r0 = 0.1

        Rt = R0 + r0*(s-0.5)
        Zt = Z0 + r0*(t-0.5)

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

end
