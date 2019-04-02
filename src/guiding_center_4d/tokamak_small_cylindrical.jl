"""
Analytic axisymmetric medium-size tokamak equilibrium in cartesian coordinates.
"""
module TokamakSmallCylindrical

    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCylindrical

    export initial_conditions_barely_passing, initial_conditions_barely_trapped,
           initial_conditions_deeply_passing, initial_conditions_deeply_trapped

    export toroidal_momentum

    μ = 2.5E-6

    equ = AxisymmetricTokamakCylindrical(1., 1., 2.)
    load_equilibrium(equ; target_module=TokamakSmallCylindrical)


    function initial_conditions_barely_passing()
        μ  = 2.448E-6
        return [1.05, 0., 0., 8.117E-4]
    end

    function initial_conditions_barely_trapped()
        μ  = 2.250E-6
        return [1.05, 0., 0., 7.610E-4]
    end

    function initial_conditions_deeply_passing()
        μ  = 2.448E-6
        return [1.05, 0., 0., 1.623E-3]
    end

    function initial_conditions_deeply_trapped()
        μ  = 2.250E-6
        return [1.05, 0., 0., 4.306E-4]
    end

    function f_loop(t)
        R0 = 1.0
        Z0 = 0.0
        φ0 = 0.0
        u0 = 4E-4
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
        u0 = 4E-4
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

end
