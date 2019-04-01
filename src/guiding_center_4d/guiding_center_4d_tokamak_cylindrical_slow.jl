"""
Slow barely passing particle in analytic, axisymmetric tokamak equilibrium.
"""
module TokamakSlowBarelyPassing

    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCylindrical
    export toroidal_momentum

    const μ  = 2.448E-6
    const qᵢ = [1.05, 0., 0., 8.117E-4]

    equ = AxisymmetricTokamakCylindrical(1., 1., 2.)
    load_equilibrium(equ; target_module=TokamakSlowBarelyPassing)

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")

    function toroidal_momentum(t,q)
        R(t,q) * ϑ3(t,q)
    end

end


"""
Slow barely trapped particle in analytic, axisymmetric tokamak equilibrium.
"""
module TokamakSlowBarelyTrapped

    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCylindrical
    export toroidal_momentum

    const μ  = 2.250E-6
    const qᵢ = [1.05, 0., 0., 7.610E-4]

    equ = AxisymmetricTokamakCylindrical(1., 1., 2.)
    load_equilibrium(equ; target_module=TokamakSlowBarelyTrapped)

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")

    function toroidal_momentum(t,q)
        R(t,q) * ϑ3(t,q)
    end

end


"""
Slow deeply passing particle in analytic, axisymmetric tokamak equilibrium.
"""
module TokamakSlowDeeplyPassing

    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCylindrical
    export toroidal_momentum

    const μ  = 2.448E-6
    const qᵢ = [1.05, 0., 0., 1.623E-3]

    equ = AxisymmetricTokamakCylindrical(1., 1., 2.)
    load_equilibrium(equ; target_module=TokamakSlowDeeplyPassing)

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")

    function toroidal_momentum(t,q)
        R(t,q) * ϑ3(t,q)
    end

end


"""
Slow deeply trapped particle in analytic, axisymmetric tokamak equilibrium.
"""
module TokamakSlowDeeplyTrapped

    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCylindrical
    export toroidal_momentum

    const μ  = 2.250E-6
    const qᵢ = [1.05, 0., 0., 4.306E-4]

    equ = AxisymmetricTokamakCylindrical(1., 1., 2.)
    load_equilibrium(equ; target_module=TokamakSlowDeeplyTrapped)

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")

    function toroidal_momentum(t,q)
        R(t,q) * ϑ3(t,q)
    end

end


"""
Slow particles on a phasespace loop in analytic, axisymmetric tokamak equilibrium.
"""
module TokamakSlowLoop

    export guiding_center_4d_loop_ode, guiding_center_4d_loop_iode
    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCylindrical

    const μ  = 2.5E-6


    include("guiding_center_4d_common.jl")

    equ = AxisymmetricTokamakCylindrical(1., 1., 2.)
    load_equilibrium(equ; target_module=TokamakSlowLoop)

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

    function f_loop(i, n)
        f_loop(i/n)
    end

    function get_initial_conditions(n)
        q₀ = zeros(4, n)

        for i in 1:n
            q₀[:,i] = f_loop(i, n)
        end

        return q₀
    end


    function guiding_center_4d_loop_ode(n)
        q₀ = get_initial_conditions(n)
        guiding_center_4d_ode(q₀; periodic=false)
    end

    function guiding_center_4d_loop_iode(n)
        q₀ = get_initial_conditions(n)
        guiding_center_4d_iode(q₀; periodic=false)
    end

end


"""
Slow particles on a phasespace surface in analytic, axisymmetric tokamak equilibrium.
"""
module TokamakSlowSurface

    export guiding_center_4d_surface_ode, guiding_center_4d_surface_iode
    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCylindrical

    const μ  = 2.5E-6


    include("guiding_center_4d_common.jl")

    equ = AxisymmetricTokamakCylindrical(1., 1., 2.)
    load_equilibrium(equ; target_module=TokamakSlowSurface)

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

    function f_surface(i, j, nx, ny)
        f_surface(i/nx, j/ny)
    end

    function get_initial_conditions(nx, ny)
        q₀ = zeros(4, nx*ny)

        for j in 1:ny
            for i in 1:nx
                q₀[:,nx*(j-1)+i] = f_surface(i, j, nx, ny)
            end
        end

        return q₀
    end


    function guiding_center_4d_surface_ode(nx, ny)
        q₀ = get_initial_conditions(nx, ny)
        guiding_center_4d_ode(q₀; periodic=false)
    end

    function guiding_center_4d_surface_iode(nx, ny)
        q₀ = get_initial_conditions(nx, ny)
        guiding_center_4d_iode(q₀; periodic=false)
    end

end
