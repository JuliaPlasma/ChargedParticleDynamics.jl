"""
Fast barely passing particle in analytic, axisymmetric tokamak equilibrium.
"""
module TokamakFastBarelyPassing

    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCylindrical
    export toroidal_momentum

    # Δt=2.5, nt=50

    const μ  = 1E-2
    const qᵢ = [2.5, 0., 0., 3.425E-1]

    equ = AxisymmetricTokamakCylindrical(2., 5., 2.)
    load_equilibrium(equ; target_module=TokamakFastBarelyPassing)

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")

    function toroidal_momentum(t,q)
        R(t,q) * ϑ3(t,q)
    end

end


"""
Fast barely trapped particle in analytic, axisymmetric tokamak equilibrium.
"""
module TokamakFastBarelyTrapped

    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCylindrical
    export toroidal_momentum

    # Δt=3, nt=100

    const μ  = 1E-2
    const qᵢ = [2.5, 0., 0., 3.375E-1]

    equ = AxisymmetricTokamakCylindrical(2., 5., 2.)
    load_equilibrium(equ; target_module=TokamakFastBarelyTrapped)

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")

    function toroidal_momentum(t,q)
        R(t,q) * ϑ3(t,q)
    end

end


"""
Fast deeply passing particle in analytic, axisymmetric tokamak equilibrium.
"""
module TokamakFastDeeplyPassing

    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCylindrical
    export toroidal_momentum

    # Δt=2.5, nt=25

    const μ  = 1E-2
    const qᵢ = [2.5, 0., 0., 5E-1]

    equ = AxisymmetricTokamakCylindrical(2., 5., 2.)
    load_equilibrium(equ; target_module=TokamakFastDeeplyPassing)

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")

    function toroidal_momentum(t,q)
        R(t,q) * ϑ3(t,q)
    end

end


"""
Fast deeply trapped particle in analytic, axisymmetric tokamak equilibrium.
"""
module TokamakFastDeeplyTrapped

    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCylindrical
    export toroidal_momentum

    # Δt=5, nt=50

    const μ  = 1E-2
    const qᵢ = [2.5, 0., 0., 1E-1]

    equ = AxisymmetricTokamakCylindrical(2., 5., 2.)
    load_equilibrium(equ; target_module=TokamakFastDeeplyTrapped)

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")

    function toroidal_momentum(t,q)
        R(t,q) * ϑ3(t,q)
    end

end


"""
Fast particles on a phasespace loop in analytic, axisymmetric tokamak equilibrium.
"""
module TokamakFastLoop

    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCylindrical

    const μ  = 1E-3

    equ = AxisymmetricTokamakCylindrical(2., 5., 2.)
    load_equilibrium(equ; target_module=TokamakFastLoop)

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

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")
    include("guiding_center_4d_loop.jl")

end


"""
Fast particles on a phasespace surface in analytic, axisymmetric tokamak equilibrium.
"""
module TokamakFastSurface

    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCylindrical

    const μ  = 1E-3

    equ = AxisymmetricTokamakCylindrical(2., 5., 2.)
    load_equilibrium(equ; target_module=TokamakFastSurface)

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
    include("guiding_center_4d_surface.jl")

end
