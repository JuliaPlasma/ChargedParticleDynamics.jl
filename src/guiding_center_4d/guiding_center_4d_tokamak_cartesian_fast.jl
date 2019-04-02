"""
Fast barely passing particle in analytic, axisymmetric tokamak equilibrium in cartesian coordinates.
"""
module TokamakCartesianFastBarelyPassing

    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCartesian

    # Δt=2.5, nt=50

    const μ  = 1E-2
    const qᵢ = [2.5, 0., 0., 3.425E-1]

    equ = AxisymmetricTokamakCartesian(2., 5., 2.)
    load_equilibrium(equ; target_module=TokamakCartesianFastBarelyPassing)

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")

end

"""
Fast barely trapped particle in analytic, axisymmetric tokamak equilibrium.
"""
module TokamakCartesianFastBarelyTrapped

    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCartesian

    # Δt=3, nt=100

    const μ  = 1E-2
    const qᵢ = [2.5, 0., 0., 3.375E-1]

    equ = AxisymmetricTokamakCartesian(2., 5., 2.)
    load_equilibrium(equ; target_module=TokamakCartesianFastBarelyTrapped)

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")

end


"""
Fast deeply passing particle in analytic, axisymmetric tokamak equilibrium.
"""
module TokamakCartesianFastDeeplyPassing

    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCartesian

    # Δt=2.5, nt=25

    const μ  = 1E-2
    const qᵢ = [2.5, 0., 0., 5E-1]

    equ = AxisymmetricTokamakCartesian(2., 5., 2.)
    load_equilibrium(equ; target_module=TokamakCartesianFastDeeplyPassing)

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")

end


"""
Fast deeply trapped particle in analytic, axisymmetric tokamak equilibrium.
"""
module TokamakCartesianFastDeeplyTrapped

    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCartesian

    # Δt=5, nt=50

    const μ  = 1E-2
    const qᵢ = [2.5, 0., 0., 1E-1]

    equ = AxisymmetricTokamakCartesian(2., 5., 2.)
    load_equilibrium(equ; target_module=TokamakCartesianFastDeeplyTrapped)

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")

end


"""
Fast particles on a phasespace loop in analytic, axisymmetric tokamak equilibrium.
"""
module TokamakCartesianFastLoop

    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCartesian

    const μ  = 1E-3

    equ = AxisymmetricTokamakCartesian(2., 5., 2.)
    load_equilibrium(equ; target_module=TokamakCartesianFastLoop)

    function f_loop(s)
        X0 = 1.75
        Y0 = 0.0
        Z0 = 0.0
        u0 = 0.5
        rx = 0.1
        ry = 0.1

        Xs = X0 + rx*cos(2π*s)
        Zs = Z0 + ry*sin(2π*s)

        qs = [Xs, Y0, Zs, u0]

        return qs
    end

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")
    include("guiding_center_4d_loop.jl")

end


"""
Fast particles on a phasespace surface in analytic, axisymmetric tokamak equilibrium.
"""
module TokamakCartesianFastSurface

    using ElectromagneticFields: load_equilibrium, periodicity, AxisymmetricTokamakCartesian

    const μ  = 1E-3

    equ = AxisymmetricTokamakCartesian(2., 5., 2.)
    load_equilibrium(equ; target_module=TokamakCartesianFastSurface)

    function f_surface(s,t)
        X0 = 1.75
        Y0 = 0.0
        Z0 = 0.0
        u0 = 0.5
        r0 = 0.1

        Xt = X0 + r0*(s-0.5)
        Zt = Z0 + r0*(t-0.5)

        qt = [Xt, Y0, Zt, u0]

        return qt
    end

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")
    include("guiding_center_4d_surface.jl")

end
