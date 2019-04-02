module GuidingCenter4dThetaPinch

    using ElectromagneticFields: load_equilibrium, periodicity, ThetaPinch

    const μ  = 2.5E-6

    equ = ThetaPinch(1.)
    load_equilibrium(equ; target_module=GuidingCenter4dThetaPinch)

    function f_loop(t)
        φ0 = 0.0
        u0 = 4E-4
        r0 = 0.5
        r1 = 0.3

        Rt = r0*cos(2π*t)
        Zt = r1*sin(2π*t)

        qt = [Rt, Zt, φ0, u0]

        return qt
    end

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")
    include("guiding_center_4d_loop.jl")

end
