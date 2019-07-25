module GuidingCenter4dThetaPinch

    using ElectromagneticFields: load_equilibrium, periodicity, ThetaPinch


    equ = ThetaPinch(1.)
    load_equilibrium(equ; target_module=GuidingCenter4dThetaPinch)

    μ_loop() = 2.5E-6

    function f_loop(t)
        μ  = 2.5E-6
        Y0 = 0.0
        u0 = 4E-4
        r0 = 0.5
        r1 = 0.3

        Xt = r0*cos(2π*t)
        Zt = r1*sin(2π*t)

        qt = [Xt, Y0, Zt, u0]

        return qt
    end

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")
    include("guiding_center_4d_loop.jl")

end
