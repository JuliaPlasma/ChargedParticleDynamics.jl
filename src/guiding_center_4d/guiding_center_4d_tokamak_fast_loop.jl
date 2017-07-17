module TokamakFastLoop

    using GeometricIntegrators
    using MagneticEquilibria

    export guiding_center_4d_ode_poincare_invariant_1st,
           guiding_center_4d_iode_poincare_invariant_1st,
           hamiltonian, toroidal_momentum, u, α

    const R₀ = 2
    const B₀ = 5
    const q  = 2
    const μ  = 1E-2

    load_equilibrium(MagneticEquilibria.AxisymmetricTokamak(R₀, B₀, q); target_module=TokamakFastLoop)

    include("guiding_center_4d_RZphi.jl")


    function f_loop(s)
        R0 = 1.75
        Z0 = 0.0
        φ0 = 0.0
        u0 = 0.5
        rx = 0.75
        ry = 0.75

        Rs = R0 + rx*cos(2π*s)
        Zs = Z0 + ry*sin(2π*s)

        qs = [xs, ys, φ0, u0]

        return qs
    end

    # function f_loop(i, n)
    #     f_loop(i/n)
    # end
    #
    # function get_initial_conditions(n)
    #     q₀ = zeros(4,n)
    #
    #     for i in 1:size(q₀,2)
    #         q₀[:,i] = f_loop(i,n)
    #     end
    #
    #     return q₀
    # end


    function guiding_center_4d_ode_poincare_invariant_1st(Δt, nloop, ntime, nsave, DT=Float64)
        PoincareInvariant1st(guiding_center_4d_ode, f_loop, α, Δt, 4, nloop, ntime, nsave, DT)
    end

    function guiding_center_4d_iode_poincare_invariant_1st(Δt, nloop, ntime, nsave, DT=Float64)
        PoincareInvariant1st(guiding_center_4d_iode, f_loop, α, Δt, 4, nloop, ntime, nsave, DT)
    end

end
