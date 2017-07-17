module SymmetricLoop

    using GeometricIntegrators
    using MagneticEquilibria

    export guiding_center_4d_ode_poincare_invariant_1st,
           guiding_center_4d_iode_poincare_invariant_1st,
           hamiltonian, toroidal_momentum, u, α

    const B₀ = 1.
    const μ  = 1E-2

    load_equilibrium(SymmetricQuadratic(B₀); target_module=SymmetricLoop)

    include("guiding_center_4d_xyz.jl")


    function f_loop(s)
        z0 = 0.0
        u0 = 0.5
        u1 = 0.01
        rx = 0.5
        ry = 0.3
        rz = 0.2

        xs = rx*cos(2π*s)
        ys = ry*sin(2π*s)
        # zs = rz*sin(2π*s)
        # us = u0 + u1 * cos(2π*s)

        qs = [xs, ys, z0, u0]
        # qs = [xs, ys, zs, u0]
        # qs = [xs, ys, zs, us]

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
