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
        rx = 0.5
        ry = 0.3
        z0 = 0.0
        z1 = 0.1
        u0 = 5E-1
        u1 = 5E-2

        xs = rx*cos(2π*s)
        ys = ry*sin(2π*s)
        zs = z0 + z1 * sin(2π*s)
        us = u0 + u1 * cos(2π*s)

        qs = [xs, ys, zs, us]

        return qs
    end


    function guiding_center_4d_ode_poincare_invariant_1st(Δt, nloop, ntime, nsave, DT=Float64)
        guiding_center_4d_ode_init(q₀) = guiding_center_4d_ode(q₀; periodic=false)
        PoincareInvariant1st(guiding_center_4d_ode_init, f_loop, α, Δt, 4, nloop, ntime, nsave, DT)
    end

    function guiding_center_4d_iode_poincare_invariant_1st(Δt, nloop, ntime, nsave, DT=Float64)
        guiding_center_4d_iode_init(q₀) = guiding_center_4d_iode(q₀; periodic=false)
        PoincareInvariant1st(guiding_center_4d_iode_init, f_loop, α, Δt, 4, nloop, ntime, nsave, DT)
    end

end
