module TokamakFastLoop

    using GeometricIntegrators
    using MagneticEquilibria

    export guiding_center_4d_ode_poincare_invariant_1st,
           guiding_center_4d_iode_poincare_invariant_1st,
           hamiltonian, toroidal_momentum, u, α

    const R₀ = 2
    const B₀ = 5
    const q  = 2
    const μ  = 1E-3

    load_equilibrium(MagneticEquilibria.AxisymmetricTokamak(R₀, B₀, q); target_module=TokamakFastLoop)

    include("guiding_center_4d_coords_RZphi.jl")


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


    function guiding_center_4d_ode_poincare_invariant_1st(Δt, nloop, ntime, nsave, DT=Float64)
        guiding_center_4d_ode_init(q₀) = guiding_center_4d_ode(q₀; periodic=false)
        PoincareInvariant1st(guiding_center_4d_ode_init, f_loop, α, Δt, 4, nloop, ntime, nsave, DT)
    end

    function guiding_center_4d_iode_poincare_invariant_1st(Δt, nloop, ntime, nsave, DT=Float64)
        guiding_center_4d_iode_init(q₀) = guiding_center_4d_iode(q₀; periodic=false)
        PoincareInvariant1st(guiding_center_4d_iode_init, f_loop, α, Δt, 4, nloop, ntime, nsave, DT)
    end

end
