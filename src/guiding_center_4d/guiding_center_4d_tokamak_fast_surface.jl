module TokamakFastSurface

    using GeometricIntegrators
    using MagneticEquilibria

    export guiding_center_4d_ode_poincare_invariant_2nd,
           guiding_center_4d_iode_poincare_invariant_2nd,
           hamiltonian, toroidal_momentum, u, β

    const R₀ = 2
    const B₀ = 5
    const q  = 2
    const μ  = 1E-3

    load_equilibrium(MagneticEquilibria.AxisymmetricTokamakCylindrical(R₀, B₀, q); target_module=TokamakFastSurface)

    include("guiding_center_4d_coords_RZphi.jl")


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


    function guiding_center_4d_ode_poincare_invariant_2nd(Δt, nx, ny, ntime, nsave, DT=Float64)
        guiding_center_4d_ode_init(q₀) = guiding_center_4d_ode(q₀; periodic=false)
        PoincareInvariant2nd(guiding_center_4d_ode_init, f_surface, ω, Δt, 4, nx, ny, ntime, nsave, DT)
    end

    function guiding_center_4d_iode_poincare_invariant_2nd(Δt, nx, ny, ntime, nsave, DT=Float64)
        guiding_center_4d_iode_init(q₀) = guiding_center_4d_iode(q₀; periodic=false)
        PoincareInvariant2nd(guiding_center_4d_iode_init, f_surface, ω, Δt, 4, nx, ny, ntime, nsave, DT)
    end

end
