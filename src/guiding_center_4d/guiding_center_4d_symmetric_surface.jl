module SymmetricSurface

    using GeometricIntegrators
    using MagneticEquilibria

    export guiding_center_4d_ode_poincare_invariant_2nd,
           guiding_center_4d_iode_poincare_invariant_2nd,
           hamiltonian, toroidal_momentum, u, α, ω

    const B₀ = 1.
    const μ  = 1E-2

    load_equilibrium(SymmetricQuadratic(B₀); target_module=SymmetricSurface)

    include("guiding_center_4d_xyz.jl")


    function f_surface(s,t)
        z0 = 0.0
        u0 = 5E-1
        r0 = 0.1

        x  = 2r0*(s-0.5)
        y  = 2r0*(t-0.5)

        q  = [x, y, z0, u0]

        return q
    end

    function guiding_center_4d_ode_poincare(q₀)
        guiding_center_4d_ode(q₀; periodic=false)
    end

    function guiding_center_4d_iode_poincare(q₀)
        guiding_center_4d_iode(q₀; periodic=false)
    end
    function guiding_center_4d_ode_poincare_invariant_2nd(Δt, nx, ny, ntime, nsave, DT=Float64)
        PoincareInvariant2nd(guiding_center_4d_ode_poincare, f_surface, ω, Δt, 4, nx, ny, ntime, nsave, DT)
    end

    function guiding_center_4d_ode_poincare_invariant_2nd_trapezoidal(Δt, nx, ny, ntime, nsave, DT=Float64)
        PoincareInvariant2ndTrapezoidal(guiding_center_4d_ode_poincare, f_surface, ω, Δt, 4, nx, ny, ntime, nsave, DT)
    end

    function guiding_center_4d_iode_poincare_invariant_2nd(Δt, nx, ny, ntime, nsave, DT=Float64)
        PoincareInvariant2nd(guiding_center_4d_iode_poincare, f_surface, ω, Δt, 4, nx, ny, ntime, nsave, DT)
    end

    function guiding_center_4d_iode_poincare_invariant_2nd_trapezoidal(Δt, nx, ny, ntime, nsave, DT=Float64)
        PoincareInvariant2ndTrapezoidal(guiding_center_4d_iode_poincare, f_surface, ω, Δt, 4, nx, ny, ntime, nsave, DT)
    end

end
