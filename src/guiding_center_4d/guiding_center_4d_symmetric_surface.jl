module SymmetricSurface

    using GeometricIntegrators
    using GeometricIntegrators.Utils
    using MagneticEquilibria

    export guiding_center_4d_ode_poincare_invariant_2nd,
           guiding_center_4d_iode_poincare_invariant_2nd,
           hamiltonian, toroidal_momentum, u, α, ω

    const B₀ = 1.
    const μ  = 1E-2

    load_equilibrium(SymmetricQuadratic(B₀); target_module=SymmetricSurface)

    include("guiding_center_4d_xyz.jl")


    function f_surface(s,t)
        r0 = 0.5
        z0 = 0.0
        z1 = 0.1
        u0 = 5E-1
        u1 = 1E-2

        x  = r0*(s-0.5)
        y  = r0*(t-0.5)
        z  = z0 + z1 * cos(2π*s) * cos(2π*t)
        u  = u0 + u1 * sin(2π*s) * sin(2π*t)

        q  = [x, y, z, u]

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
        PoincareInvariant2nd(guiding_center_4d_iode_poincare, f_surface, ω, (D²ϑ₁, D²ϑ₂, D²ϑ₃, D²ϑ₄), Δt, 4, nx, ny, ntime, nsave, DT)
    end

    function guiding_center_4d_iode_poincare_invariant_2nd_trapezoidal(Δt, nx, ny, ntime, nsave, DT=Float64)
        PoincareInvariant2ndTrapezoidal(guiding_center_4d_iode_poincare, f_surface, ω, (D²ϑ₁, D²ϑ₂, D²ϑ₃, D²ϑ₄), Δt, 4, nx, ny, ntime, nsave)
    end

end
