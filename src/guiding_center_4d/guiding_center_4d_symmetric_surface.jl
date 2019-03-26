"""
Second Poincaré invariant for a guiding center particle in an axisymmetric
magnetic field of the form ``B(x,y,z) = B_0 (1 + x^2 + y^2) e_z``.

The surface for the invariant is initialized by
```math
q (\\tau) = \\begin{pmatrix}
r_0 (\\sigma - 0.5) \\\\
r_0 (\\tau   - 0.5) \\\\
z_0 + z_1 \\cos (2\\pi \\sigma) \\cos (2\\pi \\tau) \\\\
u_0 + u_1 \\sin (2\\pi \\sigma) \\sin (2\\pi \\tau) \\\\
\\end{pmatrix}
```
with parameters
```math
B_0 = 1, \\quad
r_0 = 0.5, \\quad
z_0 = 0.0, \\quad
z_1 = 0.1, \\quad
u_0 = 0.5, \\quad
u_1 = 0.01, \\quad
\\mu = 0.01 .
```
"""
module SymmetricSurface

    using GeometricIntegrators
    using GeometricIntegrators.Utils
    using ElectromagneticFields

    export guiding_center_4d_ode_poincare_invariant_2nd,
           guiding_center_4d_iode_poincare_invariant_2nd,
           hamiltonian, toroidal_momentum, u, α, ω

    const B₀ = 1.
    const μ  = 1E-2

    load_equilibrium(SymmetricQuadratic(B₀); target_module=SymmetricSurface)

    include("guiding_center_4d_coords_xyz.jl")


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
