"""
First Poincaré invariant for a guiding center particle in an axisymmetric
magnetic field of the form ``B(x,y,z) = B_0 (1 + x^2 + y^2) e_z``.

The loop for the invariant is initialized by
```math
q (\\tau) = \\begin{pmatrix}
r_x \\cos (2\\pi \\tau) \\\\
r_y \\sin (2\\pi \\tau) \\\\
z_0 + z_1 \\sin (2\\pi \\tau) \\\\
u_0 + u_1 \\cos (2\\pi \\tau) \\\\
\\end{pmatrix}
```
with parameters
```math
B_0 = 1, \\quad
r_x = 0.5, \\quad
r_y = 0.3, \\quad
z_0 = 0.0, \\quad
z_1 = 0.1, \\quad
u_0 = 0.5, \\quad
u_1 = 0.05, \\quad
\\mu = 0.01 .
```
"""
module SymmetricLoop

    using GeometricIntegrators
    using ElectromagneticFields

    export guiding_center_4d_ode_poincare_invariant_1st,
           guiding_center_4d_iode_poincare_invariant_1st,
           hamiltonian, toroidal_momentum, u, α

    const B₀ = 1.
    const μ  = 1E-2

    load_equilibrium(SymmetricQuadratic(B₀); target_module=SymmetricLoop)

    include("guiding_center_4d_coords_xyz.jl")


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
