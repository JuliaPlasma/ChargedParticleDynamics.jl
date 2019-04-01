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

    using ElectromagneticFields

    const μ  = 1E-2

    load_equilibrium(SymmetricQuadratic(1.); target_module=SymmetricLoop)

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

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")
    include("guiding_center_4d_loop.jl")

end


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

    using ElectromagneticFields

    const μ  = 1E-2

    load_equilibrium(SymmetricQuadratic(1.); target_module=SymmetricSurface)

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

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")
    include("guiding_center_4d_surface.jl")

end
