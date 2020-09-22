@doc raw"""
First Poincaré invariant for a guiding center particle in an θ-pinch magnetic
field of the form ``B(x,y,z) = B_0 \, e_z``.

The loop for the first Poincaré invariant is initialized by
```math
q (\tau) = \begin{pmatrix}
r_x \cos (2\pi \tau) \\
y_0 \\
r_z \sin (2\pi \tau) \\
\end{pmatrix}
```
with parameters
```math
B_0 = 1, \quad
r_x = 0.5, \quad
r_z = 0.3, \quad
y_0 = 0.0, \quad
u_0 = 0.5, \quad
\mu = 2.5 \times 10^{-6} .
```
"""
module GuidingCenter4dThetaPinch

    using ElectromagneticFields.ThetaPinch

    const equ = ThetaPinch.init(1.)

    μ_loop() = 2.5E-6

    function f_loop(t)
        μ  = 2.5E-6
        Y0 = 0.0
        u0 = 4E-4
        r0 = 0.5
        r1 = 0.3

        Xt = r0*cos(2π*t)
        Zt = r1*sin(2π*t)

        qt = [Xt, Y0, Zt, u0]

        return qt
    end

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")
    include("guiding_center_4d_loop.jl")
    include("guiding_center_4d_diagnostics.jl")

end
