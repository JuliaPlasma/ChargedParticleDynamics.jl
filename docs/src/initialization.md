# Initialization

In the following, we will discuss how to initialise the various models from a common set of initial data, namely
- the guiding center position $X$ in physical coordinates,
- the particle energy $W$ in eV,
- the particle mass $m$ in kg,
- the particle charge number,
- the gyro angle $\theta \in [0, 2\pi]$,
- the pitch angle $\alpha \in [0, \pi / 2]$.

For charged particles we need to convert the energy into a velocity vector.
For the Pauli particle we need to compute the coordinates of the guiding centre $X$ of the particle and construct a velocity vector that only holds the component of the velocity that is parallel to the magnetic field, while the perpendicular component goes into the magnetic moment $\mu$.
For the guiding center, we need the same coordinate transformation as for the Pauli particle, but we only need the absolute value of the velocity in direction of the magnetic field $u$ and the magnetic moment $\mu$.

Make sure to read about the [Normalization](@ref) before proceeding.


## Velocities

Set $e = e' \hat{e}$ with $e' = 1.602 176 634 \cdot 10^{-19}$ and $\hat{e} = \mathrm{C}$ and similarly $m = m' \hat{m}$ with $\hat{m} = \mathrm{kg}$. Some common values for the mass are

| Particle   | Normalised Mass $m'$          |
| ---------- | ----------------------------- |
| electron   | $9.1093837015 \cdot 10^{-31}$ |
| proton     | $1.6726219237 \cdot 10^{-27}$ |
| deuteron   | $3.3435837724 \cdot 10^{-27}$ |
| alpha      | $6.6446573357 \cdot 10^{-27}$ |

The absolute value of the velocity is obtained from the particle energy by $\vert v \vert^2 = 2 W / m$.
In normalised units, we have $W = W' \hat{W}$ and $v = v' \hat{v}$ where $\hat{W} = e \mathrm{V}$ and $\hat{v} = \hat{l} \hat{\omega}_c = \hat{l} \hat{B} e / m$.

Let us compute the normalised kinetic energy,
```math
\frac{\vert v' \vert^2}{2}
= \frac{W' \hat{W}}{m \hat{v}^2}
= \frac{W' e \mathrm{V} m^2}{m \hat{l}^2 \hat{B}^2 e^2} .
```

Usually, we have $\hat{B} = T = \mathrm{kg} / \mathrm{C}\,\mathrm{s}$. Recall that $\mathrm{V} = \mathrm{kg}\,\mathrm{m}^2 / \mathrm{C}\,\mathrm{s}^2$ and let us set $\hat{l} = l_0 \mathrm{m}$, then
```math
\frac{\vert v' \vert^2}{2} = \frac{m' W'}{e' l_0^2} .
```

Thus for the absolute value of the velocity we have
```math
\vert v' \vert = \frac{1}{l_0} \sqrt{ 2 W' \, \frac{m'}{e'} } .
```

The pitch angle $\alpha$ determines the distribution of the kinetic energy into the perpendicular and parallel components by
```math
\begin{aligned}
W_{\perp}' &= W' \sin \alpha , \qquad
W_{\parallel}' &= W' \, (1 - \sin \alpha) .
\end{aligned}
```
Note that this is not the textbook pitch-angle convention
``W_{\perp} = W \sin^{2} \alpha``, ``W_{\parallel} = W \cos^{2} \alpha``; it is however
self-consistent with the inverse relation ``\alpha = \arcsin (W_{\perp}' / W')`` used by
`InitialConditionsGC` to recover the pitch angle from ``(u, \mu)``.
and accordingly
```math
\begin{aligned}
\vert v_{\perp}' \vert &= \frac{1}{l_0} \sqrt{ 2 W_{\perp}' \, \frac{m'}{e'} } , \qquad
\vert v_{\parallel}' \vert &= \frac{1}{l_0} \sqrt{ 2 W_{\parallel}' \, \frac{m'}{e'} } .
\end{aligned}
```

For each equilibrium, the ElectromagneticFields.jl package generates an orthonormal triad
$(a, b, c)$, where $b$ is the unit vector of the magnetic field and $a$ and $c$ span the plane
perpendicular to it.
The triad is exposed in three different component representations, which coincide only in
Cartesian coordinates:

| Functions | Components |
| ------------------------------ | ---------------- |
| `a(t,x)`, `b(t,x)`, `c(t,x)` | covariant |
| `aₚ(t,x)`, `bₚ(t,x)`, `cₚ(t,x)` | physical |
| `a⃗(t,x)`, `b⃗(t,x)`, `c⃗(t,x)` | contravariant |

The velocity vector $v$ is assembled in physical coordinates, that is from `aₚ`, `bₚ` and `cₚ`, and
subsequently transformed to contravariant coordinates with the inverse Jacobian `DF̄`.
The function `bₚ` returns the unit vector of the magnetic field, thus
```math
v_{\parallel}' = \vert v' \vert \, \sqrt{1 - \sin \alpha} \, b .
```
The functions `aₚ` and `cₚ` return two unit vectors that span the plane perpendicular to the magnetic field, thus
```math
v_{\perp}' = \vert v' \vert \, \sqrt{\sin \alpha} \, ( - a \, \sin \theta - c \, \cos \theta ) ,
```
where $\theta$ is the gyro angle.
The magnetic moment is computed as
```math
\mu = \frac{ \vert v_{\perp}' \vert^2 }{2 \vert B' \vert} = \frac{ \vert v' \vert^2 \sin \alpha }{2 \vert B' \vert} .
```


## Gyro-radius Vector and Position

In order to compute the particle position, we need to construct the gyro radius vector $\rho$, which is given by
```math
\rho = \frac{b \times v}{\omega_c} = \frac{b \times \hat{v} v'}{\hat{\omega}_c \omega_c'}  = \hat{l} \, \frac{b \times v'}{\vert B' \vert} .
```
Here we used $\hat{v} / \hat{\omega}_c = \hat{l}$, which follows from $\hat{v} = \hat{l} \hat{\omega}_c$, and recall that $\hat{l} = l_0 \mathrm{m}$.
Note that
```math
\omega_c
= \frac{e \vert B \vert}{m}
= \frac{e \vert \hat{B} B' \vert}{m}
= \frac{e \vert \hat{B} \vert}{m} \vert B' \vert
= \hat{\omega}_c \omega_c' ,
```
as
```math
\hat{\omega}_c = \frac{e \vert \hat{B} \vert}{m} ,
```
so that
```math
\omega_c' = \vert B' \vert .
```
The normalized gyro radius vector reads
```math
\rho' = \frac{\rho}{\hat{l}} = \frac{b \times v'}{\vert B' \vert} ,
```
so that the normalised particle position is
```math
x' = X' + \rho' .
```
The gyro phase $\theta$ is defined as the angle, measured in the clockwise sense, between $a$ and $\rho$, so that
```math
\frac{b \times v'}{\vert v' \vert \, \sin \alpha} = a \, \cos \theta - c \, \sin \theta .
```


## Guiding Center Coordinates

Instead of the particle energy $W$ and pitch angle $\alpha$, we can also compute all initial conditions from the guiding center initial data, i.e., parallel velocity $u = \vert v_{\parallel} \vert$, magnetic moment $\mu$ and gyro angle $\theta$.

Start by computing the absolute value of the perpendicular velocity by
```math
\vert v_\perp \vert = \sqrt{ 2 \mu \vert B \vert } ,
```
and the total velocity by
```math
\vert v \vert = \sqrt{ \vert v_{\parallel} \vert^2 + \vert v_\perp \vert^2 } .
```

With that, we can compute the perpendicular and total energy,
```math
\begin{aligned}
W_{\perp}' &= \frac{e l_0^2 \vert v_\perp \vert^2}{2m} , \qquad
W' &= \frac{e l_0^2 \vert v \vert^2}{2m} .
\end{aligned}
```
The pitch angle is obtained as
```math
\alpha = \arcsin \frac{W_{\perp}'}{W'} .
```
With that, the rest of the quantities can be computed in the same way as above.


## Initial Conditions

ChargedParticleDynamics.jl provides the `InitialConditions` module for the computation of the above quantities.
Three functions are provided that return the initial conditions for the different models.
For charged particles the initial conditions are $(x', v')$, for the Pauli particle we have $(X', v_{\parallel}', \mu)$ and for the guiding center we use $(X', u', \mu)$ with $u' = \vert v_{\parallel}' \vert$.

!!! warning "The equilibrium modules do not use this"
    Everything on this page describes machinery that is available but not wired up. Each
    equilibrium module carries its `u` and `μ` as hard-coded literals rather than deriving them
    from a position, a pitch angle and an energy through `InitialConditions`. The only callers of
    these functions are this page, the [ITER Equilibrium in Cylindrical Coordinates](@ref) example
    and `scripts/guiding_center_3d.jl`.

    The literals are also not reproducible from the functions here. For the small tokamak in
    cartesian coordinates, whose Pauli module declares `xᵢ = [1.05, 0, 0]` and
    `vᵢ = [2.1e-3, 4.3e-4, 0]`, the parallel velocity ``v \cdot b`` comes out at `0.00042987`
    against a hard-coded `0.00045136`, and ``\vert v_{\perp} \vert^{2} / 2B`` at `2.31459e-6`
    against a hard-coded `2.310e-6` — the latter being a nearby but *different* number rather than
    the same one rounded. Where the shipped values came from is not recorded anywhere in the
    repository, so the per-equilibrium docstrings say so rather than inventing a derivation.

    Routing the equilibrium modules through this layer is the fix, and is the task *One
    initial-conditions layer for all models* in `TODO.md`. It would make `u` and `μ` derived
    quantities that a new equilibrium gets for free, at the cost of shifting the shipped values by
    the percentages above.

```@autodocs
Modules = [ChargedParticleDynamics]
Order   = [:type, :function]
```


## Example

As an example, let us consider a deuteron in an ITER-like analytical equilibrium (obtained from `ElectromagneticFields.SolovevITER`). The guiding center position is $[7, 0, 0]$, the energy is $1 \, \mathrm{MeV}$, and the pitch angle is $\pi / 2$.

```@eval
using CairoMakie
using LaTeXStrings
using Markdown

using ChargedParticleDynamics
using ChargedParticleDynamics: md

import ElectromagneticFields: code
import ElectromagneticFields.Solovev: SolovevEquilibrium
import ElectromagneticFields.SolovevITER: init
equ = init()
@eval $(code(equ))

X₀ = from_cartesian(0, [7.0, 0, 0])
E₀ = 1E6
θ₀ = 0.
α₀ = π/4
m₀ = md

ic0 = InitialConditions(X₀, θ₀, α₀, E₀, m₀, 1, aₚ, bₚ, cₚ, b⃗, B, ḡ, DF̄, J; l₀=R₀)

np = 12
nr = 100
nz = 120

xgrid = LinRange( 0.5,   1.5, nr)
ygrid = LinRange(-0.75, +0.75, nz)
# the poloidal flux is the covariant component A₃ = R A_φ = ψ itself, not A₃ / R — see the note on
# `plot_fieldlines` in the audit
fieldlines = [A₃(0, x, y, 0.0) for x in xgrid, y in ygrid]

px = zeros(np)
py = zeros(np)

for i in 1:np
    ics = charged_particle(InitialConditions(X₀, 2π*i/np, α₀, E₀, m₀, 1, aₚ, bₚ, cₚ, b⃗, B, ḡ, DF̄, J; l₀=R₀))
    px[i] = ics[1][1]
    py[i] = ics[1][2]
end

boxw = 0.008
boxh = 0.008

fg = Figure(size=(800, 400))

axl = Axis(fg[1, 1], aspect=1, xlabel=L"R/R_0", ylabel=L"Z/R_0",
           limits=((0.5, 1.5), (-0.75, +0.75)))
axr = Axis(fg[1, 2], aspect=1, xlabel=L"R/R_0", ylabel=L"Z/R_0",
           limits=(X₀[1] .+ (-boxw, +boxw), X₀[2] .+ (-boxh, +boxh)))

contour!(axl, xgrid, ygrid, fieldlines, levels=50, colormap=:reds)
poly!(axl, Rect(X₀[1] - 0.05, X₀[2] - 0.05, 0.1, 0.1), color=(:white, 0.5), strokecolor=:black, strokewidth=1)

for ax in (axl, axr)
    scatter!(ax, px, py, color=:blue)
    scatter!(ax, [X₀[1]], [X₀[2]], color=:red)
end

save("initial_conditions.png", fg)

Markdown.parse("""
| Fields     | Value |
|:-----------|:------|
| `x`   | $(ic0.x)  |
| `X`   | $(ic0.X)  |
| `ρ`   | $(ic0.ρ)  |
| `v⃗`   | $(ic0.vvec)  |
| `v∥`  | $(ic0.vpar)  |
| `v⟂`  | $(ic0.vper)  |
| `v`   | $(ic0.v)  |
| `u`   | $(ic0.u)  |
| `μ`   | $(ic0.μ)  |
| `θ`   | $(ic0.θ)  |
| `α`   | $(ic0.α)  |
| `ω`   | $(ic0.ω)  |
| `M`   | $(ic0.mass)  |
| `E`   | $(ic0.energy)  |
| `C`   | $(ic0.charge)  |
""")
```

Below we plot the guiding center position (in red) and the particle position (in blue) in the poloidal plane for various gyro angles $\in [0, 2\pi]$. The right panel zooms into the white box around the guiding center.

![](initial_conditions.png)
