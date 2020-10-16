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


Set $e = e' \hat{e}$ with $e' = 1.602 176 634 \cdot 10^{-19}$ and $\hat{e} = \mathrm{C}$ and similarly $m = m' \hat{m}$ with $\hat{m} = \mathrm{kg}$. Some common values for the mass are

| Particle   | Normalised Mass $m'$          |
| ---------- | ----------------------------- |
| electron   | $9.1093837015 \cdot 10^{-31}$ |
| proton     | $1.6726219237 \cdot 10^{-27}$ |
| deuteron   | $3.3435837724 \cdot 10^{-27}$ |
| alpha      | $6.6446573357 \cdot 10^{-27}$ |

The absolute value of the velocity is obtained from the particle energy by $\vert v \vert^2 = 2 W / m$.
In normalised units, we have $W = W' \hat{W}$ and $v = v' \hat{v}$ where $\hat{W} = e \mathrm{V}$ and $\hat{v} = \hat{l} \hat{\omega}_c = \hat{l} \hat{B} e / m$.

Let us compute the normalise kinetic energy,
```math
\frac{\vert v' \vert^2}{2}
= \frac{W' \hat{W}}{m \hat{v}^2}
= \frac{W' e \mathrm{V} m^2}{m \hat{l}^2 \hat{B}^2 e^2} .
```

Usually, we have $\hat{B} = T = \mathrm{kg} / \mathrm{C} \mathrm{s}$. Recall that $\mathrm{V} = \mathrm{kg} \mathrm{m}^2 / \mathrm{C} \mathrm{s}^2$ and let us set $\hat{l} = l_0 \mathrm{m}$, then
```math
\frac{\vert v' \vert^2}{2}
= \frac{m' W'}{e' l_0^2} .
```

Thus for the absolute value of the velocity we have
```math
\vert v' \vert
= \frac{1}{l_0} \sqrt{ 2 W' \, \frac{m'}{e'} } .
```

The pitch angle $\alpha$ determines the distribution of the kinetic energy into the perpendicular and parallel components by
```math
\begin{aligned}
\vert v_{\perp}' \vert &= \vert v' \vert \sin \alpha , \\
\vert v_{\parallel}' \vert &= \vert v' \vert - \vert v_{\perp}' \vert = \vert v' \vert \, (1 - \sin \alpha) .
\end{aligned}
```

The ElectromagneticFields.jl package provides three functions `a⃗(t,x)`, `b⃗(t,x)`, `c⃗(t,x)` that can be used to construct the velocity vector $v$.
The function `b⃗` returns the unit vector of the magnetic field, thus
```math
v_{\parallel}' = \vert v' \vert \, (1 - \sin \alpha) \, b .
```
The functions `a⃗` and `c⃗` return two unit vectors that span the plane perpendicular to the magnetic field, thus
```math
v_{\perp}' = \vert v' \vert \, \sin \alpha ( - a \, \sin \theta - c \, \cos \theta ) ,
```
where $\theta$ is the gyro angle.
The magnetic moment is computed as
```math
\mu = \frac{ \vert v_{\perp}' \vert^2 }{2 \vert B' \vert} = \frac{ \vert v' \vert^2 \sin \alpha }{2 \vert B' \vert} .
```

In order to compute the particle position, we need to construct the gyro radius vector $\rho$, which is given by
```math
\rho = \frac{b \times v}{\omega_c} = \frac{b \times \hat{v} v'}{\omega_c}  = l_{0} \, b \times v' ,
```
and the normalized gyro radius vector by
```math
\rho' = \frac{\rho}{l_0} = b \times v' ,
```
so that the normalised particle position is
```math
x' = X' + \rho' .
```
The gyro phase $\theta$ is defined as the angle, measured in the clockwise sense, between $a$ and $\rho$, so that
```math
\frac{b \times v'}{\vert v' \vert \, \sin \alpha} = a \, \cos \theta - c \, \sin \theta .
```


## Initial Conditions

ChargedParticleDynamics.jl provides the `InitialConditions` module for the computation of the above quantities.
Three functions are provided that return the initial conditions for the different models.
For charged particles the initial conditions are $(x', v')$, for the Pauli particle we have $(X', v_{\parallel}', \mu)$ and for the guiding center we use $(X', u', \mu)$ with $u' = \vert v_{\parallel}' \vert$.

```@autodocs
Modules = [ChargedParticleDynamics]
Order   = [:type, :function]
```


## Example

As an example, let us consider a deuteron in an ITER-like analytical equilibrium (obtained from `ElectromagneticFields.SolovevITER`). The guiding center position is $[7, 0, 0]$, the energy is $1 \, \mathrm{MeV}$, and the pitch angle is $\pi / 2$.

```@eval
using LaTeXStrings
using Markdown
using Plots

p = palette(:tab10)
rectangle(x, y, w, h) = Shape(x .+ [-w/2,-w/2,+w/2,+w/2], y .+ [-h/2,+h/2,+h/2,-h/2])

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
α₀ = π/2
m₀ = md

ic0 = InitialConditions(X₀, θ₀, α₀, E₀, m₀, 1, a⃗, b⃗, c⃗, B; l=R₀)

np = 12
nr = 100
nz = 120

xgrid = LinRange( 0.5,   1.5, nr)
ygrid = LinRange(-0.75, +0.75, nz)
fieldlines = zeros((nr,nz))

for i in 1:nr
    for j in 1:nz
        fieldlines[i,j] = A₃(0, xgrid[i], ygrid[j], 0.0) / xgrid[i]
    end
end

px = zeros(np)
py = zeros(np)

for i in 1:np
    ics = charged_particle(InitialConditions(X₀, 2π*i/np, α₀, E₀, m₀, 1, a⃗, b⃗, c⃗, B; l=R₀))
    px[i] = ics[1][1]
    py[i] = ics[1][2]
end

boxw=0.001
boxh=0.001

plot(aspectratio=1, layout=(1,2), legend=:none, size=(800,600), xguide=L"R/R_0", yguide=L"Z/R_0")
plot!(subplot=1, rectangle(X₀[1], X₀[2], 0.1, 0.1), opacity=.5, color=:white)

plot!(subplot=1, xlims=(0.5,1.5), ylims=(-0.75,+0.75))
plot!(subplot=2, xlims=(X₀[1] .+ [-boxw,+boxw]), ylims=(X₀[2] .+ [-boxh,+boxh]))

contour!(subplot=1, xgrid, ygrid, fieldlines', levels=50)
scatter!(subplot=1, [X₀[1]], [X₀[2]], color=p[4])

scatter!(subplot=2, px, py, color=p[1])
scatter!(subplot=2, [X₀[1]], [X₀[2]], color=p[4])

savefig("initial_conditions.png")

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

Below we plot the guiding center position (in blue) and the particle position (in red) in the poloidal plane for various gyro angles $\in [0, 2\pi]$.

![](initial_conditions.png)
