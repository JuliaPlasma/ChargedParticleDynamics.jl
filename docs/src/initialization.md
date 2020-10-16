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



## Position

In order to compute the guiding centre position, we need to construct the gyro radius vector $\rho$, which is given by
```math
\rho = \frac{b \times v}{\omega_c} ,
```
then the guiding center position is
```math
X = x - \rho .
```
The gyro phase $\theta$ is then defined as the angle, measured in the clockwise sense, between $a$ and $\rho$, so that
```math
\frac{b \times v}{\vert v' \vert \, \sin \alpha} = a \, \cos \theta - c \, \sin \theta .
```
