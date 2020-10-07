# Normalization

In order to normalize the equations of the various models implemented in this package, we start at
the level of the Lagrangian. In the following, this is explained exemplary for the charged particle
Lagrangian, however, it generalizes straightforwardly to other systems like the Pauli particle or
the guiding center system.


## Charged Particle Lagrangian

Consider the phasespace Lagrangian
```math
L (q, \dot{q}, v) = ( mv + A (x) ) \cdot \dot{x}  - \frac{m}{2} \vert v \vert^2  - e \phi (x) ,
```

and, in full generality, introduce the following normalizations:
```math
t = \hat{t} t' , \quad
x = \hat{x} x' , \quad
v = \hat{v} v' , \quad
A = \hat{A} A' , \quad
\phi = \hat{\phi} \phi' , \quad
L = \hat{L} L ' .
```

This leads us to
```math
L'
= \frac{L}{\hat{L}}
= \frac{1}{\hat{L}} \, \big( m \hat{v} v' + e \hat{A} A' \big) \cdot \frac{\hat{x}}{\hat{t}} \dot{x}' - \frac{m \hat{v}}{\hat{L}} \frac{\vert v' \vert^2}{2} - \frac{e \hat{\phi}}{\hat{L}} \phi' .
```

Let us choose the following normalizations (and note that others are possible and may be more appropriate, depending on the problem at hand):
```math
\begin{aligned}
\hat{A} &= \hat{l} \hat{B} , &
\hat{\phi} &= \hat{l} \hat{E} , &
\hat{L} &= m \hat{v}^2 = e \hat{\phi} ,
\end{aligned}
```
with the characteristic length $l$.
The normalized Lagrangian becomes
```math
L' = \bigg( \frac{\hat{x}}{\hat{t} \hat{v}} \, v' + \underbrace{\frac{e \hat{B}}{m}}_{\hat{\omega}_c} \frac{\hat{x} \hat{l}}{ \hat{t} \hat{v}^2} \ A' \bigg) \cdot \dot{x}' - \frac{\vert v' \vert^2}{2} - \frac{\hat{l}}{\hat{x}} \, \phi' ,
```
where $\omega_c = e B / m$ is the gyration frequency.
This suggests to set
```math
\begin{aligned}
\hat{t} &= \hat{\omega}_c^{-1} , &
\hat{x} &= \hat{l} , &
\hat{v} &= \frac{\hat{x}}{\hat{t}} = \hat{l} \hat{\omega}_c ,
\end{aligned}
```
so that
```math
L' = ( v' + A' ) \cdot \dot{x}' - \frac{\vert v' \vert^2}{2} - \phi' .
```
We thus have obtained the normalized Lagrangian.


## Alternative Normalization

Often, especially when simulating an ensemble of particles, it is more appropriate to normalize the velocity to the thermal velocity and to choose different normalizations for $v$ and $\dot{x}$, specifically
```math
\begin{aligned}
\hat{v} &= v_{\mathrm{th}} , &
\hat{t} &= \omega_c^{-1} , &
\hat{A} &= \hat{l} \hat{B} , &
\hat{\phi} &= \hat{l} \hat{E} , &
\hat{L} &= m \hat{v}^2 = e \hat{\phi} = \hat{W} ,
\end{aligned}
```
where $W = k_B T$ denotes the thermal energy and $v_{\mathrm{th}} = \sqrt{ 2 W / m}$.
Consequently,
```math
\begin{aligned}
\hat{x} &= \hat{v} \hat{t} = \frac{v_{\mathrm{th}}}{\omega_c} = \rho_{\mathrm{th}} ,
\end{aligned}
```
and the normalized Lagrangian becomes
```math
L' = \bigg( v' + \frac{\hat{l}}{\rho_{\mathrm{th}}} \, A' \bigg) \cdot \dot{x}' - \frac{\vert v' \vert^2}{2} - \frac{\hat{l}}{\rho_{\mathrm{th}}} \, \phi' .
```
