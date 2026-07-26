
# Guiding Center Dynamics in 4D

Guiding centre dynamics is a reduced version of charged particle dynamics, where the motion of the particle in a strong magnetic field ``B`` is reduced to the motion of the guiding centre, that is the centre of the gyro motion of the particle about a magnetic field line.
The dynamics of the guiding centre can be described in terms of only four coordinates (as compared to six for the full motion of the charged particle), the position of the guiding centre ``r = (x,y,z)`` and the parallel velocity ``u``, where parallel refers to the direction of the magnetic field.


## Lagrangian Formulation

The simplest form of the guiding centre equations can be obtained from the Lagrangian
```math
\begin{aligned}
L &= \left( A (r) + u b (r) \right) \cdot \dot{r} - H(r,u) , &
H &= \tfrac{1}{2} u^{2} + \mu \vert B (r) \vert ,
\end{aligned}
```
where ``b = B / \vert B \vert`` is the unit vector of the magnetic field ``B = \nabla \times A`` with ``A`` the magnetic vector potential and ``\mu`` is the magnetic moment.
The Euler-Lagrange equations are computed as
```math
\begin{aligned}
\nabla \vartheta^{T} ( r(t), u(t) ) \cdot \dot{r} (t) - \dot{\vartheta} ( r(t), u(t) ) &= \nabla H ( r(t), u(t) ) , \\
b ( r(t) ) \cdot \dot{r} (t) &= u (t) ,
\end{aligned}
```
with ``\vartheta(r,u) = A (r) + u \, b (r)`` and the gradient ``\nabla`` denoting the derivative with respect to ``r``.


## Hamiltonian Formulation

Computing the time derivative of ``\vartheta``, the Euler-Lagrange equations can be rewritten in an explicit form as
```math
\begin{aligned}
\dot{r} (t) &= \dfrac{u (t) \, \beta (r (t))}{b (r (t)) \cdot \beta (r (t))} + \dfrac{B (r (t))}{B (r (t)) \cdot \beta (r (t))} \times \nabla H (r (t),u (t)) , \\
\dot{u} (t) &= - \dfrac{\beta (r (t))}{b (r (t)) \cdot \beta (r (t))} \cdot \nabla H (r (t), u (t)) ,
\end{aligned}
```
where ``\beta = \nabla \times \vartheta``.
This constitutes a noncanonical Hamiltonian system of the form
```math
\dot{q} = \Omega^{-T} (q) \nabla H(q) ,
```
with ``q = (x,y,z,u)`` and the symplectic matrix ``\Omega`` given by
```math
\Omega_{ij} = \dfrac{\partial \vartheta_{i}}{\partial q^{j}} - \dfrac{\partial \vartheta_{j}}{\partial q^{i}} ,
```
which is the convention used throughout this package and shared with
[GeometricProblems.jl](https://github.com/JuliaGNI/GeometricProblems.jl).
With it the Euler-Lagrange equations of ``L = \vartheta \cdot \dot{q} - H`` read
``\Omega (q) \, \dot{q} = - \nabla H (q)``.


## Modules

```@docs
ChargedParticleDynamics.GuidingCenter4d
```

Each equilibrium is its own module, and every one of them `include`s the same
`guiding_center_4d_common.jl` and `guiding_center_4d_equations.jl` — plus
`guiding_center_4d_loop.jl` and `guiding_center_4d_surface.jl` where the equilibrium defines a
Poincaré loop or surface. The model's functions are therefore documented once, below, under
`TokamakSmallCylindrical`, and hold verbatim for all eleven. What differs between the modules is the
chart, the equilibrium parameters and the initial conditions, which is what these docstrings record:

```@docs
ChargedParticleDynamics.GuidingCenter4d.SymmetricField
ChargedParticleDynamics.GuidingCenter4d.ThetaPinchField
ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCartesian
ChargedParticleDynamics.GuidingCenter4d.TokamakSmallToroidal
ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCartesian
ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCylindrical
ChargedParticleDynamics.GuidingCenter4d.TokamakIterCylindrical
ChargedParticleDynamics.GuidingCenter4d.SolovevSymmetricField
ChargedParticleDynamics.GuidingCenter4d.SolovevIter
ChargedParticleDynamics.GuidingCenter4d.SolovevIterXpoint
```

```@autodocs
Modules = [ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCylindrical]
```
