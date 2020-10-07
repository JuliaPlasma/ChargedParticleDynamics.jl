
# Guiding Center Dynamics in 4D

Guiding centre dynamics is a reduced version of charged particle dynamics, where the motion of the particle in a strong magnetic field ``B`` is reduced to the motion of the guiding centre, that is the centre of the gyro motion of the particle about a magnetic field line.
The dynamics of the guiding centre can be described in terms of only four coordinates (as compared to six for the full motion of the charged particle), the position of the guiding centre ``r = (x,y,z)`` and the parallel velocity ``u``, where parallel refers to the direction of the magnetic field.


## Lagrangian Formulation

The simplest form of the guiding centre equations can be obtained from the Lagrangian
```math
\begin{aligned}
L &= (A (r) + u b (r)) \cdot \dot{r} - H(r,u) , &
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
\Omega_{ij} = \dfrac{\partial \vartheta_{j}}{\partial q^{i}} - \dfrac{\partial \vartheta_{i}}{\partial q^{j}} .
```


## Modules

```@autodocs
Modules = [ChargedParticleDynamics.GuidingCenter4d.SymmetricField,
           ChargedParticleDynamics.GuidingCenter4d.ThetaPinchField,
           ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCartesian,
           ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCylindrical,
           ChargedParticleDynamics.GuidingCenter4d.TokamakSmallToroidal,
           ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCartesian,
           ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCylindrical,
           ChargedParticleDynamics.GuidingCenter4d.TokamakIterCylindrical,
           ChargedParticleDynamics.GuidingCenter4d.SolovevIter,
           ChargedParticleDynamics.GuidingCenter4d.SolovevIterXpoint]
```
