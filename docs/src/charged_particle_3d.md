
# Charged Particles in 3D

The motion of charged particles in an electromagnetic field ``(E,B)`` is
governed by the Lorentz force,
```math
\ddot{x} (t) = \frac{e}{m} \big[ E(x(t)) + \dot{x} (t) \times B(x(t)) \big] ,
```
where ``m`` and ``e`` denote the particle's mass and charge, respectively.


## Canonical Formulation

The canonical form of the equations can be obtained from the Hamiltonian
```math
H (x,p) = \frac{1}{2m} (p-A(x))^2 + e \phi(x),
```
as
```math
\begin{aligned}
\dot{x} (t) &= \frac{\partial H}{\partial p} (x(t),p(t)) = \frac{1}{m} (p(t) - A(x(t))), \\
\dot{p} (t) &= - \frac{\partial H}{\partial x} (x(t),p(t)) = \frac{1}{m} \nabla A(x(t)) \cdot (p(t)-A(x(t))) - \nabla \phi(x(t)) ,
\end{aligned}
```
where the fields ``(E,B)`` are related to the potentials ``(\phi, A)`` by
```math
\begin{aligned}
E (x) &= - \nabla \phi (x) , &
B (x) &= \nabla \times A (x) .
\end{aligned}
```


## Noncanonical Formulation

The noncanonical form of the equations can be obtained from the phasespace Lagrangian
```math
\begin{aligned}
L (x,\dot{x},v,\dot{v}) &= (e A(x) + mv) \cdot \dot{x} - H(x,v) , &
H(x,v) &= \frac{m}{2} v^2 + e \phi(x),
\end{aligned}
```
as
```math
\begin{aligned}
\dot{x} (t) &= v (t) , &
\dot{v} (t) &= \frac{e}{m} \big[ \nabla A (x(t)) \cdot \dot{x}(t) - \dot{A} (x(t)) - \nabla \phi(x(t)) \big] .
\end{aligned}
```
Computing the time derivative of ``A`` and using the relation between the potentials ``(\phi, A)`` and the fields ``(E,B)``, this can be rewritten as
```math
\begin{aligned}
\dot{x} (t) &= v (t) , &
\dot{v} (t) &= \frac{e}{m} \big[ E(x(t)) + v (t) \times B(x(t)) \big] .
\end{aligned}
```
This constitutes a noncanonical Hamiltonian system of the form
```math
\dot{z} (t) = \Omega^{-T} (z(t)) \nabla H(z(t)) ,
```
with ``z = (x,v)`` and the symplectic matrix ``\Omega`` given by
```math
\Omega = \frac{1}{m} \begin{pmatrix}
  \mathbb{0} & \mathbb{1} \\
- \mathbb{1} & e \hat{B} \\
\end{pmatrix} ,
```
and
```math
\hat{B} = \begin{pmatrix}
0 & -B_3 & B_2 \\
B_3 & 0 & - B_1 \\
- B_2 & B_1 & 0 \\
\end{pmatrix} .
```


## Modules

```@autodocs
Modules = [ChargedParticleDynamics.ChargedParticle3d.SingularField,
           ChargedParticleDynamics.ChargedParticle3d.SymmetricField,
           ChargedParticleDynamics.ChargedParticle3d.ThetaPinchCanonical,
           ChargedParticleDynamics.ChargedParticle3d.ThetaPinchNoncanonical,
           ChargedParticleDynamics.ChargedParticle3d.TokamakSmallCartesian,
           ChargedParticleDynamics.ChargedParticle3d.TokamakSmallCylindrical,
           ChargedParticleDynamics.ChargedParticle3d.TokamakSmallToroidal,
           ChargedParticleDynamics.ChargedParticle3d.TokamakSmallNoncanonical,
           ChargedParticleDynamics.ChargedParticle3d.TokamakIterCylindrical,
           ChargedParticleDynamics.ChargedParticle3d.SolovevIter,
           ChargedParticleDynamics.ChargedParticle3d.SolovevIterXpoint]
```
