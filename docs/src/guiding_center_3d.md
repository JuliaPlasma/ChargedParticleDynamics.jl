
# Guiding Center Dynamics in 3D

The [four-dimensional guiding centre model](guiding_center_4d.md) evolves the guiding centre position ``r = (x,y,z)`` together with the parallel velocity ``u``.
The three-dimensional model instead treats the parallel velocity as a dependent variable and evolves only the position together with its canonically conjugate momentum, so that the guiding centre dynamics takes the form of a *constrained canonical* Hamiltonian system.
The implementation follows Xinjie Li, Ruili Zhang and Jian Liu, *Approximately symplectic
Runge-Kutta methods for the guiding center dynamics*, Journal of Computational Physics **563**,
115069 (2026), [doi:10.1016/j.jcp.2026.115069](https://doi.org/10.1016/j.jcp.2026.115069), and its
companion Ruili Zhang, Jian Liu, Tong Liu, Wenxiang Li and Xiaogang Wang, *Canonical Hamiltonian
guiding center theory and classical intrinsic magnetic moment*, Frontiers of Physics **21**(2),
026200 (2026).


## Canonical Formulation

The starting point is the guiding centre one-form ``\vartheta(r,u) = A (r) + u \, b (r)``, where ``A`` is the magnetic vector potential, ``B = \nabla \times A`` and ``b = B / \vert B \vert``.
Identifying the momentum with the one-form,
```math
p = A (r) + u \, b (r) ,
```
the parallel velocity is recovered from ``p`` and ``r`` as
```math
v (r,p) = p - A (r) ,
\qquad
u (r,p) = b (r) \cdot v (r,p) ,
```
and the Hamiltonian becomes
```math
H (r,p) = \tfrac{1}{2} \, \vert v (r,p) \vert^{2} + \mu \, \vert B (r) \vert + \varphi (r) ,
```
with the magnetic moment ``\mu`` and the electrostatic potential ``\varphi``.

By construction ``v`` is parallel to ``b``, which is not preserved by the unconstrained canonical flow.
The two independent components of ``v \times b = 0`` are therefore imposed as constraints,
```math
\begin{aligned}
g_{1} (r,p) &= b_{1} (r) \, v_{2} (r,p) - b_{2} (r) \, v_{1} (r,p) = 0 , \\
g_{2} (r,p) &= b_{1} (r) \, v_{3} (r,p) - b_{3} (r) \, v_{1} (r,p) = 0 ,
\end{aligned}
```
so that the dynamics is governed by
```math
\begin{aligned}
\dot{r} (t) &= + \dfrac{\partial H}{\partial p} + \lambda_{1} \dfrac{\partial g_{1}}{\partial p} + \lambda_{2} \dfrac{\partial g_{2}}{\partial p} , \\
\dot{p} (t) &= - \dfrac{\partial H}{\partial r} - \lambda_{1} \dfrac{\partial g_{1}}{\partial r} - \lambda_{2} \dfrac{\partial g_{2}}{\partial r} ,
\end{aligned}
```
where the Lagrange multipliers ``\lambda_{1}`` and ``\lambda_{2}`` are determined algebraically by requiring ``\dot{g}_{1} = \dot{g}_{2} = 0``.
Since the multipliers are explicit functions of ``(r,p)``, this is a Hamiltonian system in the sense of a [`HODEProblem`](https://juliagni.github.io/GeometricEquations.jl/stable/), constructed by `hodeproblem`.


## Canonicalised Formulation

The multipliers can alternatively be absorbed into the Hamiltonian,
```math
\bar{H} (r,p) = H (r,p) + \lambda_{1} (r,p) \, g_{1} (r,p) + \lambda_{2} (r,p) \, g_{2} (r,p) ,
```
which agrees with ``H`` on the constraint manifold ``g_{1} = g_{2} = 0``.
Differentiating ``\bar{H}`` produces additional terms proportional to the constraints,
```math
\begin{aligned}
\dot{r} (t) &= \dfrac{\partial H}{\partial p} + \lambda_{1} \dfrac{\partial g_{1}}{\partial p} + \lambda_{2} \dfrac{\partial g_{2}}{\partial p} + \dfrac{\partial \lambda_{1}}{\partial p} g_{1} + \dfrac{\partial \lambda_{2}}{\partial p} g_{2} , \\
\dot{p} (t) &= - \dfrac{\partial H}{\partial r} - \lambda_{1} \dfrac{\partial g_{1}}{\partial r} - \lambda_{2} \dfrac{\partial g_{2}}{\partial r} - \dfrac{\partial \lambda_{1}}{\partial r} g_{1} - \dfrac{\partial \lambda_{2}}{\partial r} g_{2} ,
\end{aligned}
```
which vanish for exact initial data but not for a numerical solution that drifts off the constraint manifold.
This genuinely canonical system is constructed by `hodeproblem_canonical`, whose Hamiltonian is ``\bar{H}``; it requires the second derivatives of the magnetic field and is correspondingly more expensive to evaluate.


## Usage

Each equilibrium is wrapped in its own module, which injects the magnetic field code and provides
`initial_conditions` for converting ``(r,u)`` into ``(r,p)`` as well as `hodeproblem` and `hodeproblem_canonical`:

```julia
using GeometricIntegrators
using ChargedParticleDynamics.GuidingCenter3d.Dipole3d

problem = hodeproblem(initial_conditions_dipole())
solution = integrate(problem, PartitionedGauss(1))
```


## Modules

```@autodocs
Modules = [ChargedParticleDynamics.GuidingCenter3d.Dipole3d,
           ChargedParticleDynamics.GuidingCenter3d.QuadraticPotentials3d,
           ChargedParticleDynamics.GuidingCenter3d.SymmetricField,
           ChargedParticleDynamics.GuidingCenter3d.ThetaPinchField,
           ChargedParticleDynamics.GuidingCenter3d.TokamakSmallCartesian,
           ChargedParticleDynamics.GuidingCenter3d.TokamakSmallCylindrical,
           ChargedParticleDynamics.GuidingCenter3d.TokamakSmallToroidal,
           ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCartesian,
           ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCylindrical,
           ChargedParticleDynamics.GuidingCenter3d.TokamakIterCylindrical,
           ChargedParticleDynamics.GuidingCenter3d.SolovevSymmetricField,
           ChargedParticleDynamics.GuidingCenter3d.SolovevIter,
           ChargedParticleDynamics.GuidingCenter3d.SolovevIterXpoint]
```
