
# Pauli Particles in 3D

Pauli's Hamiltonian for an electron of charge ``-1`` is
```math
H_{\text{Pauli}} = \tfrac{1}{2} \, (p - A(x))^{2} - \dfrac{\hbar}{2} \, \sigma \cdot B(x) + \phi (x) ,
```
where ``\sigma = (\sigma_{x}, \sigma_{y}, \sigma_{z})`` is the vector of Pauli matrices describing
the intrinsic magnetic moment observed in the Stern-Gerlach experiment.
In most regimes of classical physics that intrinsic moment is negligible.
Surprisingly, however, introducing a *formal* magnetic moment term ``\mu \vert B \vert`` into the
Hamiltonian of a classical particle opens the way to structure-preserving integrators for guiding
centre dynamics.
This is the *classical Pauli particle*,
```math
H (x,p) = \tfrac{1}{2} \, (p - A(x))^{2} + \mu \vert B(x) \vert + \phi (x) ,
```
with corresponding Lagrangian
```math
L (x, \dot{x}) = \tfrac{1}{2} \vert \dot{x} \vert^{2} + A(x) \cdot \dot{x} - \mu \vert B(x) \vert - \phi (x) .
```

The construction follows Jianyuan Xiao and Hong Qin, *Slow manifolds of classical Pauli particle
enable structure-preserving geometric algorithms for guiding center dynamics*,
[arXiv:2006.03818](https://arxiv.org/abs/2006.03818) (2020).


## Relation to the Other Models

The Pauli particle sits between the [charged particle](charged_particle_3d.md) and the
[guiding centre](guiding_center_4d.md):

| Model | Lagrangian |
|---|---|
| charged particle | ``L = \tfrac{1}{2} \vert \dot{x} \vert^{2} + A(x) \cdot \dot{x} - \phi (x)`` |
| classical Pauli particle | ``L = \tfrac{1}{2} \vert \dot{x} \vert^{2} + A(x) \cdot \dot{x} - \mu \vert B(x) \vert - \phi (x)`` |
| guiding centre | ``L = ( A(x) + u \, b(x) ) \cdot \dot{x} - \tfrac{1}{2} u^{2} - \mu \vert B(x) \vert - \phi (x)``, with ``u = \dot{x} \cdot b(x)`` |

It differs from the charged particle only by the ``\mu \vert B \vert`` term, so it retains the full
six-dimensional phasespace and a *regular* Lagrangian — unlike the guiding centre model, whose
Lagrangian is degenerate.
That is the point: the Pauli particle can be integrated with ordinary symplectic or variational
partitioned Runge-Kutta methods, with no projection onto a constraint manifold.


## Slow Manifolds

Solutions of the classical Pauli particle with (almost) no gyration,
``\dot{x} \times b(x) \approx 0``, can be viewed as slow manifolds of the Pauli dynamics.
If
```math
\mu' = \dfrac{\vert \dot{x} \times b(x) \vert^{2}}{2 \vert B(x) \vert} \approx 0 ,
```
the gyro radius is close to zero and the particle position nearly coincides with its guiding centre.
From this perspective, guiding centre dynamics *is* the slow manifold of the classical Pauli
particle, and a structure-preserving integrator for the latter inherits the good long-time
behaviour on the former.

See also Joshua W. Burby, *Guiding center dynamics as motion on a formal slow manifold in loop
space*, Journal of Mathematical Physics **61**, 012703 (2020),
[doi:10.1063/1.5119801](https://doi.org/10.1063/1.5119801).


## Implementation

The models are written in the same curvilinear form as the other families, with the one-form
```math
\vartheta_{i} = g_{ij} (x) \, v^{j} + A_{i} (x)
```
and the Hamiltonian
```math
H = \tfrac{1}{2} \, g_{ij} (x) \, v^{i} v^{j} + \mu \vert B(x) \vert + \varphi (x) ,
```
where ``g_{ij}`` is the (diagonal) metric of the coordinate system and ``\mu`` is carried as a
problem parameter.
Mass and charge are normalised to one; see [Normalization](@ref).

`initial_conditions(x₀, v₀)` splits a velocity vector into the part parallel to the magnetic field
and the magnetic moment of the perpendicular part,
```math
u_{0} = v_{0} \cdot b (x_{0}) , \qquad
v_{\parallel} = u_{0} \, b (x_{0}) , \qquad
v_{\perp} = v_{0} - v_{\parallel} , \qquad
\mu = \dfrac{\vert v_{\perp} \vert^{2}}{2 \vert B (x_{0}) \vert} ,
```
and returns `(q = x₀, v = v∥, params = (μ = μ,))`, which every problem constructor of this module
takes directly — `podeproblem(initial_conditions(x₀, v₀))`.
See [Initialization](@ref) for the route from an energy and a pitch angle instead.

The problem is available in partitioned (`podeproblem`), Hamiltonian
(`hodeproblem`) and implicit (`iodeproblem`) form.


## Usage

```julia
using GeometricIntegrators
using ChargedParticleDynamics.PauliParticle3d.TokamakSmallCylindrical

problem  = hodeproblem(initial_conditions_pauli())
solution = integrate(problem, PartitionedGauss(1))
```

For the small axisymmetric tokamak in cylindrical coordinates, with
```math
A (R, Z, \varphi)
= \frac{B_{0} R_{0}}{2} \, \frac{Z}{R} \, \mathrm{d} R
- \frac{B_{0} R_{0}}{2} \, \ln \bigg( \frac{R}{R_{0}} \bigg) \, \mathrm{d} Z
- \frac{B_{0} r^{2}}{2 q_{0}} \, \mathrm{d} \varphi ,
\qquad B_{0} = 1 , \quad R_{0} = 1 , \quad q_{0} = 2 ,
```
the reference initial condition of the charged particle,
```math
x_{0} = [ 1.05, \, 0, \, 0 ] ,
\qquad
v_{0} = [ 2.1 \times 10^{-3}, \, 0, \, 4.3 \times 10^{-4} ] ,
```
corresponds to
```math
v_{\parallel} \approx [ 0, \, 1.074 \times 10^{-5}, \, 4.297 \times 10^{-4} ] ,
\qquad
\mu \approx 2.315 \times 10^{-6} ,
```
which is the magnetic moment `default_parameters()` returns for that equilibrium, and to the
guiding centre initial condition ``q_{0} = [x_{0}, u_{0}]`` with
``u_{0} = v_{0} \cdot b (x_{0}) \approx 4.299 \times 10^{-4}``.
The three models can therefore be compared directly on the same physical orbit — and
`scripts/study_model_agreement.jl` measures how directly.

!!! note "Which components these are"
    The two velocities above are quoted in **physical** components, ``\sqrt{g_{ii}} \, v^{i}``, which is
    what makes their third entries read ``4.3 \times 10^{-4}`` rather than the contravariant
    ``4.095 \times 10^{-4}`` the module actually stores — the chart has ``g_{33} = R^{2}`` and
    ``R = 1.05``. And ``u_{0}`` is the contraction ``v_{0} \cdot b`` with the *covariant* `b`, not
    ``\vert v_{\parallel} \vert``: those agree to three digits here only because the parallel velocity is
    almost entirely toroidal. The three triads `b`, `bₚ` and `b⃗` are covariant, physical and
    contravariant and are not interchangeable; see [Initialization](@ref).

    Both signs above were negative in this page before `ElectromagneticFields` 0.7.0, and one of the two
    changes was a correction rather than an update. ``v_{\parallel}`` and ``\mu`` are *invariant* under
    the field reversal — ``v_{\parallel} = (v \cdot b) \, \vec{b}`` flips twice — so their signs were
    simply wrong. ``u_{0}`` genuinely was negative in this chart, because ``b`` pointed the wrong way;
    it is positive now, and agrees with the cartesian chart of the same tokamak, which it did not
    before. See "The reversed field of the left-handed charts" in [Findings](@ref).


## Modules

```@docs
ChargedParticleDynamics.PauliParticle3d
```

Each equilibrium is its own module, and every one of them `include`s the same
`pauli_particle_3d.jl`. The model's functions are therefore documented once, below, under
`TokamakSmallCylindrical`, and hold verbatim for all eight. What differs between the modules is the
chart, the equilibrium parameters and the initial conditions, which is what these docstrings record:

```@docs
ChargedParticleDynamics.PauliParticle3d.SymmetricField
ChargedParticleDynamics.PauliParticle3d.ThetaPinchField
ChargedParticleDynamics.PauliParticle3d.TokamakSmallCartesian
ChargedParticleDynamics.PauliParticle3d.TokamakSmallToroidal
ChargedParticleDynamics.PauliParticle3d.TokamakIterCylindrical
ChargedParticleDynamics.PauliParticle3d.SolovevIter
ChargedParticleDynamics.PauliParticle3d.SolovevIterXpoint
```

```@autodocs
Modules = [ChargedParticleDynamics.PauliParticle3d.TokamakSmallCylindrical]
```
