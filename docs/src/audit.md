
# Model Audit

This page records what the models in this package do and do not represent: the simplifications
they make, the coordinate systems they are valid in, the places where they are singular, and the
state of their documentation. It is meant to be read before using a model for anything the test
suite does not already cover.

The models were checked against

* Xinjie Li, Ruili Zhang, Jian Liu, *Approximately symplectic Runge-Kutta methods for the guiding
  center dynamics*, Journal of Computational Physics **563**, 115069 (2026),
  [doi:10.1016/j.jcp.2026.115069](https://doi.org/10.1016/j.jcp.2026.115069);
* Ruili Zhang, Jian Liu, Tong Liu, Wenxiang Li, Xiaogang Wang, *Canonical Hamiltonian guiding
  center theory and classical intrinsic magnetic moment*, Frontiers of Physics **21**(2), 026200
  (2026);
* Jianyuan Xiao, Hong Qin, *Slow manifolds of classical Pauli particle enable structure-preserving
  geometric algorithms for guiding center dynamics*,
  [arXiv:2006.03818](https://arxiv.org/abs/2006.03818) (2020);

and the interfaces against
[GeometricEquations.jl](https://github.com/JuliaGNI/GeometricEquations.jl) and
[GeometricProblems.jl](https://github.com/JuliaGNI/GeometricProblems.jl).


## Simplifications

### Mass and charge are normalised to one

Every model sets ``m = e = 1``. This is not an omission but the normalisation derived in
[Normalization](@ref), and it is the same convention the reference papers use ("the mass and charge
of the original particle are normalized to unity"). The consequence is that `docs/src/charged_particle_3d.md`
carries ``m`` and ``e`` symbolically while the code does not: read those formulas with both set
to one.

The initialisation helpers in `src/utils/initial_conditions.jl` *do* carry physical mass and charge,
because they convert an energy in eV into a normalised velocity; see [Initialization](@ref).

### The strongly-magnetised normalisation is for reference only

[Normalization](@ref) derives a second, ``\epsilon``-ordered form of the Lagrangian for strongly
magnetised plasmas,
```math
L' = \bigg( v' + \frac{\hat{l}}{\hat{\rho}_{\mathrm{th}}} \, A' \bigg) \cdot \dot{x}'
     - \frac{\vert v' \vert^{2}}{2}
     - \frac{e \hat{\phi}}{m v_{\mathrm{th}}^{2}} \, \phi' ,
```
but no model exposes the two prefactors, and none is expected to. Every model uses
``\hat{\rho}_{\mathrm{th}} = \hat{l}``, ``e \hat{\phi} = m v_{\mathrm{th}}^{2}``; the drift-kinetic
and gyrokinetic orderings are recorded on that page for reference, not offered as a choice.

### The 4D guiding centre has no electrostatic potential

`GuidingCenter4d` uses ``H = \tfrac{1}{2} u^{2} + \mu \vert B \vert``, with no ``\varphi`` and no
``E`` in ``\nabla H``. `GuidingCenter3d` includes both, and both reference papers carry
``\Phi(X)``. The two guiding centre families therefore describe different systems whenever the
equilibrium has a non-zero electrostatic potential. This is the largest remaining physics gap.

### Metrics are diagonal

!!! warning "Orthogonal coordinate systems only"
    All models contract with ``g^{11}, g^{22}, g^{33}`` (or ``g_{11}, g_{22}, g_{33}``) only.
    **Off-diagonal metric components are silently ignored** — a model handed a coordinate system
    with a non-diagonal metric returns wrong answers with no indication that it has done so.

    This covers every coordinate system `ElectromagneticFields` currently provides — cartesian,
    cylindrical and toroidal are all orthogonal — and is a deliberate choice rather than an
    oversight: carrying the full metric would complicate the code substantially, since the
    one-forms, Hamiltonians, two-forms and every hand-written derivative would gain the
    off-diagonal sums. It will be generalised when a non-orthogonal coordinate system is actually
    needed.

### The 3D charged particle in noncanonical form mixes conventions

`charged_particle_3d_noncanonical.jl` builds its one-form, Hamiltonian, two-form and Jacobian from
the metric, but `charged_particle_3d_v` is the plain cartesian Lorentz force
``\dot{v} = E + v \times B``, with no metric or Christoffel terms, and `dH` differentiates a
Hamiltonian without the metric derivative terms that `hamiltonian` actually has. The
`charged_particle_3d_sode` splitting is cartesian for the same reason.

**Only cartesian configurations of this formulation are self-consistent.** Of the four equilibria
that use it, `SingularField`, `SymmetricField` and `ThetaPinchNoncanonical` are cartesian and
therefore fine; `TokamakSmallNoncanonical` is toroidal and is the one that is not.
The canonical formulation (`charged_particle_3d_canonical.jl`) is fully curvilinear and correct.

### The magnetic moment is a fixed parameter

``\mu`` is carried as a problem parameter in all guiding centre and Pauli models, which is correct
for these formulations. The companion paper of Zhang et al. instead derives the intrinsic magnetic
moment from the canonical structure; that treatment is not implemented here.


## Known singularities

### The 3D guiding centre constraint pair

The constrained canonical formulation replaces the parallel-velocity constraint by two of the three
independent components of ``v \times b = 0``. The JCP paper labels them

```math
\begin{aligned}
g^{1} &= b_{1} v_{3} - b_{3} v_{1} , &
g^{2} &= b_{2} v_{3} - b_{3} v_{2} , &
g^{3} &= b_{1} v_{2} - b_{2} v_{1} ,
\end{aligned}
```

and shows (§4.1–4.2) that the Lagrange multipliers carry
``\{ g^{1}, g^{2} \} = b_{3} [ B + (p - A) \cdot (\nabla \times b) ]`` in their denominator for the
pair ``(g^{1}, g^{2})``, and ``b_{1} [ \cdots ]`` for the pair ``(g^{1}, g^{3})``.

**This implementation uses the pair ``(g^{3}, g^{1})``**, so the multipliers are singular wherever
``b_{1} = 0``.

**Seven of the eleven equilibria have ``b_{1} = 0`` exactly at the initial condition the package
ships**, so the multipliers are infinite and `hode` cannot be started there at all — the Newton
solver hits a NaN on the first step. This is not confined to the curvilinear cases:
`TokamakSmallCartesian` fails too, because its initial condition sits at ``y = z = 0`` where
``B_{x} = 0``. What decides it is where the initial condition lies relative to the field, not the
coordinate system.

The four that do work — `Dipole3d`, `QuadraticPotentials3d`, `TokamakMediumCartesian` and
`SolovevIterXpoint` — are exactly the ones the test suite and the scripts exercise. The remaining
3D test blocks are commented out for this reason, and `scripts/guiding_center_3d.jl` displaces its
initial condition off the midplane by `sqrt(eps())` to work around it. Even `SolovevIterXpoint` is
only marginal, at ``b_{1} \approx 5.9 \times 10^{-3}``.

Choosing ``(g^{1}, g^{2})`` instead puts ``b_{3} = b_{\varphi}`` in the denominator — the largest
component of a tokamak field, and non-zero on the midplane — which would make most of these
equilibria usable. The choice is hard-coded, with the alternative commented out beside it in
`guiding_center_3d_equations.jl`. See [Findings](@ref) for the measurements and `TODO.md` for the
proposed fix.

### The two-form of the 4D guiding centre

``\det \Omega = (b \cdot \beta)^{2}`` with ``\beta = \nabla \times \vartheta``, so the vector field
is singular where the parallel component of ``\beta`` vanishes. This is intrinsic to the
formulation rather than to the implementation, and is far from the region the supplied initial
conditions explore.


## Interface conventions

### The symplectic two-form

All models use
```math
\Omega_{ij} = \frac{\partial \vartheta_{i}}{\partial z^{j}} - \frac{\partial \vartheta_{j}}{\partial z^{i}} ,
```
on the phase space ``z``, which is the convention `GeometricProblems` documents for
`MasslessChargedParticle`, `LotkaVolterra2d` and `PointVortices`. With it, the Euler-Lagrange
equations of ``L = \vartheta \cdot \dot{z} - H`` read ``\Omega \, \dot{z} = - \nabla H``.

Note that no `GeometricIntegrators` integrator currently evaluates the two-form of an `LODE` — the
only call site is commented out — so an error there would not show up in a simulation. The
`ω` of every model is checked for antisymmetry in `test/structure_tests.jl` instead.

### Callback signatures

Every callback follows `GeometricEquations`: the output array comes first, then `t`, then the
state, then any velocity/multiplier slots, then `params`.

```julia
ϑ(p, t, q, v, params)    f(f, t, q, v, params)    g(g, t, q, v, λ, params)
ω(Ω, t, q, v, params)    v̄(v, t, q, p, params)    f̄(f, t, q, v, params)
```

### Periodicity

`GeometricEquations` recognises periodicity only as an `(xmin, xmax)` tuple of arrays. A bare
vector of periods is accepted by the constructor and then silently ignored, which is worth knowing
when adding a model.

### Parameters

Every equilibrium module provides `default_parameters(::Type{T} = Float64)` returning a named
tuple, following `GeometricProblems`. For the charged particle models, which take no parameters —
the electromagnetic field is injected as code rather than passed in — it returns an empty named
tuple, so that every problem in the package can be constructed the same way.


## Documentation coverage

| Model | State |
|---|---|
| [Normalization](@ref) | Complete. Derives the normalisation used by every model. |
| [Initialization](@ref) | Complete. Note the non-standard pitch-angle convention flagged there. |
| [Charged Particles in 3D](@ref) | Complete for both formulations. Reads ``m``, ``e`` symbolically; they are one in the code. |
| [Guiding Center Dynamics in 3D](@ref) | Complete: constraints, multipliers, both `hode` and `hode_canonical`. |
| [Guiding Center Dynamics in 4D](@ref) | Complete, but the Hamiltonian shown omits ``\varphi``, matching the code (see above). |
| [Pauli Particles in 3D](@ref) | Complete. |
| [Gyrokinetic Guiding Centre Dynamics in 4D](@ref) | Complete. Note that its independent variable is the rescaled time, not the physical time. |
| Equilibrium modules | Mostly a one-line docstring. `symmetric_quadratic.jl` and `theta_pinch.jl` document their Poincaré loop and surface in full; the rest state neither the coordinate system, the equilibrium parameters, nor the provenance of the hard-coded ``\mu`` and ``u`` values of their `initial_conditions_*`. |


## Open work

* Add ``\varphi`` and ``E`` to the 4D guiding centre Hamiltonian and its gradient. The same gap
  applies to `GyroKinetics4d`, where the notes give the zero Larmor radius Hamiltonian with the
  ``\phi`` and ``A_{\parallel}`` terms but only ``H_{0}`` is implemented.
* Make `charged_particle_3d_noncanonical` consistently curvilinear, or restrict it to cartesian
  equilibria.
* Make the 3D guiding centre constraint pair selectable.
* Provide more than the one equilibrium for `GyroKinetics4d`, and write the coordinate-transforming
  `IRK` integrator that its `coordinate_transformations.jl` is waiting for.
* Rename the problem constructors to the `GeometricProblems` scheme (`odeproblem`, `iodeproblem`,
  `lodeproblem`, …) and the keyword arguments to `timespan`/`timestep`.
