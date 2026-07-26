
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

### The gyrokinetic Hamiltonian is the equilibrium one

Every model now carries the electrostatic potential: ``H = \tfrac{1}{2} u^{2} + \mu \vert B \vert
+ \varphi``, with ``-E`` in ``\nabla H``. `GuidingCenter4d` and `GyroKinetics4d` used to omit it.

What remains is one level up. `GyroKinetics4d` implements ``H_{0}``, while the notes it follows
give the zero Larmor radius Hamiltonian ``H^{\mathrm{zlr}}``, which carries ``A_{\parallel}`` and
``\langle \tilde{A}_{\parallel} \rangle^{2}`` terms on top of it and changes ``\beta`` and
``\gamma`` as well. See `TODO.md`.

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

### The noncanonical charged particle's splitting is cartesian only

The formulation itself is now fully curvilinear: `charged_particle_3d_v` is derived from
``\Omega \dot{z} = - \nabla H`` using the same metric-carrying ``\omega`` as the one-form,
Hamiltonian and Jacobian beside it, rather than being the plain cartesian Lorentz force it used to
be. `dH` and `charged_particle_3d_iode_f` carry the metric-derivative terms to match. All four
equilibria — including the toroidal `TokamakSmallNoncanonical` — are self-consistent.

!!! warning "`sodeproblem` requires a trivial metric"
    The Boris splitting does not generalise. At frozen position the curvilinear kick is *quadratic*
    in ``v``, so the Cayley transform is no longer its exact flow and the splitting has no exact
    solution map. `sodeproblem` therefore throws an `ArgumentError` unless the metric is trivial,
    rather than silently integrating the cartesian system in a chart where that is a different
    model. In practice this excludes `TokamakSmallNoncanonical`; use `odeproblem` or `iodeproblem`
    there.

### `plot_fieldlines` is restricted to axisymmetric cylindrical equilibria

It contours the poloidal flux ``\psi = R A_{\varphi}`` as `equ.A₃(0, x, y, 0) / x`, which reads the
first two coordinates as ``(R, Z)`` and the third as the toroidal angle. That is meaningless for a
cartesian or toroidal equilibrium, and it used to draw a plausible-looking wrong picture for them
rather than failing. It now rejects any chart it is not valid in, and so does
`plot_trajectory_poloidal(R, Z, equ)`, which forwards its equilibrium to it. The predicate is
`is_axisymmetric_cylindrical`, decided from the coordinate ranges the field code injects.

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
ships**, so the multipliers are infinite and `hodeproblem` cannot be started there at all — the Newton
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
on the phasespace ``z``, which is the convention `GeometricProblems` documents for
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
| [Guiding Center Dynamics in 3D](@ref) | Complete: constraints, multipliers, both `hodeproblem` and `hodeproblem_canonical`. |
| [Guiding Center Dynamics in 4D](@ref) | Complete. |
| [Pauli Particles in 3D](@ref) | Complete. |
| [Gyrokinetic Guiding Centre Dynamics in 4D](@ref) | Complete. Note that its independent variable is the rescaled time, not the physical time. |
| Equilibrium modules | Each states its coordinate system and equilibrium parameters. The provenance of the hard-coded ``\mu`` and ``u`` of the `initial_conditions_*` is not recorded anywhere in the repository, and the docstrings say so rather than inventing one. |


## Interface

Problem constructors follow the `GeometricProblems` scheme — `odeproblem`, `iodeproblem`,
`lodeproblem`, `podeproblem`, `sodeproblem`, `hodeproblem`, with a suffix only where two
formulations coexist in one module (`hodeproblem` vs `hodeproblem_canonical`). They take

```julia
odeproblem(q₀; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP,
               parameters = default_parameters())
```

and additionally accept the named tuple that every `initial_conditions_*` returns —
`(q = …, params = …)`, plus `p` or `v` where the model has one — so that
`odeproblem(initial_conditions_deeply_passing())` carries that condition's own ``\mu``.


## Open work

* Make the 3D guiding centre constraint pair selectable. This is the sharpest remaining limitation:
  seven of the eleven equilibria have ``b_{1} = 0`` at their shipped initial condition and cannot be
  started at all.
* Write the coordinate-transforming `IRK` integrator that `coordinate_transformations.jl` is waiting
  for. The 0.x-era sketch is on file at
  `src/gyro_kinetics_4d/irk_with_coordinate_transformation.jl` as the starting point; it is not
  included by any module and does not compile against the current `GeometricIntegrators`.
* Give `GyroKinetics4d` the full zero Larmor radius Hamiltonian, not just ``H_{0}``.
* Route the equilibrium modules through one initial-conditions layer. Each carries its ``u`` and
  ``\mu`` as hard-coded literals, although `src/utils/initial_conditions.jl` already provides
  `InitialConditions` to derive them from a position, a pitch angle and an energy, together with
  per-model converters — and no equilibrium module calls any of it. The literals are not
  reproducible from those functions either, and their provenance is not recorded; see
  [Initialization](@ref).
* `lodeproblem_formal_lagrangian` is a verbatim copy of `lodeproblem` rather than the formal
  Lagrangian of the equations of motion in noncanonical Hamiltonian form.

See `TODO.md` for the details of each.
