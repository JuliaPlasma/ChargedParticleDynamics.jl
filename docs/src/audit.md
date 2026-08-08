
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
``\{ g_{1}, g_{2} \} = \pm b_{m} [ B + (p - A) \cdot (\nabla \times b) ]`` in their denominator, with
``m`` the index of the constraint the pair leaves out: ``+b_{3}`` for ``(g^{1}, g^{2})``, ``+b_{1}``
for ``(g^{3}, g^{1})``, ``-b_{2}`` for ``(g^{2}, g^{3})``. Each pair is therefore singular on a
different surface. The sign is a consequence of how each pair is ordered rather than of anything
physical — see `constraint_pair` — but it is spelled out because it means ``\lambda_{\mathrm{o}}`` and
``b_{m}`` need not share a sign in the tabulated values.

The implementation once used ``(g^{3}, g^{1})`` alone, and **eight of the eleven equilibria have
``b_{1} = 0`` exactly at the initial condition the package ships**, so the multipliers were infinite
and `hodeproblem` could not be started there at all. This was not confined to the curvilinear cases:
both cartesian tokamaks fail too, because their initial conditions sit at ``y = z = 0`` where
``B_{x} = 0``. What decides it is where the initial condition lies relative to the field, not the
coordinate system — and if anything a cartesian chart is the harder case for this pair, since there the
toroidal direction is split across ``b_{1}`` and ``b_{2}`` instead of sitting in ``b_{3}`` alone.

**All three pairs are now implemented and selectable** through the `constraints` keyword, and each
equilibrium declares a `default_constraints` that is regular at its own initial conditions. Nine of
the eleven equilibria integrate as a result, and the two that do not — `SymmetricField` and
`ThetaPinchField` — ship Poincaré-invariant loops rather than point initial conditions and never had
integration tests. `SolovevSymmetricField` is the case that made all three necessary rather than two:
``b = e_{2}`` there, so ``b_{1}`` and ``b_{3}`` vanish together and only ``(g^{2}, g^{3})`` survives.

Two residual limitations are worth knowing about, neither of them removable by choosing a pair.

* **A pair that is regular at the initial condition can still be badly conditioned along the orbit.**
  `Dipole3d` is regular for all three, but its orbit takes ``b_{1}`` and ``b_{2}`` through zero, and
  the two pairs that divide by them lose the orbit entirely past ``t = 30``. Its default is chosen on
  the conditioning along the orbit for that reason. `TokamakMediumCartesian`'s is chosen on conditioning
  too, at the initial condition rather than along the orbit: ``b_{2} = 0.9923`` is the largest component
  there, giving ``\lambda_{\mathrm{o}} = -3.85`` where ``(g^{1}, g^{2})`` gives ``0.48``. The other nine
  are chosen on regularity alone. See `TODO.md`.
* **`hodeproblem_canonical` can need a smaller step than the module declares.** On
  `TokamakMediumCartesian` it cannot hold the thousand-step example at ``\Delta t = 0.1`` on any pair
  regular there, and takes a thousand steps at ``\Delta t = 0.01`` instead. On `SolovevSymmetricField`
  it cannot produce a usable orbit at any step tried. Both are properties of that formulation rather
  than of the constraint pair, and in neither case does it raise — it returns a trajectory that has left
  the device, so it has to be screened on magnitude.
* **The omitted constraint is conserved ``1/b_{m}`` times worse than the retained pair**, since
  ``b_{1} g^{2} - b_{2} g^{1} + b_{3} g^{3} = 0`` determines it from them pointwise. This is why
  `compute_constraints` reports all three components rather than the two the pair retains.

The compact form of the paper's Eq. (29), `hodeproblem_compact`, is free of the ``b_{m}``
denominator altogether and does not depend on the pair; see
[Guiding Center Dynamics in 3D](@ref) for the derivation, which had to be redone
coordinate-generally because Eq. (29) as printed is cartesian. See [Findings](@ref) for the
measurements.

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

### Default time spans and steps

Every equilibrium module declares `DEFAULT_TIMESTEP` and `DEFAULT_TIMESPAN`, and the two are chosen
so that the default problem runs **1000 steps** in every model. That makes the shipped example
directly comparable between models, coordinate systems and formulations, and it is what the test
suite integrates: with the exceptions listed below, every test runs at its module's own defaults
rather than at numbers that exist only in the test file.

A declared default is not, however, a guarantee of a resolved orbit, and two limits are worth
knowing before using one for anything quantitative.

**A default is chosen per model, not per formulation.** Where a model offers several formulations,
they do not all tolerate the same step and span, and the default follows the majority rather than the
weakest. The known cases are `GuidingCenter3d.TokamakMediumCartesian`, where
`hodeproblem_canonical` loses `initial_conditions_deeply_passing` beyond ``t \approx 20`` while
`hodeproblem` and `hodeproblem_compact` hold the full span, and
`GuidingCenter3d.SolovevSymmetricField`, where `hodeproblem` and `hodeproblem_canonical` lose the
orbit beyond ``t \approx 10`` at *any* step size and only `hodeproblem_compact` survives the declared
example. The tests give those two formulations a shorter span explicitly.

**A step that is too coarse can lose an orbit silently.** The nonlinear solver frequently does not
warn while the relative energy error is of order one, so neither the solver's silence nor the
`isa GeometricSolution` assertion the tests make is evidence that a run is physical. The defaults
recorded here were each checked against energy conservation, and against the constraints for the 3D
guiding centre.

`GuidingCenter3d.Dipole3d` is deliberately not at its most accurate step: `Δt = 0.1` is chosen to
*exhibit* the constraint drift that distinguishes the three formulations rather than to minimise it.

### Nonlinear solver tolerances

Every test module, script and example that integrates an implicit method passes
`f_abstol = 1E-12`. `f_abstol` is an *absolute* bound on the residual, and the variational and Pauli
residuals carry the momentum equation ``p = \vartheta(q)``, so their round-off floor is
``\|\vartheta\| \, \varepsilon`` — about ``3.3 \times 10^{-15}`` for the ITER-scale Solov'ev
equilibria and ``5.7 \times 10^{-14}`` for `SolovevSymmetricField`. A tolerance below that floor is
unsatisfiable, and the solver spends its whole iteration budget failing to meet it.

The library default does not clear those floors. Since `GeometricIntegratorsBase` 0.5.1 it scales
with the size of the stage system,

```julia
f_abstol = max(8, solversize(method, problem)) * eps(datatype(problem))
```

which comes to ``1.8 \times 10^{-15}`` for the 4D guiding centre and Pauli solves
(`solversize = 8`) and ``2.7 \times 10^{-15}`` for the 3D `hodeproblem` under `PartitionedGauss(2)`
(`solversize = 12`). Both are below the ITER floor, so the explicit `f_abstol` is required rather
than merely conservative. Caller options are merged into that default, not substituted for it, so
passing one option does not disturb the others.

`f_reltol` is deliberately left at its `SimpleSolvers` default of ``\sqrt{\varepsilon}``: its
relative term ``f_{reltol} \|F(x_0)\|`` is what lets a large-magnitude solve converge at all.

See [Findings](@ref) for the measurements behind these numbers.

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

* Choose each 3D guiding centre equilibrium's `default_constraints` on the conditioning of the pair
  along its orbit rather than on regularity at the initial condition. Two of the eleven are currently
  chosen on conditioning: `Dipole3d`, along the orbit, where the difference is four orders of
  magnitude, and `TokamakMediumCartesian`, at the initial condition. The other nine are on regularity
  alone.
* Decide whether `GyroKinetics4d` should rescale by the physical ``B^{\star}_{\parallel}`` rather than
  by `ωabs`, which is ``\det(DF) \, B^{\star}_{\parallel}`` and therefore negative in the left-handed
  charts, so that the rescaled time ``s`` always runs with physical time. Nothing is broken as it
  stands; see `TODO.md`.
* Place the Pauli initial conditions on the *first-order* slow manifold rather than the lowest-order
  one. `initial_conditions(x, u, μ)` sets ``v = u \vec{b}``, which omits the drift velocity, and that
  omission — not the difference between the models — is what separates the Pauli orbit from the
  guiding centre one; see [Findings](@ref).
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
