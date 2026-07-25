# Changelog

All notable changes to ChargedParticleDynamics.jl are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [Unreleased]

### Added

- **Guiding centre dynamics in 3D.** A new `GuidingCenter3d` module implements the guiding
  centre equations as a constrained canonical Hamiltonian system in the position `r` and its
  conjugate momentum `p = A(r) + u b(r)`, with the parallel velocity recovered as a dependent
  variable. The two independent components of `v × b = 0` are imposed as constraints `g₁` and
  `g₂`, whose Lagrange multipliers `λ₁` and `λ₂` are determined algebraically. Constructed by
  `hode`. Available for thirteen equilibria: dipole, quadratic potentials, symmetric quadratic,
  theta pinch, symmetric Solov'ev, ITER-like Solov'ev with and without X-point, the small tokamak
  in cartesian, cylindrical and toroidal coordinates, the medium tokamak in cartesian and
  cylindrical coordinates, and the ITER-size tokamak in cylindrical coordinates.
- **Canonicalised 3D guiding centre formulation.** `hode_canonical` absorbs the multipliers into
  the Hamiltonian, `H̄ = H + λ₁ g₁ + λ₂ g₂`, which agrees with `H` on the constraint manifold but
  yields a genuinely canonical system. It requires the second derivatives of the magnetic field
  and lives in `src/guiding_center_3d/guiding_center_3d_canonical.jl`, included from all thirteen
  equilibria. `hamiltonian_canonical` is now exported alongside it.
- **Makie plotting.** A `ChargedParticlePlots` package extension, loaded when `Makie` and
  `LaTeXStrings` are available, provides `plot_fieldlines`, `plot_invariant`,
  `plot_invariant_error`, `plot_components`, `plot_trajectory_cartesian`,
  `plot_trajectory_cylindrical`, `plot_trajectory_poloidal`, `plot_trajectory_projection`,
  `plot_trajectory_3d`, `plot_poincare_invariant_error`, `plot_poincare_loop`,
  `plot_poincare_surface` and `plot_poincare_trajectories`, together with their mutating
  variants. These are plain functions returning `(figure, axis)` and taking plain coordinate
  vectors, rather than `Plots.jl` recipes taking solution objects.
- `cartesian_solution` converts a solution to cartesian coordinates, returning the time series
  and the `R`, `X`, `Y` and `Z` components.
- Documentation page for the 3D guiding centre model, and a worked example integrating the
  charged particle, Pauli particle and 3D and 4D guiding centre models in an ITER-like
  equilibrium in cylindrical coordinates.
- `initial_conditions_default` for the small tokamak in cartesian coordinates.
- Example script for the Pauli particle.
- **`default_parameters(::Type{T} = Float64)` for every equilibrium module** of every model
  family, following the `GeometricProblems` convention: it returns a named tuple and honours its
  type argument. For the guiding centre and Pauli models this is the magnetic moment `μ` — for the
  Pauli models computed from the module's default `(qᵢ, vᵢ)` by splitting the velocity at `qᵢ`.
  The charged particle models take no parameters, since the electromagnetic field is injected as
  code rather than passed in, and return an empty named tuple, so that every problem in the package
  can now be constructed the same way. Where a module already had a `const parameters`, it is now
  defined as `default_parameters()`; the seventeen modules that had none previously raised an
  `UndefVarError` when a constructor was called without an explicit parameter argument.
- **A documentation page for the Pauli particle** (`docs/src/pauli_particle_3d.md`), which was
  previously a heading and an `@autodocs` block. It derives the classical Pauli particle from
  Pauli's Hamiltonian, tabulates its Lagrangian against those of the charged particle and the
  guiding centre, explains the slow-manifold relation to guiding centre dynamics, and records the
  reference initial condition shared by the three models.
- **A model audit page** (`docs/src/audit.md`) recording, per model, the simplifications made (the
  `m = e = 1` normalisation, the missing electrostatic potential of the 4D guiding centre, the
  diagonal-metric restriction, the cartesian-only noncanonical charged particle), the known
  singularities (the `b₁ = 0` constraint pair of the 3D guiding centre), the interface conventions,
  and the state of the documentation.
- **`test/structure_tests.jl`**, a set of cheap tests over the whole cross product of model
  families and equilibria, which the integration tests cannot cover because each integration is
  expensive. It asserts that every `initial_conditions_*` and every problem constructor of every
  module runs, that every module has a type-respecting `default_parameters`, that the analytic
  second derivatives of the 3D guiding centre Hamiltonian and the 4D `ḡ` agree with central
  differences, that the projection forces are linear in the multiplier, and that the two-forms are
  antisymmetric. Every one of these corresponds to a bug listed under *Fixed* that no existing test
  would have caught.
- References to the literature the models follow, on the index page and in the audit.
- **A findings page** (`docs/src/findings.md`) and **`TODO.md`**, together with the five
  `scripts/study_*.jl` that reproduce the former. The findings page records the numerical studies
  made during the audit that characterise the models rather than merely verify a line of code: the
  conditioning of the 3D guiding centre constraint pair, which expression is the conserved toroidal
  momentum, the exact proportionality between the gyrokinetic and the guiding centre vector fields,
  the volume preservation of the gyrokinetic splitting measured against step size, and the
  conservation levels of each family side by side. `TODO.md` carries the open work.
- **Diagnostics for the 3D guiding centre model.** `guiding_center_3d_diagnostics.jl` was an empty
  file, included by all thirteen equilibria, so the model had no `compute_energy`,
  `compute_energy_error` or `compute_toroidal_momentum` at all. Alongside them,
  `compute_constraints` returns the two constraints ``g₁`` and ``g₂`` as data series: they vanish
  along the exact flow, so their magnitude is the drift off the manifold on which `p = ϑ`, which is
  the quantity the approximately symplectic methods are designed to keep bounded.
- **A documentation page for `GyroKinetics4d`** (`docs/src/gyro_kinetics_4d.md`), which had none
  and was not listed in the index. It derives the gyrokinetic characteristics, the rescaled time
  and the ``\beta``/``\gamma`` potentials from the notes the module follows, explains the six-way
  splitting and why it is volume preserving, and records what is and is not implemented.

### Changed

- Updated to the current `GeometricEquations` problem interface and the current
  `GeometricSolutions` diagnostics interface. `compute_energy`, `compute_toroidal_momentum` and
  their error variants now take the problem parameters, following the new `compute_invariant`
  and `compute_invariant_error` signatures; the convenience methods taking a `GeometricSolution`
  obtain them from the problem. `compute_momentum_error`, `compute_one_form` and
  `compute_error_drift` are re-exported from `GeometricSolutions`.
- **Updated to GeometricEquations 0.21, GeometricIntegrators 0.16, GeometricSolutions 0.6 and
  PoincareInvariants 0.5**, the versions the companion packages and `GeometricExamples` are built
  on.
- **The 4D guiding centre Poincaré integral invariants were rewritten for the
  `PoincareInvariants` 0.5 interface.** The previous `guiding_center_4d_loop.jl` and
  `guiding_center_4d_surface.jl` were still written against the 0.3-era API
  (`PoincareInvariant1st(init, f, ϑ, Δt, D, n, ntime, nsave, DT)` and its `2nd`,
  `2ndTrapezoidal` siblings), which had already been removed in 0.4, and they built the ensemble of
  initial conditions by hand. They now provide the pieces the new interface composes:
  `guiding_center_4d_poincare_invariant_1st(N)` and `..._2nd(N)` return a noncanonical `FirstPI`
  or `SecondPI` built from the guiding centre one- and two-form; `guiding_center_4d_loop_ode`,
  `..._iode`, `..._lode` (and the `surface` counterparts) return the base problem at the loop's or
  surface's magnetic moment; and `guiding_center_4d_loop_ensemble(prob, pinv)` samples the
  parameterisation into a `GeometricEquations.EnsembleProblem` to integrate and pass to `compute!`.
  Verified on the medium tokamak in cylindrical coordinates: the invariant is independent of the
  sampling to eight digits and its drift falls as `Δt²`, i.e. it is limited by the integrator
  rather than by the quadrature. The `vode` variants are gone with
  `guiding_center_4d_vode_formal_lagrangian`, which no longer exists; use the `lode` ones.
- `ϑ` and `ω` of the 4D guiding centre gained `(out, t, q, params)` methods, so they can be handed
  to `PoincareInvariants` — and anything else following the `GeometricEquations` form convention —
  without a wrapper. Neither depends on the parameters.
- **The plotting extension is now covered by the test suite** (`test/plots_tests.jl`). It was
  previously never loaded by the tests, since they pull in neither `Makie` nor `LaTeXStrings`, which
  is how the two breaks listed under *Fixed* went unnoticed. The new testset exercises the
  diagnostics, the trajectory and field-line plots, and the full Poincaré path from problem through
  ensemble and integration to `compute!` and the plots.

  Note that `plot_invariant` now has to be qualified as `ChargedParticleDynamics.plot_invariant`
  whenever `PoincareInvariants` is loaded too, which the guiding centre Poincaré workflow requires:
  `PoincareInvariants` 0.5 exports a `plot_invariant` of its own, for the error of a Poincaré
  invariant over an `EnsembleSolution`. The two are unrelated functions sharing a name.
- Initial conditions and coordinate helpers moved from `src/initialization` to `src/utils`.
- Minimum Julia version raised from 1.6 to 1.10.
- CI now tests Julia 1.10, 1.11, 1.12 and nightly.
- The documentation is built from its own environment (`docs/Project.toml`) rather than the
  package environment.
- Various initial conditions, time steps and integration intervals of the 3D guiding centre
  equilibria were adjusted.

### Fixed

- **`compute_energy`, `compute_energy_error`, `compute_toroidal_momentum` and
  `compute_toroidal_momentum_error` no longer dispatched on a `GeometricSolution`.** Their
  low-level methods annotate the time axis as `::TimeSeries`, but since GeometricSolutions 0.6
  restructured the solution to hold a vector of states, `sol.t` is a `ScalarDataSeries` — so the
  convenience methods, which forward `sol.t`, raised a `MethodError` for every model. The
  annotation is now `Union{TimeSeries, ScalarDataSeries}`, the union
  `GeometricSolutions.compute_invariant` itself accepts.
  `src/guiding_center_4d/guiding_center_4d_diagnostics.jl`,
  `src/charged_particle_3d/charged_particle_3d_canonical.jl`,
  `src/charged_particle_3d/charged_particle_3d_noncanonical.jl`.
- **`plot_poincare_invariant_error` raised a `MethodError`** for the output of
  `PoincareInvariants.compute!`: it went through `GeometricSolutions.compute_relative_error`, which
  takes a `ScalarDataSeries`, while `compute!` returns a plain `Vector`. The relative error is now
  computed in place, as the `O(Δt²)`-corrected variant beside it always did.
  `ext/plots/poincare_invariants.jl`.
- **The test suite asserted `@test_nowarn` around every integration**, which is not a usable
  criterion for these models: the stiffer guiding centre and Pauli particle orbits routinely
  exhaust the Newton solver's iteration budget at the requested accuracy, and `SimpleSolvers` warns
  once per occurrence — over thirty thousand times in a single testset. Which orbits trip it depends
  on the floating-point details of the platform, so the suite failed on the 4D symmetric Solov'ev
  equilibrium on Linux and on the Pauli particle in cartesian coordinates on Windows, and on a
  third equilibrium again after the dependency update. Relaxing the tolerance does not help: the
  warning count stays in the same order of magnitude at `1E-13` and rises for some problems. The
  tests now assert that `integrate` returns a `GeometricSolution`, and the warnings are filtered
  and counted by the new `test/quiet_solver_warnings.jl`.
- `compute_energy` and friends threw a `MethodError` for every model, because the imports had
  been switched from `GeometricProblems.Diagnostics` to `GeometricSolutions` without updating
  the call sites for the added `parameters` argument.
- `compute_momentum_error` in the 4D guiding centre diagnostics still called into
  `GeometricProblems`, and `compute_one_form` and `compute_error_drift` were exported without
  being defined.
- `convert_coordinates_RZphi_to_xyz` used the removed `SDataSeries`; it is replaced by
  `cartesian_solution`.
- `toroidal_momentum` referred to the undefined `ϑ3` instead of `ϑ₃` in the symmetric Solov'ev
  equilibrium and in the small tokamak in cartesian, cylindrical and toroidal coordinates, in
  both the 3D and the 4D guiding centre models.
- The medium tokamak in cartesian coordinates exported `toroidal_momentum` without defining it.
- `plot_trajectory_3d_figure` declared positional instead of keyword arguments, so any keyword
  passed to `plot_trajectory_3d` raised a `MethodError`.
- `plot_invariant_error` divided by `log10(0)` for an identically vanishing error series, and
  had an unbalanced bracket in its axis label.
- The test suite aborted before the 4D guiding centre, gyrokinetic and Pauli particle tests,
  which imported `MidpointExtrapolation` from `GeometricIntegratorsBase` — no longer part of the
  test environment. Once reachable again, those tests revealed solvers that did not converge at
  the default time step of some equilibria, and a degenerate Hermite extrapolation of the initial
  guess for the Pauli theta pinch, whose momentum is constant along the initial trajectory.
- The example block plotting the ITER magnetic field was missing its `@` and was never executed;
  the initial conditions figure of the initialization page is reactivated using Makie.
- Numerous fixes to the Pauli particle Solov'ev equilibria, the 3D guiding centre equations and
  the scripts and examples.
- **The mixed second derivatives of the 3D guiding centre Hamiltonian were wrong in every
  curvilinear coordinate system.** The metric term of `∂²H/∂qᵢ∂qⱼ` is
  `vₖ (∂gᵏᵏ/∂qᵢ ∂vₖ/∂qⱼ + ∂gᵏᵏ/∂qⱼ ∂vₖ/∂qᵢ)`. For `i = j` the two summands coincide and collapse to
  a factor of two, and the factor of two had been written for `i ≠ j` as well. The Hessian was
  therefore not symmetric — `d²Hdq₁dq₂` and `d²Hdq₂dq₁` differed by 65 % on the ITER Solov'ev
  X-point equilibrium — and `hode_canonical`, which differentiates the Lagrange multipliers and so
  depends on it, was wrong for every equilibrium except the cartesian ones. The derivatives are now
  checked against central differences in cartesian, cylindrical and toroidal coordinates.
- **`d²Hdq₃dq₁` and `d²Hdq₃dq₂` raised an `UndefVarError`**, referring to `d²v₁dq₃dq₁` and
  siblings that are defined nowhere. They were unreachable only because the multiplier derivatives
  happen to use the transposed variants. They have been removed along with thirteen other unused
  transposed duplicates, which is what let the two halves of the Hessian drift apart in the first
  place; the file is 150 lines shorter.
- **Two of the four `∂²ϑ/∂x_k∂x_j` blocks of the 4D guiding centre were wrong**, and none of them
  was reachable. `D²ϑd₂` wrote its first block to row 2 instead of row 1, so that block was dead
  code and row 1 was never assigned, and used `A₂` where `A₃` was required; `D²ϑd₃` had the same
  row error and used `d²A₃dx₂` where `d²A₂dx₃` was required. Separately, the `ḡ₁`–`ḡ₄` helpers
  called them with the output matrix in the last argument slot rather than the first, which raised
  a `MethodError`. `ḡ` is now checked against central differences of the one-form gradient.
- **`initial_conditions_*` raised a `MethodError` in six equilibrium modules** of the 3D guiding
  centre — `merge(nt, params = ...)` passes `params` as a keyword argument, and no such `merge`
  method exists; it must be `merge(nt, (params = ...,))`, as the other modules already did. This
  affected `SolovevIter`, `SolovevSymmetricField`, `TokamakIterCylindrical`,
  `TokamakSmallCylindrical` and `TokamakSmallToroidal`, i.e. every one of those equilibria was
  unusable.
- **The projection force `g = λ ⋅ ∇ϑ` was not linear in `λ`.** The Pauli particle version was built
  from `v` and never referred to `λ` at all, so the projection did not depend on what it was
  projecting; the canonical charged particle version used `λⱼ²` in the metric term. Both now use
  `(∂Aⱼ/∂qⁱ + ∂gⱼⱼ/∂qⁱ vʲ) λʲ`, and linearity is asserted in the tests.
- **The Lagrangian of the canonical charged particle omitted the `A·v` term**, so it was not the
  Legendre dual of its own one-form `ϑ = g v + A` and `charged_particle_3d_lode` described a
  different system than `charged_particle_3d_pode`.
- **`charged_particle_3d_lode` raised an `UndefVarError` in all eight canonical modules**: it
  passes `ω`, which was defined only in the noncanonical file. The canonical formulation now has
  its own two-form — the Lagrangian is regular, so it is the canonical form on the six-dimensional
  phase space.
- **`ω` and `dϑ` of the noncanonical charged particle referred to `dϑ₁dx₂` and siblings that are
  defined nowhere**, and both had their output argument in the wrong slot for the
  `GeometricEquations` convention, so the `LODEProblem` built from `ω` could never be evaluated.
  The derivatives are now defined, including the metric terms, and both functions take the output
  first. Their velocity block is `∂ϑᵢ/∂vʲ = gᵢⱼ`; it had been hard-coded to the identity in `ω` and
  to zero in `dϑ`.
- **`charged_particle_3d_sode` could never have run**: its velocity map referred to undefined `q`,
  `t` and `Δt`, called `eye`, removed from Julia in 0.7, and built a rotation matrix whose signs
  were not those of `v × B`; the constructor passed an undefined `x₀`. It is now the Cayley
  transform of `v̇ = v × B`, which is the rotation the implicit midpoint rule produces and
  preserves `|v|` exactly.
- **The `κ` wrappers of `guiding_center_4d_dg` all put the output array in the wrong slot** and so
  raised a `MethodError`, and the problem was built with `qᵢ` as its own initial momentum instead
  of the `κ`-modified one-form evaluated there.
- **`guiding_center_4d_iode_λ` called `guiding_center_4d_λ₀`, which does not exist**, and the
  underlying `guiding_center_4d_λ` passed a bare `μ` to a function that unpacks it from a parameter
  named tuple, while `guiding_center_4d_λᵢ` called a seven-argument method that was never defined.
  None of this had ever run.
- **`guiding_center_4d_loop_ode` and its siblings passed `parameters` as a keyword argument**
  although it is the second positional argument of the underlying constructors, referred to
  `guiding_center_4d_vode_formal_lagrangian`, which does not exist, and `..._loop_vode` initialised
  itself from the `iode` variant. (Superseded by the `PoincareInvariants` 0.5 rewrite above.)
- **`charged_particle_3d_periodicity` returned a bare vector of periods.** `GeometricEquations`
  recognises periodicity only as an `(xmin, xmax)` tuple of arrays, so it was accepted by the
  constructor and then silently ignored, and no coordinate was ever wrapped. It also hard-coded the
  third coordinate as an angle, which would have wrapped `z` in the cartesian equilibria once the
  periodicity started taking effect; the range now comes from the coordinate system itself, via the
  `rangemin`/`rangemax` injected with the field code, as it already did for the guiding centre
  models.
- **`charged_particle(ics; noncanonical = true)` returned a four-element vector**, concatenating
  the position with the scalar `|v|` rather than with the velocity vector.
- **`fix_initial_momentum` dropped the inverse metric** when recovering the parallel velocity,
  contradicting the `u` defined four lines above it in every curvilinear coordinate system. It now
  calls `u` rather than spelling the contraction out a second time.
- The two-form convention in the 4D guiding centre documentation was stated with the opposite sign
  to the one the code implements and `GeometricProblems` documents.
- The parallel and perpendicular velocities in the initialization documentation were missing the
  square roots implied by the energy split two equations above them, which the code has always
  had. The non-standard pitch-angle convention is now flagged there.
- The reference for the 3D guiding centre model was cited under the wrong title; both it and its
  companion paper are now cited in full, with DOIs.
- The docstring of the dipole equilibrium was a copy of the one for the quadratic potentials.
- **`toroidal_momentum` was not conserved, because of a spurious factor of `R`.** Fifteen of the
  sixteen 3D and 4D guiding centre equilibria defined it as `R(t,q) * ϑ₃(t,q)`, but `ϑ₃` is already
  the covariant φ-component and is the canonical momentum on its own: over 10³ time units on the
  small tokamak the relative variation is 2E-13 for `ϑ₃` and 3E-3 for `R ϑ₃`. In cartesian
  coordinates the third coordinate is `z` rather than an angle, so neither expression is right and
  the momentum is the generator of rotation about the z-axis, `x ϑ₂ - y ϑ₁`. In the 3D model, whose
  state is `(q, p)` with `p = ϑ`, the definition also indexed `q[4]` of a three-component vector
  and raised a `BoundsError` on every call. All sixteen are corrected and the conservation is now
  asserted in the tests.
- **`GyroKinetics4d` integrated the wrong system.** Its problem constructors substituted a
  coordinate transformation `q̃ = ω₀ q` into the vector field without applying the corresponding
  Jacobian, and evaluated the `hamiltonian` invariant on the transformed state. The transformation
  itself was a literal reading of the `R̃ = B*∥ R` of the notes, which denotes the time
  reparametrisation `dt = B*∥ ds` that the vector field already carries rather than a change of
  variables — and freezing `B*∥` at the initial condition, as `params.ω₀` does, would not
  reproduce that substitution in any case. The constructors now work in the untransformed state.
  The transformation functions are retained as utilities, documented, and no longer applied.

  With that removed, the module is verified: its vector field agrees with the 4D guiding centre
  model to machine precision after multiplication by `B*∥`, the six subsystems sum exactly to the
  full field, and the Strang composition holds the determinant of the one-step map at 1 for every
  step size tested while a first-order method's volume error grows linearly with the step. Its
  integration tests, all of which were commented out, are re-enabled.
- The default time step of `GyroKinetics4d` was that of the physical time, a factor of
  `B*∥ ≈ 2E2` too large for the rescaled time the module actually integrates in.
- The callbacks of `GyroKinetics4d` — `dH`, `ϑ`, both methods of `ω`, `β` and `γ` — all took their
  output array in the third argument slot rather than the first, so none of them conformed to the
  `GeometricEquations` convention.
- `guiding_center_4d_periodicity` of `GyroKinetics4d` returned a bare vector of periods, silently
  ignored, and computed it by calling an ambiguous `periodicity`; it now derives the range from the
  coordinate system like the other models.

### Removed

- The `Plots`, `RecipesBase`, `GeometricProblems`, `GeometricIntegrators`,
  `GeometricIntegratorsBase` and `ForwardDiff` dependencies, together with the `Plots.jl`
  recipes and the PyPlot-based Poincaré invariant plots they supported. `GeometricIntegrators`
  remains a test dependency.
- The unused `Documenter`, `EulerLagrange` and `Printf` dependencies, and the `SimpleSolvers`
  test dependency.
- **`guiding_center_4d_iode_dec128` and the `DecFP` dependency of the examples.** Support for
  number types other than `Float64` will be handled uniformly rather than by a problem constructor
  per precision. Note that `examples/guiding_center_4d_poincare_invariant_2nd_quad.jl` still uses
  `Dec128`; it is an untracked 2017 script written against a `GeometricIntegrators` API
  (`IntegratorVPRKpSymplectic`, `getTableauVPGLRK`, `set_config`) that no longer exists, and was
  left in place rather than deleted.
- Fifteen unused transposed second derivatives of the 3D guiding centre Hamiltonian, two of which
  raised an `UndefVarError`, and the four `D²ϑ₁`–`D²ϑ₄` Hessians of the 4D guiding centre, whose
  only consumers were the `PoincareInvariants` trapezoidal variants that the 0.5 rewrite removed.
- The dead Travis badge from the documentation index.


## [0.1.0] - 2024-12-13

Initial release.


[Unreleased]: https://github.com/JuliaPlasma/ChargedParticleDynamics.jl/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/JuliaPlasma/ChargedParticleDynamics.jl/releases/tag/v0.1.0
