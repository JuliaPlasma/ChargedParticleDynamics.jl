# Changelog

All notable changes to ChargedParticleDynamics.jl are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [Unreleased]

**This release renames every problem constructor and their keyword arguments.** There are no
deprecation shims; see *Changed* below for the mapping.

### Added

- **All three constraint pairs of the 3D guiding centre model, and a `constraints` keyword to
  select between them.** Only ``(g^3, g^1)`` was implemented, taken from an earlier version of Li,
  Zhang & Liu; its Lagrange multipliers divide by `b₁`, which vanishes on the midplane of an
  axisymmetric equilibrium — where seven of the eleven equilibria that ship an initial condition sit.
  `hodeproblem` and `hodeproblem_canonical` now take `constraints = :g31 | :g12 | :g23`, and every
  equilibrium declares a `default_constraints` that is regular at its own initial conditions.
  `SolovevSymmetricField` needed all three to exist: `b = e₂` there, so two of the pairs are singular
  at once.

  Five of the seven commented-out test blocks in `test/guiding_center_3d_tests.jl` are re-enabled as
  a result, and neither `scripts/guiding_center_3d.jl` nor `docs/src/examples/iter_cylindrical.md` has
  to displace its initial condition off the midplane by `sqrt(eps())` any longer. The 3D guiding
  centre test file goes from four integration blocks to ten, and correspondingly from about four
  minutes to about ten in CI — the cost of covering nine equilibria instead of four, most of it
  compilation of one right-hand side specialisation per (equilibrium, formulation, pair).

  The constraints and all of their first and second derivatives are now written once against a
  generic index pair in the new `src/guiding_center_3d/guiding_center_3d_constraints.jl`, rather than
  by hand per pair, which shrank `guiding_center_3d_canonical.jl` from 722 lines to 361 rather than
  tripling it.
- **`hodeproblem_compact`, the compact form of the 3D guiding centre equations** — Eq. (29) of Li,
  Zhang & Liu, which drops the constraint-proportional terms from the Hamilton-Dirac right-hand side.
  That equation as printed is cartesian and assumes `H = ½Σ(pᵢ-Aᵢ)²`, so it is derived
  coordinate-generally here instead of ported; the derivation is at the top of
  `guiding_center_3d_compact.jl` and written up in `docs/src/guiding_center_3d.md`, and is checked
  against the paper's expressions on the three cartesian equilibria in `test/structure_tests.jl`.

  In that form it is *independent of the constraint pair* and free of every `bₘ` denominator, which
  is what `constraints = :parallel`, its default, selects. Passing a pair symbol instead reproduces
  the paper literally, `(pₘ-Aₘ)/bₘ` and its singularity included.
- **Seven more equilibria for `GyroKinetics4d`**, which previously shipped only the ITER-like
  Solov'ev equilibrium with X-point: `SolovevIter`, `TokamakIterCylindrical`,
  `TokamakMedium{Cartesian,Cylindrical}`, `TokamakSmall{Cartesian,Cylindrical,Toroidal}`. Each
  mirrors the `GuidingCenter4d` module of the same equilibrium and reuses its initial conditions.
  Their default time steps are the 4D guiding centre's divided by `B*∥` at their own initial
  condition, because the module integrates in the rescaled time with `dt = B*∥ ds`; the factor is
  computed per equilibrium and asserted in the tests, since getting it wrong by `B*∥ ≈ 2E2` is the
  bug the audit found in the one shipped module.

  Three of the eleven `GuidingCenter4d` equilibria have no counterpart. `SymmetricField` and
  `ThetaPinchField` carry only a Poincaré loop and surface parameterisation and no point initial
  condition, and this model has no loop/surface machinery. `SolovevSymmetricField` cannot have one:
  `SolovevSymmetric.@code` injects the equilibrium parameters `α` and `β` into the module and `β`
  collides with the vector potential `β` of `gc_common.jl`.
- **The electrostatic potential in the 4D guiding centre and gyrokinetic models.** Both used
  `H = ½u² + μ|B|`; they now carry `+ φ` in the Hamiltonian and `- Eᵢ` in its gradient, matching
  `GuidingCenter3d`, the charged particle and the Pauli particle. `GyroKinetics4d` additionally
  carries `- dEⱼdxᵢ` in its second derivatives. No shipped equilibrium of either model has a
  non-zero `φ` — the only `ElectromagneticFields` configurations that do are `QuadraticPotentials`,
  the three Penning traps and `EzCosZPerturbation`, and of those only `QuadraticPotentials` is used
  here, by `GuidingCenter3d` — so `φ`, `E` and `dE` are all *exactly* zero and no result changed.
  It becomes a behaviour change the moment such an equilibrium is added.
- **Named-tuple constructors.** Every problem constructor accepts the named tuple that
  `initial_conditions_*` returns directly: `odeproblem(initial_conditions_deeply_passing())`, with
  that condition's own `μ` travelling in `params`.
- **Test blocks for `Dipole3d` and `QuadraticPotentials3d`** in the 3D guiding centre suite. Neither
  was blocked by the `b₁ = 0` singularity — `b₁` is -0.41 and 2E-3 at their initial conditions — and
  neither had a test block at all, so these took the 3D model's dynamics from two equilibria of eleven
  to four, before the selectable constraint pair took it to nine. `QuadraticPotentials3d` is the only
  shipped equilibrium with a non-zero potential, so it is also the only integration test that
  exercises the `φ` and `E` terms on a field where they do not vanish.
- **An integration test for the κ-dependent "dg" formulation** (`iodeproblem_dg`) at κ = 0, ½ and 1.
  Its `ḡ` was checked against finite differences but it had never been integrated, and κ = 0 — the
  default — reduces it to the plain `iodeproblem`, so the κ terms were entirely unexercised.
- **A test that the noncanonical charged particle's vector field is the one its own `ω` and `H`
  generate**, comparing it against `-Ω⁻¹∇H` and `dH` against central differences of `hamiltonian`.
- **A test sweeping every `GyroKinetics4d` equilibrium** for integration and volume preservation.
  The determinant of the one-step map stays at 1 to ~1E-9 in all eight, including the cylindrical
  and toroidal charts: the formulation is intrinsic — `ϑᵢ = Aᵢ + u bᵢ` in covariant components and
  `Ω = dϑ` need no metric, and the Liouville measure is `√det Ω = B*∥` in any chart — so the
  volume-preserving splitting carries over unchanged.
- **`src/gyro_kinetics_4d/irk_with_coordinate_transformation.jl`**, the 0.x-era sketch of the
  coordinate-transforming Runge-Kutta integrator, tracked as reference material for the port onto
  the current `IRK`. It had existed only in one working tree, which is why the 0.2.0 entry below
  records it as not being in the repository. It is not `include`d by any module and does not
  compile against the current `GeometricIntegrators`; a header comment says so, and `TODO.md`
  records what of it survives the port and what has to be rewritten. Renamed from
  `firk_with_coordinate_transformation.jl` to match the current integrator naming.
- Docstrings for the twelve charged particle and Pauli particle equilibrium modules that had none,
  and the coordinate system and equilibrium parameters throughout. The provenance of the hard-coded
  `μ` and `u` of the `initial_conditions_*` is not recorded anywhere in the repository; the
  docstrings say so rather than inventing one. `ThetaPinchNoncanonical` and
  `TokamakSmallNoncanonical` were missed in that sweep and are covered now, as are the five model
  modules `ChargedParticle3d`, `PauliParticle3d`, `GuidingCenter3d`, `GuidingCenter4d` and
  `GyroKinetics4d`, none of which had one. Every module in the package now carries a docstring.
- **A check that every module is documented and reaches the manual**, in `docs/make.jl`. It walks
  every submodule of the package and warns when one has no docstring of its own, or is named in no
  `@docs`/`@autodocs` block. Adding an equilibrium and forgetting the page is now a warning instead
  of silence — which is how the seven new `GyroKinetics4d` equilibria came to be documented in the
  source and absent from the manual without anything noticing.

  Documenter's own missing-docs check cannot serve this purpose here, and the reason is structural.
  Enabling it means passing `modules`, and `Documenter.DocumentBlueprint` keeps
  `submodules(modules; ignore = checkdocs_ignored_modules)` in one field used for two jobs: the
  bindings that *must* be documented, and the filter on what a `@docs`/`@autodocs` block is
  *allowed* to render. Because every equilibrium module carries its own copy of the shared model
  docstrings, passing all of them reports 382 docstrings as missing when each is excluded by design,
  and narrowing the set with `checkdocs_ignored_modules` makes the `@docs` blocks that name those
  modules fail with "no docs found" instead. The check above does the job directly.

### Changed

- **The 3D guiding centre right-hand side is 2.7× faster**, incidentally to the constraint rewrite:
  the old code evaluated `λ₁` and `λ₂` once per component of `v` and of `f`, so the two multipliers
  and the six Poisson-bracket terms behind them were computed three times per call. A hundred steps of
  `Dipole3d` with `PartitionedGauss(2)` went from 1383 ms to 517 ms, and `hodeproblem_canonical` from
  11.9 s to 10.1 s. Allocations and compile time are unchanged to the byte and the second.
- **`compute_constraints` returns three series where it returned two**, and `g₁`, `g₂` in the
  equilibrium modules now mean the paper's ``g^1``, ``g^2`` rather than the two members of the one
  hard-coded pair, which were ``g^3`` and ``g^1``. Anything that read `c.g₁`/`c.g₂` still reads a
  constraint that vanishes along the flow, but not the same one.

  In the same vein, `λₒ(t, q, p)`, `λ₁(t, q, p, params)` and `λ₂(t, q, p, params)` — the short forms
  the scripts in `scripts/` use — take the constraint pair as a trailing argument now and default it to
  the equilibrium's own `default_constraints()`. For the eight modules whose default is no longer
  ``(g^3, g^1)`` those calls therefore return the multipliers of a different pair than before. Pass
  `constraint_pair(:g31)` explicitly to get the old quantity.
- **Every problem constructor is renamed to the `GeometricProblems` scheme.** No deprecation shims.

  | Module | Was | Is |
  |---|---|---|
  | `GuidingCenter3d` | `hode`, `hode_canonical` | `hodeproblem`, `hodeproblem_canonical` |
  | `GuidingCenter4d` | `guiding_center_4d_{ode,iode,iode_λ,lode,dg,formal_lagrangian}` | `odeproblem`, `iodeproblem`, `iodeproblem_λ`, `lodeproblem`, `iodeproblem_dg`, `lodeproblem_formal_lagrangian` |
  | `GuidingCenter4d` | `guiding_center_4d_{loop,surface}_{ode,iode,lode}` | `{loop,surface}_{ode,iode,lode}problem` |
  | `GuidingCenter4d` | `guiding_center_4d_{loop,surface}_ensemble`, `..._poincare_invariant_{1st,2nd}` | `{loop,surface}_ensemble`, `poincare_invariant_{1st,2nd}` |
  | `GyroKinetics4d` | `guiding_center_4d_{ode,sode}` | `odeproblem`, `sodeproblem` |
  | `ChargedParticle3d` | `charged_particle_3d_{ode,pode,iode,lode,sode}` | `odeproblem`, `podeproblem`, `iodeproblem`, `lodeproblem`, `sodeproblem` |
  | `PauliParticle3d` | `pauli_particle_3d_{pode,hode,iode}` | `podeproblem`, `hodeproblem`, `iodeproblem` |

  Internal callbacks (`charged_particle_3d_v`, `guiding_center_4d_ϑ`, …) keep their prefixes.
- **`tspan`/`tstep` → `timespan`/`timestep`**, and the module constants `const tspan`/`const Δt` →
  `const DEFAULT_TIMESPAN`/`const DEFAULT_TIMESTEP` in all 35 equilibrium modules. This removes the
  `tspan=tspan` shadowing hazard, where the keyword default and a module constant shared a name —
  and which, in the two `hode(x₀, parameters; tspan=tspan, …)` forwarding methods, would have
  silently pinned the forwarded span to the module default had the rename not been done together.
- **`parameters` is a keyword argument** defaulting to `default_parameters()`, matching
  `GeometricProblems`. It was the second positional argument.
- **Every `initial_conditions_*` returns the same named tuple**, `(q = …, params = …)` plus `p` or
  `v` where the model has one. The 4D guiding centre, gyrokinetic, Pauli and charged particle
  families returned bare tuples of two or three elements. `test/structure_tests.jl` no longer
  branches on `ics isa NamedTuple`; it asserts the shape instead.
- **`charged_particle_3d_sode` now throws for a curvilinear equilibrium.** See *Fixed* below.
- The two 4D guiding centre testsets for `SymmetricField` and `ThetaPinchField` ran **zero tests**:
  their contents were commented out because the calls passed sampling counts positionally —
  `loop_odeproblem(nl)` — which the `PoincareInvariants` 0.5 rewrite replaced with keyword-only
  constructors. They are restored and now exercise the loop and surface problems.
- `test/quiet_solver_warnings.jl` no longer calls the suppressed-warning count a "tripwire".
  Nothing asserts on it and nothing should: which orbits exhaust the solver's iteration budget is
  platform-dependent by construction, which is why the warnings are filtered in the first place, so
  any threshold tight enough to catch a regression would also flake. Both comments now say that
  plainly.
- "phase space" and "phase-space" are normalised to **"phasespace"** throughout the documentation,
  source comments, tests and scripts — the one-word spelling the package uses, which until now
  appeared only in `docs/src/charged_particle_3d.md` and `docs/src/normalization.md`. Comments and
  prose only; no code changes.
- The `*_periodicity` helpers pass the origin to `rangemin`/`rangemax` rather than the `±Inf` they
  had just filled `xmin`/`xmax` with. Those two take an evaluation point only for uniformity with
  the other generated field functions — a coordinate's range is a property of the chart, and
  `minx¹`…`maxx³` are baked in as literals — so nothing changes, but the code no longer reads as if
  the range were being asked for at the point at infinity.
- **The model pages document the shared model code once instead of once per equilibrium.** Every
  equilibrium module `include`s the same model source, so each of the fifty-odd modules carries its
  own copy of every shared docstring — and each page listed all of its modules in one `@autodocs`
  block, rendering the same functions ten or thirteen times over. `gyro_kinetics_4d.md` crossed
  Documenter's 200 KiB hard `size_threshold` and failed the build outright once the seven new
  equilibria were added to it; `guiding_center_4d.md` was 4 KiB short of the same fate. Each page now
  renders the shared API once, from a `TokamakSmallCylindrical` anchor, and lists the remaining
  modules in a `@docs` block, which emits each module's own docstring and nothing else. The pages
  drop from 204/196/136/116 KiB to well under the warning threshold and no page is near the limit.
- The example page `docs/src/examples/iter_cylindrical.md` asked for `f_abstol = 1E-15` with a
  1000-iteration budget. That tolerance is below the round-off floor of the ITER-scale residual,
  `‖ϑ‖ eps`, so the solver had no reachable stopping criterion and the line search ran to its limit
  on every step — the same defect the scripts were fixed for. It now uses the tolerances the test
  suite does.
- **The 3D guiding centre right-hand side is about twenty times cheaper**, in two independent
  changes, neither of which alters a trajectory by a single bit.

  `MidpointExtrapolation(5)` is no longer requested. It is not the default for either method —
  `PartitionedGauss` and `VPRKGauss` both declare `default_iguess() = HermiteExtrapolation()` — and
  nothing recorded why the 3D tests asked for it; the one documented reason is the Pauli theta pinch,
  whose momentum is constant along the initial trajectory, and that call site keeps it. It was **70 %
  of wall clock**, for 3.0–3.5× across all thirteen equilibria and all three formulations, with
  Newton converging in one iteration either way.

  The right-hand side no longer re-enters the generated field code once per term. `db₁dx₁` on the
  ITER Solov'ev X-point is 1905 statements containing 108 separate evaluations of `log(x₁)` and costs
  549 ns, against a `g¹¹` that folds to a constant, because `ElectromagneticFields` 0.6.2 emitted SymEngine's expanded tree
  with no common subexpression elimination. One evaluation of the Hamilton-Dirac right-hand side
  called it more than seventy times: each of the twelve Poisson-bracket terms re-entered from the
  top, and `λ₁` and `λ₂` each recomputed the `λₒ` they share. `fieldvalues` evaluates each injected
  function once into a `FieldValues`, which is passed where the coordinate vector would go —
  everything downstream is already generic in that argument, so the accessors pick it up by dispatch.
  `SecondFieldValues` does the same for the ninety-nine second derivatives `hodeproblem_canonical`
  needs, kept separate so that the other two formulations do not pay for derivatives they never read.

  `guiding_center_3d_v` goes from 24.7 µs to 3.6 µs — 6.9×, matching the end-to-end per-step ratio
  exactly. At the `ElectromagneticFields` 0.6.2 these two changes were developed against, and holding
  the initial guess as the suite used to set it, `hodeproblem` on `SolovevIterXpoint` goes from
  5.87 ms/step to 0.230 and `hodeproblem_canonical` from 19.9 to 1.57. (The 0.6.3 requirement below
  takes those to 0.082 and 0.39; [Findings](docs/src/findings.md) measures all four combinations.) A
  new block in `test/structure_tests.jl` pins every cached accessor to its uncached counterpart, so a
  missing or misindexed field cannot pass as a merely different trajectory. `src/guiding_center_3d/`.
- **The `3D guiding centre diagnostics` block asks for `f_abstol = 1E-12`** like the rest of the
  suite, rather than restating `GeometricIntegratorsBase`'s default of `8eps() = 1.8E-15`. It did so
  in order to measure what an unconfigured `integrate` delivers, but that is below the round-off floor
  of the ITER-scale residual — `quiet_solver_warnings.jl` and `study_solver_tolerances.jl` both say
  so — and four steps of the Solov'ev X-point ran to the iteration cap chasing it, at the crossing of
  `b₁ = 0` where `λₒ ∝ b₁` collapses. Those four converge at no budget at all: 50, 100, 200 and 1000
  iterations all terminate at the cap and all four give a bit-identical final state. The library
  defaults are still exercised by the calls that genuinely pass no options. Energy error is identical
  at both tolerances, and `structure_tests.jl` now emits no suppressed solver warnings at all.
  `test/structure_tests.jl`, `test/quiet_solver_warnings.jl`.
- **`ElectromagneticFields` 0.6.3 is now required**, up from 0.6. It eliminates common subexpressions
  in the field functions it generates, which this package injects into 32 equilibrium modules and
  evaluates on every right-hand side. `db₁dx₁` for the ITER Solov'ev X-point goes from 549 ns to 57
  and `d²b₁dx₁dx₁` from 1733 to 80 — the second derivatives, which `hodeproblem_canonical` needs,
  gain the most. The compat bound is tightened rather than left at `0.6` because the timings in
  [Findings](docs/src/findings.md) are no longer reproducible below it.

  Nothing changes in value. The pass only *names* subexpressions, and all **27456** field values this
  package injects across those 32 modules are bit-identical between 0.6.2 and 0.6.3 — checked on this
  package's own `@code(2.0, 5.0, 2.0)` and `@code_iter_xpoint()` parameterisations, which are not the
  ones `ElectromagneticFields` tests — as is every trajectory integrated from them, and every Newton
  iteration count along the way. `Project.toml`.
- **The cost tables in [Findings](docs/src/findings.md) are re-measured across four configurations**,
  because two independent changes bear on them and the previous tables could not tell them apart. The
  package's own right-hand side rewrite and the `ElectromagneticFields` release are now measured in
  both combinations, in interleaved randomised rounds on one machine. For `hodeproblem` the two are
  separable and almost exactly multiplicative — 6.9× here, 2.8× upstream, 19.5× together — and for
  `hodeproblem_canonical` they reinforce each other, 12.7× here on 0.6.2 but 15.4× on 0.6.3.

  This corrects an attribution error. The relative-cost table committed earlier credited the whole
  change to the right-hand side rewrite, but had in fact been measured with the upstream pass already
  active: `scripts/Manifest.toml` carried a `path = "../../ElectromagneticFields"` stanza with no
  matching `[sources]` entry in `scripts/Project.toml`, so `julia --project=scripts` resolved the
  dependency to a working checkout. `Pkg.status` cannot detect this — for a path dependency it prints
  the version recorded in the manifest, not the version at the path. The manifest is removed, and
  `findings.md` now says to check `Base.locate_package` instead.

  `cost()` in `scripts/study_guiding_center_3d_conditioning.jl` was also not fit to measure this: it
  timed a single `@elapsed` sample against its own admitted ~20 % spread, held `warn_iterations = 1`
  *inside* the timed region so that a fixed per-step logging cost was charged to every solve — which
  compresses exactly the ratio under test — and swallowed every exception into an indistinguishable
  `diverged`, so a misconfigured environment read as a result about the physics. It now takes the
  minimum of several `@timed` samples, asserts per sample that no compilation entered the timing,
  counts iterations in a separate untimed run, and reports what it catches.

### Removed

- **`const parameters` in all 52 equilibrium modules.** It was the second positional argument's
  default; with `parameters` now a keyword defaulting to `default_parameters()`, nothing referenced
  it any more. Call `default_parameters()` instead. Leaving a module constant with the same name as
  a keyword argument is the hazard the `tspan=tspan` note above describes.

### Fixed

- **``\partial\lambda_{2}/\partial q`` and ``\partial\lambda_{2}/\partial p`` of the 3D guiding centre
  carried the wrong sign** on the term proportional to ``\partial \{ g_{1}, g_{2} \} / \partial\cdot``.
  Both multipliers are a bracket over ``\lambda_{\mathrm{o}} = \{ g_{1}, g_{2} \}``, so both derivatives
  carry ``-\lambda \, \partial\lambda_{\mathrm{o}} / \lambda_{\mathrm{o}}`` whatever the sign of the
  numerator; ``\lambda_{1}`` had it, ``\lambda_{2}`` had the ``+`` the quotient rule gives before the
  minus of ``-\{ g_{1}, H \}`` is folded back in. The error was off by exactly
  ``2 \lambda_{2} \, \partial\lambda_{\mathrm{o}} / \lambda_{\mathrm{o}}``.

  Only `hodeproblem_canonical` uses these, and only multiplied by a constraint, so it vanished on the
  constraint manifold — which is why the canonicalised form still agreed with `hodeproblem` to 1E-8 and
  why nothing caught it. What it cost was accuracy and solver work. `TokamakMediumCartesian` now
  conserves energy to 2.2E-9 over a hundred steps where it managed 7.7E-9 before, and `Dipole3d` — the
  one equilibrium whose canonicalised solve needs more than one Newton iteration per stage — drops from
  3.8 iterations per stage to 2.6, peak 50 to 9, and from 29.6× to 18.8× the cost of the Hamilton-Dirac
  form. The right-hand side itself costs the same either way; the whole difference is the solve. The
  cost table in [Findings](docs/src/findings.md) is updated accordingly. (Those two ratios are the
  measurement of the day; the right-hand side has since been made some twenty times cheaper and the
  same row now reads differently — see the four-way table in [Findings](docs/src/findings.md).)

  It does *not* account for the canonicalised form's other weaknesses — `SolovevSymmetricField` still
  cannot complete an orbit under it, and `Dipole3d`'s conservation is unchanged — so the comparisons
  between the three formulations stand as written.

  Now pinned by a finite-difference test on all four multiplier derivatives in three charts, at a point
  displaced off the constraint manifold, since on it either sign passes.
- **`sodeproblem` of the noncanonical charged particle accepted `parameters` and dropped it.**
  Neither `SODEProblem` branch forwarded the keyword, so it — and anything the named-tuple form put
  into `ics.params` — was silently discarded, unlike in `odeproblem` beside it and in the
  `GyroKinetics4d` `sodeproblem`. Inert while these modules are parameter-free; now asserted in
  `test/structure_tests.jl`.
- **`charged_particle_3d_noncanonical` integrated a system inconsistent with its own structure.**
  Its one-form, Hamiltonian, two-form and Jacobian all carried the metric while
  `charged_particle_3d_v` was the plain cartesian Lorentz force and `dH` was the gradient of a
  different, metric-free Hamiltonian. The vector field is now derived from `Ω ż = -∇H` using the
  module's own `ω`,

      ẋⁱ = vⁱ ,    v̇ⁱ = ( -∂H/∂xⁱ + (v × β)ⁱ ) / gᵢᵢ ,

  with `β = ∇ × ϑ` the generalised magnetic field — the shape of the Lorentz force with `B` replaced
  by `β`, `E` by `-∇ₓH` and an overall `1/gᵢᵢ`. `dH` gains the metric-derivative terms and the
  `gᵢᵢ` factor on its velocity components, and `charged_particle_3d_iode_f` likewise. This reduces
  *exactly* to the old expression when the metric is trivial, so `SingularField`, `SymmetricField`
  and `ThetaPinchNoncanonical` are bit-identical; only `TokamakSmallNoncanonical`, which is
  toroidal, changes — and it changes because it was wrong.
- **`dH` of the noncanonical charged particle referred to a `dφdxᵢ` that does not exist.** The
  field-code generator emits `φ`, `Eᵢ` and `dEᵢdxⱼ`, never a derivative of the potential, so every
  component raised an `UndefVarError`. It went unnoticed because `dH` had no callers until the
  vector field above started consulting it. The potential enters as `∂φ/∂xᵢ = -Eᵢ`.
- **`charged_particle_3d_sode` is restricted to cartesian equilibria.** At frozen position the
  curvilinear kick is *quadratic* in `v`, so the Boris/Cayley map is no longer its exact flow and
  the splitting has no exact solution map. Rather than integrate the cartesian system in a chart
  where that is a different model, the constructor throws unless the metric is trivial;
  `TokamakSmallNoncanonical` no longer exports it. The other three are unaffected.
- **`plot_fieldlines` accepted equilibria it is meaningless for.** It contours `equ.A₃(0,x,y,0)/x`,
  the poloidal flux `ψ = R A_φ`, which assumes an axisymmetric cylindrical chart; for the cartesian
  and toroidal equilibria it drew a plausible-looking wrong picture rather than failing, and
  `plot_trajectory_poloidal(R, Z, equ)` forwarded any `equ` to it. Both now throw an
  `ArgumentError`. The new `is_axisymmetric_cylindrical` decides this from the coordinate ranges the
  field code injects — the chart whose *third* coordinate is angular and whose first two are not —
  so it needs no hard-coded list and classifies all eleven 4D equilibria correctly.
- **`hamiltonian_u` of the 3D guiding centre omitted `φ`**, unlike the `hamiltonian` beside it,
  which made the two incomparable for an equilibrium with a potential. `TODO.md` recorded it as
  unused; it is not — four scripts plot it against `hamiltonian` as a measure of constraint drift,
  which is exactly the comparison the omission broke. This changes the value it reports for
  `QuadraticPotentials3d`, the one equilibrium where `φ ≠ 0`.
- Docstring cross-references that could not resolve, and `lodeproblem`, `odeproblem`,
  `iodeproblem` and `iodeproblem_λ` of the 4D guiding centre, which had no docstrings at all.
- **The explanation given for the cost of the 3D `hode` path was wrong in every part.** Both
  `docs/src/findings.md` and the 0.2.0 entry below attributed it to "nested `ForwardDiff`: the second
  derivatives of the 3D Hamiltonian differentiate field functions that are themselves AD-generated,
  and the backtracking line search re-evaluates them". Profiling contradicts all three clauses. There
  is no nesting and no AD in the field layer: `ElectromagneticFields` generates its field functions
  symbolically with SymEngine at precompile time and does not depend on `ForwardDiff`. `hodeproblem`
  never touches a second derivative of the Hamiltonian — its right-hand side is first derivatives
  throughout, and the Hessian belongs to `hodeproblem_canonical` alone, hand-written in closed form
  rather than differentiated. And the line search is 8 % of wall clock against the Jacobian's 13 %:
  of the roughly ninety-eight right-hand side evaluations per step, **four** are on `Dual`.

  Only the clause about `MidpointExtrapolation(5)` survived, and it was the whole of the first layer
  rather than a multiplier on top of another one. `findings.md:280` already stated correctly that
  `ForwardDiff` is not a dependency, contradicting the claim eighty lines below it. The passage is
  rewritten with the profile behind it. The 0.2.0 entry is left as it stands, being released.

### Known issues

- `lodeproblem_formal_lagrangian` is still a verbatim copy of `lodeproblem`. It was renamed with the
  rest but not fixed: it should be the formal Lagrangian of the equations of motion in noncanonical
  Hamiltonian form, `L(z, y, ż) = yⁱ (Ωᵢⱼ(z) żʲ + ∂ᵢH(z))` on the doubled state `(z, y)`, which is
  an eight-dimensional problem rather than the four-dimensional one it builds. See `TODO.md`.
- **The equilibrium modules hard-code `u` and `μ` and their provenance is unknown.**
  `src/utils/initial_conditions.jl` already provides `InitialConditions`, which derives both from a
  position, a pitch angle and an energy, along with per-model converters — and no equilibrium module
  calls any of it; the only callers are two documentation pages and one script. The literals are not
  reproducible from those functions either: on the small tokamak in cartesian coordinates the
  parallel velocity comes out at `0.00042987` against a hard-coded `0.00045136`, and `|v⊥|²/2B` at
  `2.31459E-6` against a hard-coded `2.310E-6` — the latter a nearby but different number, not the
  same one rounded. The docstrings added in this release say the provenance is unrecorded rather
  than inventing one. Routing every module through a single initial-conditions layer is the fix and
  is filed in `TODO.md`; it would shift the shipped values by those percentages, so it wants its own
  diff.


## [0.2.0] - 2026-07-25

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
- **A findings page** (`docs/src/findings.md`) and **`TODO.md`**, together with the six
  `scripts/study_*.jl` that reproduce the former. The findings page records the numerical studies
  made during the audit that characterise the models rather than merely verify a line of code: the
  conditioning of the 3D guiding centre constraint pair, which expression is the conserved toroidal
  momentum, the exact proportionality between the gyrokinetic and the guiding centre vector fields,
  the volume preservation of the gyrokinetic splitting measured against step size, and the
  conservation levels of each family side by side. `TODO.md` carries the open work.
- **Diagnostics for the 3D guiding centre model.** `guiding_center_3d_diagnostics.jl` was an empty
  file, included by all thirteen equilibria, so the model had no `compute_energy`,
  `compute_energy_error` or `compute_toroidal_momentum` at all. Alongside them,
  `compute_constraints` returns all three constraints ``g^1``, ``g^2``, ``g^3`` as data series: they
  vanish along the exact flow, so their magnitude is the drift off the manifold on which `p = ϑ`,
  which is the quantity the approximately symplectic methods are designed to keep bounded. All three
  are reported rather than the two the problem's constraint pair retains, because the omitted one
  carries the drift of the retained pair amplified by ``1/b_m`` — a factor of 800 over a thousand
  steps on `TokamakMediumCartesian` — and is invisible otherwise. `SolovevIter`, whose
  `include` of the diagnostics was commented out, has them too now.
- **A documentation page for `GyroKinetics4d`** (`docs/src/gyro_kinetics_4d.md`), which had none
  and was not listed in the index. It derives the gyrokinetic characteristics, the rescaled time
  and the ``\beta``/``\gamma`` potentials from the notes the module follows, explains the six-way
  splitting and why it is volume preserving, and records what is and is not implemented.

### Changed

- **Five 3D guiding centre integration testsets are commented out**, leaving two:
  `TokamakMediumCylindrical`, `TokamakSmall{Cartesian,Cylindrical,Toroidal}` and
  `SolovevSymmetricField`. All five have `b₁ = 0` exactly at their shipped initial condition, where
  the multipliers of the `(g³, g¹)` constraint pair are infinite and `hode` cannot be started at
  all — see *Known singularities* in the audit. They could not have run before this branch either:
  `initial_conditions_*` raised a `MethodError` in three of them. The cross product of families and
  equilibria is covered instead by `test/structure_tests.jl`, which does not integrate. Restoring
  them is the top item in `TODO.md`. `test/guiding_center_3d_tests.jl`.
- **`pauli_particle_3d_pode` and `pauli_particle_3d_hode` default to `tstep = Δt` rather than
  `Δt / 100`**, matching `pauli_particle_3d_iode` and every other constructor in the package. Any
  caller relying on the old default now takes a hundredfold larger step.
- The test dependencies moved from `[extras]`/`[targets]` in the root `Project.toml` to
  `test/Project.toml`, which Pkg uses in preference as soon as it exists; the old section was dead
  and had drifted out of sync with it. The `GeometricIntegrators` compat bound moved with them.
- `scripts/Project.toml` no longer sources `ElectromagneticFields` from `../../ElectromagneticFields`.
  That path is outside the repository and does not exist, so `julia --project=scripts` — the command
  `docs/src/findings.md` gives for reproducing every study — failed to instantiate. The registered
  0.6 satisfies the bound.
- The two 3D guiding centre integration blocks run 100 steps rather than 1000. They only assert that
  `integrate` returns a solution, which 100 steps exercises identically, and the 3D `hode` path is
  expensive per step: its second derivatives differentiate field functions that are themselves
  AD-generated, and `MidpointExtrapolation(5)` triples that. Long-time behaviour of the 3D model is
  checked by the `3D guiding centre diagnostics` block in `test/structure_tests.jl`, which keeps its
  1000 steps and its conservation assertions. Note this cost is not specific to the ITER-like
  equilibria — the medium tokamak is slower still. `test/guiding_center_3d_tests.jl`.
- CI collects coverage on a single matrix job instead of all twelve. Instrumentation costs a factor
  of 2.4–6.7 on these models and only one job's results were ever uploaded.
  `.github/workflows/CI.yml`.
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

- **The test suite requested a solver tolerance the ITER-scale models cannot reach**, which made
  every ITER-like Solov'ev block take minutes instead of seconds and emitted over fifty thousand
  suppressed solver warnings per run. `SimpleSolvers` accepts an iterate when
  `rfₐ ≤ f_abstol + f_reltol ‖F(x₀)‖`; the tests set *both* terms to `1E-15`, which collapses the
  criterion to `1E-15` regardless of problem scale and in particular removes the relative term that
  exists so a large-magnitude solve can converge. The 4D guiding centre residual carries the
  momentum equation `p = ϑ(q)` with `‖ϑ‖ ≈ 15` in an ITER field, so no residual can go below the
  round-off floor `‖ϑ‖ eps ≈ 3.3E-15`; the solver then had no reachable stopping criterion and ran
  to its 1000-iteration limit on nearly half of all steps. The tests now ask for `f_abstol = 1E-12`,
  leave `f_reltol` alone, and cap `max_iterations` (with `warn_iterations` set to match, so that
  hitting the cap is still reported). Solutions are unchanged to ~1E-14.
  `test/guiding_center_3d_tests.jl`, `test/guiding_center_4d_tests.jl`,
  `test/pauli_particle_3d_tests.jl`, `test/gyro_kinetics_4d_tests.jl`.
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
  its own two-form, `Ωᵢⱼ = ∂ϑᵢ/∂qʲ - ∂ϑⱼ/∂qⁱ` on the three-dimensional configuration space.
  It is *not* the `6 × 6` canonical form on `(q,p)`: `GeometricEquations` sizes the buffer it hands
  the `LODE` two-form from the problem dimension, `D = length(q) = 3`, so the canonical form does
  not fit and would have raised a `BoundsError` at the first evaluation.
- **`ω` and `dϑ` of the noncanonical charged particle referred to `dϑ₁dx₂` and siblings that are
  defined nowhere**, and both had their output argument in the wrong slot for the
  `GeometricEquations` convention, so the `LODEProblem` built from `ω` could never be evaluated.
  The derivatives are now defined, including the metric terms, and both functions take the output
  first. Their velocity block is `∂ϑᵢ/∂vʲ = gᵢⱼ`; it had been hard-coded to the identity in `ω` and
  to zero in `dϑ`.
- **`charged_particle_3d_iode_g` of the noncanonical charged particle was written for a one-form
  without the metric**, and so described a different `ϑ` than the `dϑ` and `ω` beside it: the
  position block omitted the `∂gⱼⱼ/∂xⁱ vʲ` terms and the velocity block was the identity where
  `∂ϑⱼ/∂vⁱ = gᵢᵢ δᵢⱼ`. Being linear in `λ`, it passed the linearity check that caught the other
  two projection forces; it is now written as `dϑᵀ λ` through the same derivatives `dϑ` uses, and
  the tests compare all three models' `g` against central differences of their own one-form.
- **`charged_particle_3d_sode` could never have run**: its velocity map referred to undefined `q`,
  `t` and `Δt`, called `eye`, removed from Julia in 0.7, and built a rotation matrix whose signs
  were not those of `v × B`; the constructor passed an undefined `x₀`. Three further breaks were
  in the same constructor:
    * it dropped the electric field — the splitting was `ẋ = v` and `v̇ = v × B`, while the model
      is `v̇ = E + v × B`. `E` is now carried by the frozen-position half, which makes it the Boris
      push, `(1 - h/2 B̂) v₁ = (1 + h/2 B̂) v₀ + h E`; for `E = 0` that is the Cayley transform as
      before, an exact rotation preserving `|v|`. Every shipped equilibrium has `φ = 0`, so no
      result changes, but the splitting now covers the whole vector field for any that does not;
    * `SODEProblem(nothing, (fv, fx), …)` overflowed the stack. The typed constructor requires
      `v::Tuple`, and the fallback `SODEProblem(v, args...) = SODEProblem(v, nothing, args...)`
      recurses on itself when it does not match. The problem now supplies both the vector fields
      and their exact solution maps, index by index, which is the form `SODE` is happiest with;
    * the solution maps took `(q₁, t₁, q₀, t₀)`, one argument short of the
      `(q₁, t₁, q₀, t₀, params)` that `GeometricEquations` calls them with.

  The splitting is now covered by `test/structure_tests.jl`, which asserts that its two maps sum to
  `charged_particle_3d_v` — the Lie-Trotter consistency condition, and the assertion that fails if
  a term is ever dropped from one of them again.
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
- The audit and `TODO.md` named the wrong equilibria as affected by the cartesian-only limitation of
  `charged_particle_3d_noncanonical`. `SingularField` and `SymmetricField` are cartesian, as is
  `ThetaPinchNoncanonical`; the only one of the four that is curvilinear is
  `TokamakSmallNoncanonical`, which is toroidal rather than cylindrical. The caveat is narrower
  than it was documented to be.
- The audit, the findings page and `scripts/study_guiding_center_3d_conditioning.jl` all cited
  `scripts/guiding_center_3d.jl`, which was untracked and so existed only in one working tree. It is
  now part of the repository. The other dangling references were resolved by dropping the reference
  rather than adding the file; see *Removed*.
- `initial_conditions_dipole` returned a bare tuple where every other `GuidingCenter3d` module
  returns a `NamedTuple`.
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
  per precision.
- `src/guiding_center_4d/guiding_center_4d_problem.jl` and `docs/src/analytic.md`, both of which
  `TODO.md` referred to while they existed only in one working tree. The first was a four-line stub
  for a problem-type hierarchy that was never built; the second a handful of lines of commented-out
  code. Neither had content worth keeping, so the references went instead of the files becoming
  part of the repository.
- The references to `src/gyro_kinetics_4d/firk_with_coordinate_transformation.jl`, from `TODO.md`,
  the audit and the gyrokinetics page. That file is a 0.x-era sketch of a coordinate-transforming
  Runge-Kutta integrator, written against a `GeometricIntegrators` API (`Parameters`,
  `ODEIntegratorCache`, `create_nonlinear_solver`, `function_stages!`) removed several major
  versions ago — and the transformation was never finished: all three `coordinate_transformation_*`
  functions are commented out and the one that would apply the transform is an empty stub. Nothing
  in it survives the API turnover, so the two design decisions worth keeping — that the
  transformation is *implicit*, `q̃ = B*∥(q) q` solved for `q`, and that the stage vectors come in
  `(Q, Q̃)` pairs — are recorded in `TODO.md` instead.
- Fifteen unused transposed second derivatives of the 3D guiding centre Hamiltonian, two of which
  raised an `UndefVarError`, and the four `D²ϑ₁`–`D²ϑ₄` Hessians of the 4D guiding centre, whose
  only consumers were the `PoincareInvariants` trapezoidal variants that the 0.5 rewrite removed.
- The dead Travis badge from the documentation index.


## [0.1.0] - 2024-12-13

Initial release.


[Unreleased]: https://github.com/JuliaPlasma/ChargedParticleDynamics.jl/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/JuliaPlasma/ChargedParticleDynamics.jl/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/JuliaPlasma/ChargedParticleDynamics.jl/releases/tag/v0.1.0
