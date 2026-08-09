# TODO

Open work left by the correctness audit.

Supporting material: [`docs/src/findings.md`](docs/src/findings.md) has the numerical studies that
motivate several of these, and `docs/src/audit.md` the user-facing summary of the known
limitations.

**Phase 4 is done** — the constructors follow the `GeometricProblems` scheme
(`odeproblem`/`iodeproblem`/…), the keywords are `timespan`/`timestep` against
`DEFAULT_TIMESPAN`/`DEFAULT_TIMESTEP`, `parameters` is a keyword defaulting to
`default_parameters()`, and every `initial_conditions_*` returns the same named tuple. So are the
electrostatic potential, the curvilinear noncanonical charged particle, and the additional
`GyroKinetics4d` equilibria. See the CHANGELOG.


---

# Migration: SimpleSolvers 0.10 / GeometricIntegrators 0.17 — done

Released as part of 0.4.0; see the CHANGELOG. `test/`, `docs/` and `scripts/` are on GeometricIntegrators
0.17.0, SimpleSolvers 0.10.1 and GeometricIntegratorsBase 0.5.1, each bounded in the environment
that resolves it.

The suite now emits **no** nonlinear-solver warnings and `runtests.jl` asserts that per test file,
which replaces both the un-asserted count this section wanted repurposed and the `@test_nowarn`
that was abandoned as platform-dependent. The 1853 messages it started from had three causes, all
real and all fixed: `GuidingCenter4d.SolovevSymmetricField` running its `ode` orbits at a step size
where they could not converge (1616), one Pauli `iodeproblem()` left at the module default
`Δt = 500` (174), and the unconfigured calls exercising the old flat `f_abstol = 8eps()` library
default (the rest), which `GeometricIntegratorsBase` 0.5.1 fixed upstream.

Two claims in this section were wrong and are worth not repeating:

* `status`/`isstalled` are **not** exported by SimpleSolvers — deliberately, since exporting
  `status` would be hostile to `using SimpleSolvers`. They have to be qualified, and reading them
  means driving the solver a step at a time rather than calling `integrate`.
* The residue was said to be irreducible until `GeometricIntegratorsBase` changed its default. It
  did change it, in 0.5.1, to `max(8, solversize(method, problem)) * eps(datatype(problem))` — but
  that only silenced the unconfigured calls. The bulk was this repository's own two mis-set time
  steps.

What did **not** change: the `‖ϑ‖ eps` reasoning behind `f_abstol = 1E-12` still holds, because the
sized default is 1.8E-15–2.7E-15 and the ITER floor is 3.3E-15. `max_iterations = 50` is kept: the
solves that reach it are genuinely still iterating rather than stalling, so `max_stalls` does not
catch them.

## Retracted claim worth leaving retracted

- [ ] `CHANGELOG.md` and `docs/src/findings.md` record a *retracted* attribution of the 3D `hode`
      cost to the line search — profiling showed it was 8 % of wall clock against the Jacobian's
      13 %. That retraction is still correct and was not revisited for 0.10.

---

# Next

## Phase 5 — finish porting the stale subsystems

The `PoincareInvariants` layer, the `Dec128` removal and the `GyroKinetics4d` equations of motion
are done. What remains:

### A coordinate-transforming IRK integrator

**The integrator is the right home for the coordinate transformation.** The scaling
``\tilde q = \omega_0 q`` in `coordinate_transformations.jl` is a *preconditioner* for the nonlinear
solver, not a change of variables in the model, which is why it is no longer applied by the problem
constructors. It belongs inside the solver, and nowhere else.

Write it against the current `IRK` integrator: `GeometricIntegrators/src/integrators/rk/
integrators_irk.jl` — `struct IRK{TT<:Tableau,ImplicitUpdate} <: IRKMethod`, `IRKCache`,
`initmethod`, and the `components!`/`residual!` split.

The 0.x-era sketch of this is tracked at
[`src/gyro_kinetics_4d/irk_with_coordinate_transformation.jl`](src/gyro_kinetics_4d/irk_with_coordinate_transformation.jl).
It is not `include`d by any module and does not compile: it is written against a
`GeometricIntegrators` API (`Parameters`, `ODEIntegratorCache`, `CacheDict`,
`create_nonlinear_solver`, `function_stages!`, `AtomicSolutionODE`) removed several major versions
ago, and the transformation itself was never finished — all three `coordinate_transformation_*`
functions are commented out and the one that would apply the transform is an empty stub. It is on
file as a starting point, not as code to be revived.

What it gives the port:

* the `(Q, Q̃)` stage-pair cache layout of `IntegratorCacheFIRKwCT`;
* the `compute_stages!`/`function_stages!` split, which is `components!`/`residual!` in the current
  integrator;
* the commented-out `coordinate_transformation_rhs!`, which spells out the relation to be solved,
  ``b = \tilde q - B^{\star}_{\parallel}(q) \, q``.

Everything that touches the removed API has to be rewritten. Its internal names still use the old
`FIRK` spelling; rename them to `IRK` with the port. The two design decisions worth stating
explicitly:

* **The transformation is implicit.** The intended relation is ``\tilde q = B^{\star}_{\parallel}(q)
  \, q``, with ``B^{\star}_{\parallel}`` evaluated at the *unknown* ``q``, so recovering ``q`` from
  ``\tilde q`` needs its own nonlinear solve inside the step — it is not the explicit rescaling by a
  frozen ``\omega_0`` that `coordinate_transformations.jl` provides. That is the same point the
  audit makes about why freezing ``B^{\star}_{\parallel}`` at the initial condition does not
  reproduce the substitution.
* The stage vectors come in pairs, ``Q`` and ``\tilde Q``, with the vector field evaluated on the
  untransformed stage and the solve carried out on the transformed one.

Once it exists, add a test that it produces the same trajectory as plain `IRK` on the same problem
— the transformation must not change the solution, only the conditioning of the solve — and
ideally one showing the iteration count it saves.

### `lodeproblem_formal_lagrangian` is not the formal Lagrangian

`src/guiding_center_4d/guiding_center_4d_equations.jl`. It is a verbatim copy of `lodeproblem` —
the same `LODEProblem` from the same `ϑ`, `f`, `g`, `ω` and `lagrangian`. It should instead be the
**formal Lagrangian of the equations of motion in noncanonical Hamiltonian form**: writing them as
``F_i(z, \dot z) = \Omega_{ij}(z) \dot z^j + \partial_i H(z) = 0``, the formal Lagrangian is

```math
L(z, y, \dot z) = y^i \, F_i(z, \dot z)
```

on the doubled state ``(z, y)``, whose Euler-Lagrange equations return the system on variation of
``y`` and its adjoint on variation of ``z``. That is an **eight-dimensional** problem for the 4D
guiding centre, not the four-dimensional one it builds today.

`GeometricProblems.PointVortices.lodeproblem_formal_lagrangian` is **not** the template: there the
name is used for the degenerate phasespace Lagrangian ``L = \vartheta \cdot \dot q - H``, which is
exactly what `lodeproblem` already is here. Copying it would reproduce the present duplication
under a different justification.


---

# Physics

## ~~Add the remaining constraint formulations to the 3D guiding centre model~~ — done

All three constraint pairs are implemented, selected by a `constraints = :g31` / `:g12` / `:g23`
keyword on `hodeproblem` and `hodeproblem_canonical`, with a per-equilibrium `default_constraints`
that is regular at that equilibrium's own initial conditions. The constraints and every one of their
first and second derivatives are now written once against a generic index pair in
`src/guiding_center_3d/guiding_center_3d_constraints.jl`, as decided;
`guiding_center_3d_canonical.jl` lost some four hundred lines rather than tripling in size.

The compact form of Eq. (29) is implemented as well, as a third formulation, `hodeproblem_compact`.
Since that equation as printed is cartesian, it is derived coordinate-generally instead of ported;
the derivation is at the top of `guiding_center_3d_compact.jl` and written up in
[docs/src/guiding_center_3d.md](docs/src/guiding_center_3d.md). It turns out to be independent of the
constraint pair and free of every `bₘ` denominator.

Five of the seven commented-out 3D test blocks are re-enabled. The two that are not —
`SymmetricField` and `ThetaPinchField` — never had integration content: they called
`guiding_center_3d_loop_ode` and `guiding_center_3d_surface_ode`, which this package no longer
defines. See the note where those blocks used to be.

Two things this left open, neither of them about the constraint pair:

* **`hodeproblem_canonical` diverges on `SolovevSymmetricField`**, although the pair is well
  conditioned there (``\lambda_o \approx -228``) and the same orbit runs to completion under
  `hodeproblem` and `hodeproblem_compact` at 4.5e-12 and 2.4e-11. The
  ``\partial\lambda/\partial q``, ``\partial\lambda/\partial p`` terms amplify the constraint drift
  rather than ignoring it. It is the only equilibrium where this happens. Note that it does not
  *stop*: over a hundred steps it reaches ``\vert \Delta H/H \vert`` of 1.5e+228 and
  ``\max \vert g^k \vert`` of 2.6e+114 at 44.2 Newton iterations per solve against 1.1 for the other
  formulations, so it reports a destroyed trajectory as a number rather than raising. Anything that
  screens these formulations automatically has to test the magnitude, not just for an exception.
  This equilibrium is a cartesian chart, so the 0.7.0 field reversal is not involved.
* **The omitted constraint drifts `1/bₘ` times as much as the retained pair**, since
  ``b_1 g^2 - b_2 g^1 + b_3 g^3 = 0`` determines it from them. On the ITER Solov'ev X-point that is a
  factor of 170. Nothing is wrong with it, but it is one of two reasons to prefer the pair whose
  ``b_m`` is largest rather than merely non-zero.
* **`default_constraints` is chosen on regularity at the initial condition, not on conditioning along
  the orbit** — with two exceptions, both now settled. `Dipole3d` had to be moved to `:g12` because its
  orbit takes both `b₁` and `b₂` through zero and the other two pairs lose the orbit entirely by
  ``t = 30``. `TokamakMediumCartesian` moved to `:g23` when its initial condition was aligned to
  `y = 0` to match its 4D and Pauli counterparts: that puts it on the `b₁ = b_x = 0` midplane, where
  `:g31` is singular rather than merely ill-conditioned, and `b₂ = 0.9923` makes `:g23` the best
  conditioned of the two that remain (``\lambda_o = -3.85`` against `:g12`'s ``0.48``).

  The remaining nine were left where regularity put them. A sweep of ``\min |b_m|`` along each shipped
  orbit would settle them properly; note that changing one moves numbers findings.md tabulates, and
  that the `TokamakMediumCartesian` move also cost that equilibrium its role as the one place all three
  pairs could be compared — which is now `SolovevIterXpoint`, and is a better home for it, the three
  pairs agreeing there to 8E-16 against the cartesian chart's 6E-9.

## ~~The per-step cost of the 3D guiding centre model~~ — done

Never tracked as an item: `docs/src/findings.md` recorded the 3D `hodeproblem` at two orders of
magnitude more per step than the Pauli model, explained it as "nested `ForwardDiff`: the second
derivatives of the 3D Hamiltonian differentiate field functions that are themselves AD-generated,
and the backtracking line search re-evaluates them", and left it there.

That explanation was wrong in every part, and profiling says so. `ElectromagneticFields` generates
its field functions symbolically with SymEngine at precompile time and does not depend on
`ForwardDiff`; `hodeproblem` never evaluates a second derivative of the Hamiltonian, only
`hodeproblem_canonical` does, and those are hand-written closed forms rather than differentiated;
and of the roughly ninety-eight right-hand side evaluations per step, four are on `Dual`, the line
search being 8 % of wall clock against the Jacobian's 13 %.

There were three real costs. Two are fixed here. `MidpointExtrapolation(5)`, which is not the default
for either method, was 70 % of wall clock and is no longer requested. The right-hand side no longer
re-enters the generated field code once per bracket term: `fieldvalues` evaluates each injected
function once into a `FieldValues`, and `SecondFieldValues` does the same for the second derivatives
the canonicalised form needs.

The third was the expression swell in the generated code itself, and lives upstream. It is delivered:
[JuliaPlasma/ElectromagneticFields.jl#10](https://github.com/JuliaPlasma/ElectromagneticFields.jl/pull/10)
is merged and released as **0.6.3**, which `Project.toml` now requires. It takes `db₁dx₁` from 549 ns
to 57 and `d²b₁dx₁dx₁` from 1733 to 80, and changes no value: all 27456 field values this package
injects across its 32 equilibrium modules are bit-identical between 0.6.2 and 0.6.3, and so is every
trajectory integrated from them.

`hodeproblem` on `SolovevIterXpoint` goes from 1.59 ms/step to 0.082 and `hodeproblem_canonical` from
19.9 to 0.39. For `hodeproblem` the two contributions are separable and multiplicative — 6.9× from
this package, 2.8× from the release — and for the canonicalised form they reinforce each other,
because it is the one dominated by second derivatives and those gained the most. `docs/src/findings.md`
has the four-way measurement that separates them.

One trap this left behind, worth not re-creating: `scripts/Manifest.toml` had a
`path = "../../ElectromagneticFields"` stanza with no matching `[sources]` entry in
`scripts/Project.toml`, so `julia --project=scripts` silently resolved the dependency to a working
checkout — and `Pkg.status` reports the manifest's version string rather than the version at the path,
so it could not be noticed that way. One round of the measurements above was taken against a feature
branch before this was found. The manifest is deleted; check `Base.locate_package` rather than
`Pkg.status` when a timing looks surprising.

## One initial-conditions layer for all models

**The equilibrium modules should not carry `u` and `μ` as literals.** Each one should declare its
initial conditions in *physical* terms — a position, a pitch angle and an energy — and let a shared
layer derive the state each model actually needs. Today all forty-odd modules hard-code the
derived numbers instead, and where those numbers came from is no longer known; see *Provenance*
below for why that is more than an aesthetic complaint.

### Most of the layer already exists

`src/utils/initial_conditions.jl` has the pieces and nothing uses them:

* `InitialConditions` — a struct holding the physical state once: `x`, `X`, `ρ`, `vvec`, `vpar`,
  `vper`, `v`, `u`, `μ`, `θ`, `α`, `ω`, `mass`, `energy`, `charge`.
* `InitialConditions(X, θ, α, E, M, C, â, b̂, ĉ, b, B, g̅, DF̄, J; l₀)` — builds it from a gyro-centre
  position, gyro angle, pitch angle and an energy in eV, doing the parallel/perpendicular split
  against the equilibrium's own field and metric. `InitialConditionsGC(X, θ, u, μ, …)` is the
  variant that starts from `(u, μ)` instead.
* `charged_particle(ics)`, `guiding_center(ics)`, `pauli_particle(ics)` — per-model converters.

**No equilibrium module calls any of it.** The only callers are `docs/src/initialization.md`,
`docs/src/examples/iter_cylindrical.md` and `scripts/guiding_center_3d.jl`. So there is a function
that computes `u` and `μ` from physical inputs, and forty modules that hard-code them anyway.

### What is missing

* **A 3D guiding centre converter.** That model's state is `(q, p)` with `p = ϑ(q)`, so it needs
  `guiding_center_3d(ics)` alongside the three above; at present the modules reach the momentum
  through the model-local `initial_conditions(tᵢ, Qᵢ)` in `guiding_center_3d_equations.jl`.
* **A gyrokinetic converter.** `GyroKinetics4d` shares the 4D state, so this is `guiding_center`
  plus the rescaled time step — the factor `B^{\star}_{\parallel}` that its `DEFAULT_TIMESTEP`
  carries and that must not be copied between equilibria.
* **The uniform return shape.** The three converters return bare tuples. Every constructor now
  takes `(q = …, params = …)`, plus `p` or `v` where the model has one, so the converters should
  return that directly and the modules should be one line each:

  ```julia
  initial_conditions_deeply_passing() = guiding_center(InitialConditions(x₀, θ, α, E, md, 1, …))
  ```

* **The per-equilibrium field arguments.** `InitialConditions` takes nine field functions
  positionally. Since they are all injected into the module by `@code`, a small wrapper per module —
  or a macro alongside the injection — would remove that boilerplate.

### Provenance: why the literals cannot simply be documented

Two different number styles sit side by side, which is the first sign they have different origins:

```julia
const uᵢ = -0.00045135897235326736                              # 17 digits: a Float64 round-trip
initial_conditions_pauli() = (q = [1.05, 0., 0., 4.3E-4], params = (μ = 2.310E-6,))   # 4 digits
```

The first is machine output that was printed and pasted; the second is hand-rounded. Neither is
reproducible from the package's own helpers. Taking the small tokamak in cartesian coordinates,
with the `xᵢ = [1.05, 0, 0]` and `vᵢ = [2.1E-3, 4.3E-4, 0]` its Pauli module declares:

| Quantity | What the helpers give | What is hard-coded | |
|---|---|---|---|
| ``v \cdot b`` | `0.00042987` | `0.00045136` | 5 % apart |
| ``\vert v_\perp \vert^{2} / 2B`` | `2.31459E-6` | `2.310E-6` | 0.2 % apart |

The second is the telling one: `2.31459E-6` does not round to `2.310E-6` — it rounds to
`2.315E-6`. It is not the same number at lower precision, it is a nearby but different number. So
the literals are near-misses against every derivation this repository can perform, and one of the
following is true, with no way to tell which from the code:

* they were computed for a slightly different position or field configuration than the modules now
  use;
* they predate a change of normalisation or of the parallel/perpendicular split;
* they were derived outside the repository and transcribed.

The consequences are practical rather than cosmetic. Nobody can produce the corresponding numbers
for a *new* equilibrium without redoing whatever the original derivation was. Nobody can say
whether "barely passing" versus "deeply passing" reflects a computed trapped-passing boundary or a
value that was tuned until the orbit looked right. And a genuine transcription error would be
indistinguishable from the discrepancies above.

The docstrings therefore say the provenance is unrecorded rather than inventing one. **Building the
layer above resolves this by making the question moot**: the physical inputs are the declaration,
and `u` and `μ` become derived quantities that any future equilibrium gets for free. Note that
re-deriving them will shift the numbers by the percentages in the table, so test reference values
move and it wants its own diff.

## The full zero Larmor radius Hamiltonian for `GyroKinetics4d`

`H₀ = ½u² + μ|B| + φ` is implemented. The notes the module follows give the zero Larmor radius
Hamiltonian, which carries ``A_\parallel`` and ``\langle \tilde A_\parallel \rangle^2`` terms on
top of that; adding them changes `β` and `γ` as well, so it is a larger job than the potential was.


---

# Studies worth doing

## Long-time energy behaviour of the volume-preserving splitting

The gyrokinetic splitting preserves volume exactly but is not symplectic, and over a realistic span
its energy error is already five orders of magnitude worse than Gauss(2) on the same field
(3.7e-5 against 2.7e-10 — see [findings.md](docs/src/findings.md)). Whether that error is *bounded*
or *drifts* is unanswered and decides whether the scheme is preferable to a symplectic method here.

Integrate for 10⁵–10⁶ steps and compare the drift of the Strang composition against Gauss(2) and
against the 4D guiding centre model. This is also the thinnest section of the notes the module
follows.

## Constraint drift under `hode` vs `hode_canonical`

The argument of Li, Zhang & Liu is that canonical Runge-Kutta applied to their compact form keeps
`g¹` and `g²` bounded, and that this is what makes the method approximately symplectic. With
`compute_constraints` now available and the Hessian bug fixed, this can be measured directly:
constraint drift and energy error for both formulations, over long runs, at several step sizes.
Their §5 runs the same experiment on a toy field, a dipole and a tokamak — all three of which this
package has.

## The first-order slow manifold for the Pauli initial conditions

`scripts/study_model_agreement.jl` §3 shows that what separates the Pauli orbit from the guiding centre
one is almost entirely *where the initial condition is put*, not the dynamics. `initial_conditions(x, u, μ)`
sets `v = u b⃗`, the lowest-order slow manifold; the guiding centre also carries the ∇B and curvature
drifts, so the Pauli particle starts off the manifold by the drift velocity and gyrates about the right
orbit for ever after. Supplying the guiding centre velocity instead collapses the spread across
equilibria from four decades to two, and is worth 650–1080× on the ITER Solov'ev.

It is only half the correction, and the small tokamak rows prove it: there it is 2–10× *worse*,
consistently. The velocity mismatch is fixed but the *position* mismatch is not — the Pauli state is the
particle position, which differs from the guiding centre position by the gyroradius vector. The first-order
slow manifold needs both.

Worth doing as a `initial_conditions_slow_manifold(x, u, μ)` that applies both corrections, with the
drift and gyroradius written out locally rather than by depending on `GuidingCenter4d`. It would make the
Pauli family directly comparable with the guiding centre one at whatever order is wanted, and would give
the `test/model_agreement_tests.jl` threshold a reason to be tight rather than merely measured.

## Should `ωabs` carry the chart's orientation? — answered: no, and it was a bug

`ElectromagneticFields` 0.7.0 exposed that `GyroKinetics4d.ωabs` contracted the coordinate curl
`∂₂ϑ₃ - ∂₃ϑ₂`, which carries the *signed* Jacobian, so it came out `det(DF) · B*∥` — negative in the six
left-handed-chart equilibria. This section used to record that as a benign reparametrisation, on the
grounds that the splitting stays volume preserving and `v_gk = ωabs v_gc` still holds. Both are true and
neither is the point.

The factor is the phasespace Jacobian `√det Ω = J B*∥`, absorbed into the distribution function as
`f̃ = ωabs f`. It is a measure density and cannot be negative — as the name says. Carried into the vector
field, the sign made the model **chart-dependent**: the same physical particle circulated one way in the
cartesian chart of the small tokamak and the other way in its cylindrical chart. Nothing in the suite
caught it because every other assertion was a magnitude.

Resolved by applying the chart's `orientation()` — generated into each module by its `@code` call,
since `ElectromagneticFields` 0.7.1 — in `ωabs` and `v`. The trade-off
this section posed — divergence-free coordinate right-hand side *versus* chart-independent `dt` — was
not real: negating a divergence-free field leaves it divergence-free and each subsystem Hamiltonian, so
the splitting is untouched. `s` now runs with `t` in every chart, and the already-declared positive
`DEFAULT_TIMESTEP`s are correct, which under the signed factor they were not.

Rescaling by the *physical* `B*∥` instead — the other option this section listed — remains rejected, and
now for a stated reason: it is a scalar rather than a density, so the right-hand side would no longer be
a pure coordinate curl and would lose the divergence-free property that is the module's reason to exist.

`test/gyro_kinetics_4d_tests.jl` asserts `ωabs > 0` for all eight, checks each chart's generated
`orientation()` against the bare contraction this package computes, and compares the physical
toroidal velocity across all three charts of the small tokamak — the one statement that is not a
magnitude.

## Order of the composition methods for the splitting

Only `Strang()` has been exercised. The notes mention higher-order compositions; measuring the
observed order of a few would be cheap and would validate the subsystem integrators.


---

# Smaller things

* Per-equilibrium docstrings: the coordinate system and the equilibrium parameters are now stated
  everywhere. The provenance of the hard-coded `μ` and `u` is not, because it is unknown — see
  *One initial-conditions layer for all models* above, which subsumes this. If the derivation
  exists somewhere outside the repository, recording it would settle the question without the
  refactor.
* `docs/src/initialization.md` documents a non-standard pitch-angle convention
  (`W⊥ = W sin α` rather than `W sin²α`). Self-consistent and now flagged, but worth deciding
  whether to move to the textbook convention. Changing it perturbs every shipped initial condition,
  so it wants doing deliberately and on its own.
* `SolovevSymmetricField` has no `GyroKinetics4d` counterpart, and cannot have one as things stand:
  `SolovevSymmetric.@code` injects the equilibrium parameters `α` and `β` into the module, and `β`
  collides with the vector potential `β` of `gc_common.jl`. Fixing it means renaming one of the
  two — the model's `β` comes from the notes and is exported, the field's from
  `ElectromagneticFields` — so neither rename is local to this package.


---

# From the review of the audit branch

Observations from reviewing the audit that were not defects and so were not fixed with it. Recorded
here rather than lost in a pull request thread.

* **Integration coverage of `GuidingCenter3d` was four equilibria of eleven** — the remaining seven
  were blocked by the `b₁ = 0` singularity of the one implemented constraint pair rather than by
  anything the audit did. Selecting the pair per equilibrium lifted that: nine of the eleven now
  integrate, and the two that do not (`SymmetricField`, `ThetaPinchField`) ship Poincaré-invariant
  loops rather than initial conditions and never had integration tests to begin with.
* **The SciML reformat landed in the same commits as the semantic changes**, which makes several
  equilibrium modules read as whole-file rewrites; `git diff -w` is the way to review them. Worth
  keeping formatting sweeps to their own commit next time, now that `.JuliaFormatter.toml` exists
  and the one-off reflow is behind us.


---

# Resolved

Recorded so a future session does not re-open them.

**Constraint pair.** Not a mistake — the implemented variant came from an earlier version of the
paper, and the others were simply missing. All three now exist, on `Val{:g31}`-style dispatch as
settled, along with the compact form of Eq. (29); see above. There is nothing left to decide here
except which pair each equilibrium should *prefer* among those that are regular for it.

**The suppressed-warning count is not a tripwire.** It is reported by `@info` and nothing asserts
on it. A threshold was considered and rejected: which orbits exhaust the solver's iteration budget
is platform-dependent by construction, which is the whole reason the warnings are filtered, so any
bound tight enough to catch a regression would also flake. The comments in
`test/quiet_solver_warnings.jl` and `test/guiding_center_3d_tests.jl` now say so instead of
claiming otherwise.

**The noncanonical charged particle's `SODE` stays cartesian.** Making the model curvilinear made
the frozen-position kick quadratic in `v`, so the Boris/Cayley map is no longer its exact flow. The
constructor now throws for a non-trivial metric rather than silently integrating the cartesian
system; `TokamakSmallNoncanonical` is the one equilibrium affected. Generalising the splitting
would mean giving up either exactness or the Boris structure, which is a different design, not a
fix.

**`GyroKinetics4d` coordinate transformation.** It was intended as preconditioning for the
nonlinear solver. Removing it from the default problem constructors was correct; it belongs inside a
coordinate-transforming `IRK` integrator. See above.

**Normalisation.** The ε-ordered strongly-magnetised form in `docs/src/normalization.md` is for
reference only; no model exposes the `(l̂/ρ̂_th)` and `(eφ̂/mv_th²)` prefactors and none is expected
to. Noted in that page.

**Off-diagonal metrics.** Diagonal-only is fine for now. Generalising would complicate the code
substantially and will wait until there is a need. Noted prominently in `docs/src/audit.md`.
