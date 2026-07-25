# TODO

Open work left by the correctness audit. The two sections under **Next** are the priority: they are
the parts of Phases 4 and 5 of the audit plan that were deferred.

Supporting material: [`docs/src/findings.md`](docs/src/findings.md) has the numerical studies that
motivate several of these, and `docs/src/audit.md` the user-facing summary of the known
limitations.


---

# Next

## Phase 4 — finish the API modernisation

`default_parameters` (D3) and the 3D guiding centre diagnostics (E4) are done. What remains breaks
every example, script and notebook, which is why it was held back for its own diff.

### Rename the constructors to the GeometricProblems scheme

`hode`, `hode_canonical`, `ode`, `guiding_center_4d_ode`, `charged_particle_3d_iode`,
`pauli_particle_3d_pode`, … → `odeproblem`, `hodeproblem`, `iodeproblem`, `lodeproblem`,
`podeproblem`, `sodeproblem`. Keep the model-distinguishing suffix only where two formulations
coexist in one module: `hodeproblem` vs `hodeproblem_canonical`.

Seventeen constructors across four naming schemes; `GeometricProblems`'
`massless_charged_particle_common.jl` is the reference.

### Keyword and constant names

`tspan`/`tstep` → `timespan`/`timestep`, and `const tspan`/`const Δt` →
`DEFAULT_TIMESPAN`/`DEFAULT_TIMESTEP`. There are 19 instances of the literal
`tspan=tspan, tstep=Δt`, where the keyword default shadows a module constant of the same name;
renaming the constants removes that hazard as a side effect.

### Uniform `initial_conditions_*` return shape

Currently a mix of `NamedTuple` (3D guiding centre, now uniformly so) and bare tuple (4D guiding
centre, gyrokinetics, Pauli). `test/structure_tests.jl` has to branch on `ics isa NamedTuple` to
sweep across families. Settle on the `NamedTuple`.

No deprecation shims — the package is pre-1.0 with no external users. Note the renames in the
CHANGELOG.


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

There was a 0.x-era sketch of this, `firk_with_coordinate_transformation.jl`, written against a
`GeometricIntegrators` API (`Parameters`, `ODEIntegratorCache`, `create_nonlinear_solver`,
`function_stages!`) removed several major versions ago. It is not in the repository: nothing in it
survived the API turnover, and the transformation itself was never finished — all three
`coordinate_transformation_*` functions were commented out and the one that would have applied the
transform was an empty stub. The one design decision worth carrying over is recorded here instead:

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

### More equilibria for `GyroKinetics4d`

Only the ITER-like Solov'ev equilibrium with X-point is provided, where the other families cover a
dozen each. `gc_solovev_iter_xpoint.jl` is the template and the module structure now supports more.


---

# Physics

## Add the remaining constraint formulations to the 3D guiding centre model

**Context from the author:** several formulations of the constraints exist. What is implemented is
one variant taken from an *earlier version* of Li, Zhang & Liu — so this is incompleteness, not a
mistake. The task is to generalise and add the missing variants.

This also happens to be the fix for the sharpest limitation in the package. The active pair is
``(g^3, g^1)``, which puts ``b_1`` in the denominator of both Lagrange multipliers, and **seven of
the eleven equilibria have ``b_1 = 0`` exactly at their shipped initial condition** — the
multipliers are infinite and `hode` cannot be started at all. See the table in
[findings.md](docs/src/findings.md); the pair ``(g^1, g^2)`` puts ``b_3 = b_\varphi`` there
instead, which does not vanish on the midplane.

Proposed shape: a `constraints = :g31` / `:g12` / `:g23` keyword threaded through `hode` and
`hode_canonical`, defaulting to the present `:g31` for compatibility, with a per-equilibrium default
chosen so each ships with a pair that is non-singular at its own initial condition.

The obstacle is that the constraints and *all* their first and second derivatives are written out
by hand for one pair — `guiding_center_3d_canonical.jl` is 700 lines of it. Three pairs would
triple that unless the derivatives are generated or the code restructured around
`Val{:g31}`-style dispatch. Worth deciding before starting.

Once it works, re-enable the commented-out 3D test blocks in `test/guiding_center_3d_tests.jl`
(~150 lines, all of them blocked on exactly this).

## Add the electrostatic potential to the 4D guiding centre Hamiltonian

`GuidingCenter4d` uses `H = ½u² + μ|B|` with no `φ` and no `E` in `∇H`
(`guiding_center_4d_common.jl`, `hamiltonian` and `dHdx₁₋₃`). `GuidingCenter3d` includes both, and
both reference papers carry `Φ(X)`, so the two guiding centre families currently describe different
systems whenever the equilibrium has a non-zero electrostatic potential.

Two lines, mirroring `guiding_center_3d_equations.jl`. What needs care:

* Which `ElectromagneticFields` configurations have a non-zero `φ`? For those the results change,
  so it is a behaviour change for the CHANGELOG and the test reference values may move.
* `GyroKinetics4d` inherits the same gap: its notes give the zero Larmor radius Hamiltonian with the
  `φ`, `A∥` and `⟨Ã∥⟩²` terms, but only `H₀` is implemented. Decide whether it gains the full
  `H^zlr` at the same time — a larger job, since `β` and `γ` acquire extra terms.
* `hamiltonian_u` in `guiding_center_3d_equations.jl` also omits `φ`, unlike `hamiltonian` beside
  it. It is unused; fix or delete.

## Make `charged_particle_3d_noncanonical` consistently curvilinear

It pairs a metric one-form, Hamiltonian, two-form and Jacobian with a plain cartesian Lorentz force
in `charged_particle_3d_v`, and a `dH` that does not differentiate the `hamiltonian` beside it.
Of the four equilibria that use the formulation, `SingularField`, `SymmetricField` and
`ThetaPinchNoncanonical` are cartesian and so unaffected; only `TokamakSmallNoncanonical`, which is
toroidal, is not self-consistent.

Recommended: derive the vector field from `Ω q̇ = -∇H` using the (now correct) `ω`, and add the
missing metric-derivative terms to `dH` and `charged_particle_3d_iode_f`. Restricting the
formulation to cartesian equilibria instead would mean deleting one equilibrium — cheaper than it
looked while the caveat was thought to cover two, so it is now a real alternative.
`charged_particle_3d_sode_fv`/`_vv` are cartesian for the same reason and would follow.


---

# Studies worth doing

## Long-time energy behaviour of the volume-preserving splitting

The gyrokinetic splitting preserves volume exactly but is not symplectic, and over a realistic span
its energy error is already four orders of magnitude worse than Gauss(2) on the same field
(7.7e-4 against 1.2e-7 — see [findings.md](docs/src/findings.md)). Whether that error is *bounded*
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

## Order of the composition methods for the splitting

Only `Strang()` has been exercised. The notes mention higher-order compositions; measuring the
observed order of a few would be cheap and would validate the subsystem integrators.


---

# Smaller things

* `test/guiding_center_3d_tests.jl` has ~150 lines of commented-out test blocks. Every equilibrium
  they cover — `TokamakMediumCylindrical`, `TokamakSmall{Cartesian,Cylindrical,Toroidal}`,
  `SolovevSymmetricField` — has `b₁ = 0` at its initial condition, so all are blocked on the
  constraint-formulation item above. `Dipole3d` and `QuadraticPotentials3d` are *not* blocked and
  have no test block at all; they should get one now.
* Per-equilibrium docstrings. Most modules have a one-line description and state neither the
  coordinate system, nor the equilibrium parameters, nor where the hard-coded `μ` and `u` of their
  `initial_conditions_*` came from (e.g. `uᵢ = -0.00045135897235326736`).
* `docs/src/initialization.md` documents a non-standard pitch-angle convention
  (`W⊥ = W sin α` rather than `W sin²α`). Self-consistent and now flagged, but worth deciding
  whether to move to the textbook convention.
* The `κ`-dependent "dg" formulation (`guiding_center_4d_dg`) now builds and its `ḡ` is verified
  against finite differences, but it has never been integrated in a test. Add one, or record why
  not.
* `compute_invariant_error` returns a `(value, error)` tuple, which is easy to misuse. Worth a line
  in the diagnostics docstrings.


---

# From the review of the audit branch

Observations from reviewing the audit that were not defects and so were not fixed with it. Recorded
here rather than lost in a pull request thread.

* **`test/quiet_solver_warnings.jl` is described as a tripwire but cannot trip.** `runtests.jl`
  reports the suppressed-warning count with `@info` and nothing asserts on it, so a tolerance
  slipping back below a residual floor would show up only if somebody read the log. It would become
  a real tripwire with an assertion — `@test suppressed_warning_count() < N` for some generous `N`
  — but that needs a number that is stable across platforms first, and the count is
  platform-dependent by construction, which is the whole reason the warnings are filtered. Worth
  deciding: either pick a loose bound and assert it, or drop the word "tripwire" from the comment.
* **Integration coverage of `GuidingCenter3d` is two equilibria of eleven**, and the structure
  tests, which cover all of them, do not integrate. This is a consequence of the `b₁ = 0`
  singularity rather than of anything the audit did — the other nine could not be started before it
  either — but it means the 3D model's *dynamics* is exercised on `SolovevIterXpoint` and
  `TokamakMediumCartesian` alone. Fixing the constraint pair fixes the coverage too; until then it
  is worth knowing how thin it is.
* **`plot_fieldlines` assumes an axisymmetric cylindrical equilibrium.** It contours
  `equ.A₃(0, x, y, 0) / x`, i.e. the poloidal flux ψ = R A_φ, which is meaningless for the cartesian
  and toroidal equilibria. Nothing stops it being called on them — `plot_trajectory_poloidal(R, Z,
  equ)` forwards any `equ` — and it draws a plausible-looking wrong picture rather than failing.
  Either restrict it or compute ψ from the coordinate system.
* **The SciML reformat landed in the same commits as the semantic changes**, which makes several
  equilibrium modules read as whole-file rewrites; `git diff -w` is the way to review them. Worth
  keeping formatting sweeps to their own commit next time, now that `.JuliaFormatter.toml` exists
  and the one-off reflow is behind us.


---

# Resolved

Recorded so a future session does not re-open them.

**Constraint pair.** Not a mistake — the implemented variant comes from an earlier version of the
paper, and the others are simply missing. Generalising is a planned task; see above.

**`GyroKinetics4d` coordinate transformation.** It was intended as preconditioning for the
nonlinear solver. Removing it from the default problem constructors was correct; it belongs inside a
coordinate-transforming `IRK` integrator. See above.

**Normalisation.** The ε-ordered strongly-magnetised form in `docs/src/normalization.md` is for
reference only; no model exposes the `(l̂/ρ̂_th)` and `(eφ̂/mv_th²)` prefactors and none is expected
to. Noted in that page.

**Off-diagonal metrics.** Diagonal-only is fine for now. Generalising would complicate the code
substantially and will wait until there is a need. Noted prominently in `docs/src/audit.md`.
