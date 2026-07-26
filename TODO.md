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
by hand for one pair — `guiding_center_3d_canonical.jl` is 700 lines of it.

**Decided: restructure around `Val{:g31}`-style dispatch** rather than hand-writing the other two
pairs or generating the derivatives symbolically. Writing the constraints and their derivatives
once against a generic index pair shrinks that file instead of tripling it, at the cost of touching
every line of it. Do not open this again as a three-way choice.

Once it works, re-enable the commented-out 3D test blocks in `test/guiding_center_3d_tests.jl`
(~150 lines, all of them blocked on exactly this).

## The full zero Larmor radius Hamiltonian for `GyroKinetics4d`

`H₀ = ½u² + μ|B| + φ` is implemented. The notes the module follows give the zero Larmor radius
Hamiltonian, which carries ``A_\parallel`` and ``\langle \tilde A_\parallel \rangle^2`` terms on
top of that; adding them changes `β` and `γ` as well, so it is a larger job than the potential was.


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
  constraint-formulation item above.
* Per-equilibrium docstrings: the coordinate system and the equilibrium parameters are now stated
  everywhere, but the provenance of the hard-coded `μ` and `u` of the `initial_conditions_*` is
  still unknown (e.g. `uᵢ = -0.00045135897235326736`) and the docstrings say so rather than
  inventing one. If the derivation exists somewhere outside the repository, record it.
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

* **Integration coverage of `GuidingCenter3d` is four equilibria of eleven** — `SolovevIterXpoint`,
  `TokamakMediumCartesian`, and now `Dipole3d` and `QuadraticPotentials3d`. The structure tests
  cover all eleven but do not integrate. The remaining seven are blocked by the `b₁ = 0`
  singularity rather than by anything the audit did; fixing the constraint pair fixes the coverage
  too, and until then it is worth knowing how thin it is.
* **The SciML reformat landed in the same commits as the semantic changes**, which makes several
  equilibrium modules read as whole-file rewrites; `git diff -w` is the way to review them. Worth
  keeping formatting sweeps to their own commit next time, now that `.JuliaFormatter.toml` exists
  and the one-off reflow is behind us.


---

# Resolved

Recorded so a future session does not re-open them.

**Constraint pair.** Not a mistake — the implemented variant comes from an earlier version of the
paper, and the others are simply missing. Generalising is a planned task; see above. The shape is
settled too: `Val{:g31}`-style dispatch, not hand-written pairs or symbolic generation.

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
