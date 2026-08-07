# Findings

Numerical studies made while auditing the package that say something about the *models and
algorithms*, as opposed to merely confirming that a line of code is correct. The verification
checks that only did the latter became assertions in `test/structure_tests.jl` and are not
repeated here.

Each section names the script that reproduces it. They live in `scripts/` and run with

```
julia --project=scripts scripts/study_<name>.jl
```

All numbers below were produced on the state of the branch at the time of writing; they are
reproducible but not pinned by a test, so treat small differences as normal.

Every timing here depends on which `ElectromagneticFields` is resolved, because that package emits
the field functions the right-hand sides are built from and 0.6.3 made them several times faster than
0.6.2. Before quoting a number against a modified environment, check which one is actually loaded —
and do not check it with `Pkg.status`, which for a path dependency prints the version recorded in the
manifest rather than the version at the path:

```
julia --project=scripts -e 'println(Base.locate_package(Base.identify_package("ElectromagneticFields")))'
```

That resolves the load path without loading anything. A `scripts/Manifest.toml` carrying a
`path = "../../ElectromagneticFields"` stanza with no matching `[sources]` entry in
`scripts/Project.toml` is how one round of these measurements was silently taken against a working
checkout on a feature branch; deleting the manifest and re-resolving is the fix.


## Conditioning of the 3D guiding centre constraint formulations

`scripts/study_guiding_center_3d_conditioning.jl`

The constrained canonical formulation of Li, Zhang & Liu replaces the parallel-velocity constraint
by two of the three independent components of `v × b = 0`,

```
g¹ = b₁ v₃ - b₃ v₁ ,    g² = b₂ v₃ - b₃ v₂ ,    g³ = b₁ v₂ - b₂ v₁ ,
```

and the Poisson bracket in the denominator of both Lagrange multipliers carries a factor of the
`b`-component belonging to the constraint the pair leaves out: `+b₁` for `(g³, g¹)` (`constraints =
:g31`), `+b₃` for `(g¹, g²)` (`:g12`), `-b₂` for `(g², g³)` (`:g23`). The package once implemented
`(g³, g¹)` alone, so it was singular wherever `b₁ = 0`.

The minus on the third of those follows from that pair's ordering, not from the physics — reversing a
pair flips `λₒ` and both multipliers together — so only `|λₒ|` is meaningful. It is worth noting when
reading the table below, where `λₒ(:g23)` and `b₂` come out with opposite signs throughout.


### Which pairs are usable

At the shipped initial condition of each equilibrium. `λₒ = {g₁, g₂}` is the bracket both
multipliers divide by; `D` is the same quantity as the compact form computes it, as a weighted mean
over all three pairs that never divides by a single component of `b`.

| equilibrium | coords | b₁ | b₂ | b₃ | λₒ(:g31) | λₒ(:g12) | λₒ(:g23) | D | default |
|---|---|---|---|---|---|---|---|---|---|
| Dipole3d                 | cartesian   | -4.082e-01 | -8.165e-01 | 4.082e-01  | -3.402e+01 | 3.402e+01  | 6.804e+01  | 8.333e+01 | `:g12` |
| QuadraticPotentials3d    | cartesian   | 2.000e-03  | -3.000e-03 | 1.000e+00  | 2.000e-01  | 9.999e+01  | 3.000e-01  | 9.999e+01 | `:g31` |
| TokamakSmallCartesian    | cartesian   | **0**      | 9.997e-01  | 2.499e-02  | **0**      | 2.379e-02  | -9.516e-01 | 9.519e-01 | `:g23` |
| TokamakMediumCartesian   | cartesian   | -3.966e-02 | 9.914e-01  | 1.245e-01  | -1.537e-01 | 4.827e-01  | -3.843e+00 | 3.876e+00 | `:g31` |
| TokamakSmallCylindrical  | cylindrical | **0**      | -2.499e-02 | -1.050e+00 | **0**      | -1.051e+00 | 2.502e-02  | 1.001e+00 | `:g12` |
| TokamakMediumCylindrical | cylindrical | **0**      | -1.245e-01 | -2.483e+00 | **0**      | -2.596e+01 | 1.302e+00  | 1.046e+01 | `:g12` |
| TokamakIterCylindrical   | cylindrical | **0**      | 3.888e-01  | -2.303e+00 | **0**      | -8.281e+01 | -1.398e+01 | 3.595e+01 | `:g12` |
| SolovevIter              | cylindrical | **0**      | 2.387e-02  | -2.500e+00 | **0**      | -5.093e+02 | -4.862e+00 | 2.037e+02 | `:g12` |
| SolovevIterXpoint        | cylindrical | 5.926e-03  | 2.626e-02  | -2.500e+00 | 1.207e+00  | -5.093e+02 | -5.349e+00 | 2.037e+02 | `:g31` |
| SolovevSymmetricField    | cylindrical | **0**      | 1.000e+00  | **0**      | **0**      | **0**      | -2.278e+02 | 2.278e+02 | `:g23` |
| TokamakSmallToroidal     | toroidal    | **0**      | -1.498e-03 | -1.054e+00 | **0**      | -5.780e-02 | 8.213e-05  | 5.482e-02 | `:g12` |

**Seven of eleven equilibria have `b₁ = 0` exactly at the initial condition the package ships**, so
`λₒ` vanishes for `(g³, g¹)` and its multipliers are infinite. Those problems cannot be started at
all — the Newton solver hits a NaN in its direction vector on the first step — which is why the
package used to integrate on four equilibria of eleven. `SolovevSymmetricField` is the case that
made all three pairs necessary rather than two: `b = e₂` there, so `b₁` and `b₃` vanish together and
only `(g², g³)` survives.

Two things about this are worth recording:

* **It is not a curvilinear-coordinates problem.** `TokamakSmallCartesian` fails for `(g³, g¹)` too:
  its initial condition sits at `y = z = 0`, where `B_x = 0`. What decides the outcome is where the
  initial condition lies relative to the field, not the coordinate system.
* **`D` is finite everywhere**, including where two of the three pairs are singular. The three
  ratios `{cᵢ,cⱼ}(m) / bₘ` it averages agree to machine precision wherever they are all defined, in
  every chart — which is what makes the average a faithful stand-in for the singular ratio rather
  than a fudge, and is the numerical check on the derivation in
  `src/guiding_center_3d/guiding_center_3d_compact.jl`.


### The formulations agree where they are all usable

A hundred steps at min(`DEFAULT_TIMESTEP`, 0.1) with `PartitionedGauss(2)`. `Δy` is the relative
deviation of the final `(q, p)` from the compact form in its `:parallel` variant, which is used as
the reference because it is the only one of the six defined for every equilibrium.

| equilibrium | variant | Δy | \|ΔH/H\| | max\|gᵏ\| |
|---|---|---|---|---|
| TokamakMediumCartesian | hode `:g31`         | 2.243e-10 | 6.016e-11 | 8.239e-10 |
|                        | hode `:g12`         | 1.930e-09 | 6.028e-11 | 9.677e-09 |
|                        | hode `:g23`         | 5.666e-09 | 6.018e-11 | 2.638e-08 |
|                        | canonical `:g31`    | 2.248e-09 | 2.198e-09 | 1.318e-08 |
|                        | compact `:g31`      | 3.414e-11 | 5.807e-11 | 1.321e-10 |
|                        | compact `:parallel` | —         | 5.745e-11 | 1.345e-10 |

`TokamakMediumCartesian` is the only equilibrium for which no component of `b` vanishes at the
initial condition, so it is the only one where all three pairs can be compared directly. They agree
to 6e-9 — the accuracy of the nonlinear solve at this tolerance — which is the substantive check
that the three pairs parameterise the same constrained system and that the compact form is equivalent
to the Hamilton-Dirac one it was derived from. The same holds on every other equilibrium for the
pairs that are regular there; the full table is in the script's output.

Two entries in that table are not agreements. `Dipole3d` disagrees by 2e-4 between its
well-conditioned pair and its two badly conditioned ones, which is the subject of the next section.
And `SolovevSymmetricField` cannot run `hodeproblem_canonical` at all: its pair is well conditioned,
`λₒ ≈ -228` and never below -173 along the orbit, but the `∂λ/∂q`, `∂λ/∂p` terms the canonicalised
form adds amplify the constraint drift enough that Newton meets a NaN after some thirty steps. The
same orbit runs to completion under `hodeproblem` and `hodeproblem_compact`. It is the only
equilibrium where this happens, and it is a property of that formulation rather than of the
constraint pair.


### Regular is not the same as well conditioned

`Dipole3d` is where that distinction bites. All three pairs are regular at its initial condition —
`λₒ` is -34, 34 and 68 — but the orbit takes `b₁` and `b₂` through zero:

| variant | max\|gᵏ\| at t = 3 | \|ΔH/H\| at t = 3 | \|ΔH/H\| at t = 30 |
|---|---|---|---|
| hode `:g31` | 1.30e-03 | 4.06e-05 | 2e+03 |
| hode `:g23` | 1.69e-03 | 3.98e-05 | not measured |
| hode `:g12` | 1.83e-08 | 1.98e-09 | 1.98e-09 |
| canonical `:g31` | 3.27e-05 | 8.16e-07 | 7.6e+05 |
| canonical `:g12` | 5.63e-06 | 1.59e-08 | 1.66e-08 |

Four orders of magnitude at `t = 3`, and past `t = 30` the badly conditioned pair loses the orbit
altogether under both formulations — `|ΔH/H|` of 2e+03 and 7.6e+05 are not error estimates, they are
destroyed trajectories. With `:g12` both the Hamilton-Dirac and the canonicalised form hold their
levels out to `t = 300`.

`Dipole3d` is therefore the one equilibrium whose `default_constraints` is chosen on the conditioning
of the pair *along the orbit* rather than at the initial condition. The others were left where
regularity put them; see `TODO.md`.

The related effect is that the constraint the pair omits is conserved `1/bₘ` times worse than the two
it retains, since `b₁g² - b₂g¹ + b₃g³ = 0` determines it from them pointwise. Over a thousand steps
of `TokamakMediumCartesian` with `:g31` the retained pair holds to 2e-9 while the omitted `g²`
reaches 1.4e-6, the orbit having passed through `b₁ = 3.5e-4`. That is why `compute_constraints`
reports all three.


### What the formulations cost

Per step, taking the Hamilton-Dirac form as unity, across all eleven equilibria that ship initial
conditions. Three states, because two separate changes bear on this and they have to be told apart:
**A** is the code before the right-hand side was rewritten, on `ElectromagneticFields` 0.6.2; **B** is
after the rewrite, still on 0.6.2; **C** is after it, on 0.6.3, which added common-subexpression
elimination to the generated field code. **C is what the package does now.**

| formulation | A | B | C | why |
|---|---|---|---|---|
| compact, literal pair | 0.79 – 0.93 | 1.01 – 1.12 | 1.02 – 1.10 | replaces both multipliers by a single bracket ratio; no second derivatives |
| Hamilton-Dirac | 1 | 1 | 1 | two multipliers over six bracket terms |
| compact, `:parallel` | 1.19 – 1.44 | 1.03 – 1.34 | 1.11 – 1.37 | the regularised `D` averages all three brackets, so nine terms rather than three |
| canonicalised | 9.3 – 14.8 | 5.7 – 11.2 | 4.8 – 7.2 | `∂λ/∂q` and `∂λ/∂p` need every second derivative of the field |

The single most useful thing in the table is that A → B moved the compact row across one:
**the compact form is no longer cheaper than Hamilton-Dirac.** It won before by evaluating three
bracket terms where the Hamilton-Dirac form evaluates six, and each of those terms re-entered
`dbᵢdxⱼ` from the top. Now the field values are read once per call and shared, so the bracket-term
count barely registers and what shows instead is the extra work the compact form does — the parallel
velocity, and the contraction over `∂A/∂x + u ∂b/∂x` in its momentum equation. Its two variants
converged on each other for the same reason. The upstream release then moved almost nothing here,
which is the expected shape: it makes every field call cheaper by roughly the same factor, and these
are ratios.

Where it does show is the canonicalised row, which nearly halved again from B to C. That row is the
one dominated by second derivatives, and those gain more from the elimination than first derivatives
do — `d²b₁dx₁dx₁` went from 1733 ns to 80 against `db₁dx₁`'s 549 to 57.

A ratio is a statement about two right-hand sides only when both variants ask the solver for the same
work, so rows where they do not are excluded from the ranges above and belong here instead:

* `Dipole3d` needs 1.9 Newton iterations per stage for the Hamilton-Dirac form against 1.0 for
  `compact :parallel`, which is the entire reason that ratio comes out at **0.72** rather than above
  one. Its canonicalised row is 9.7 at 3.0 iterations per stage. That figure read 29.6 before the sign
  error in `∂λ₂/∂q` and `∂λ₂/∂p` was fixed, which had Newton at 3.8; the right-hand side itself costs
  the same either way.
* `SolovevSymmetricField` has no canonicalised row at all — it is the equilibrium where that
  formulation diverges with a `NaN` in the Newton direction, which `TODO.md` still records as open.

Every other row runs at exactly 1.0 iterations in every column, so the rest of the table is a clean
comparison of expressions.

The three pairs cost the same as each other to within noise — they are the same expressions with
permuted indices — so the pair should be chosen on conditioning alone. The canonicalised form remains
the one that matters: it buys exact symplecticity rather than approximate, and gives up as much as two
or three orders of magnitude of energy conservation and four or five of constraint conservation for it
at these step sizes. The spread is wide — on five of the eleven equilibria it is within a factor of two
of the Hamilton-Dirac form on both, and on `TokamakSmallCartesian` it conserves energy slightly better
— so the figure to expect is the worst case, not the typical one.

Separately, the rewrite of the constraints made the Hamilton-Dirac right-hand side **2.7× faster**:
the old code evaluated `λ₁` and `λ₂` once per component, so the multipliers and the six bracket
terms behind them were computed three times per call. A hundred steps of `Dipole3d` went from
1383 ms to 517 ms, and `hodeproblem_canonical` from 11.9 s to 10.1 s, with allocations and compile
time unchanged.


## The toroidal momentum of the guiding centre models

`scripts/study_toroidal_momentum.jl`

Fifteen of the sixteen equilibria that define `toroidal_momentum` defined it as `R(t,q) * ϑ₃(t,q)`.
Relative variation of each candidate over the orbit, Gauss(2):

| equilibrium | coords | `ϑ₃` | `R ϑ₃` | `x ϑ₂ - y ϑ₁` |
|---|---|---|---|---|
| TokamakSmallCylindrical | cylindrical | **1.80e-13** | 3.28e-03 | 6.51e-02 |
| TokamakSmallToroidal    | toroidal    | **6.21e-15** | 3.28e-03 | 1.39e+01 |
| SolovevIterXpoint       | cylindrical | **5.75e-13** | 1.22e-02 | 1.34e+00 |
| TokamakSmallCartesian   | cartesian   | 9.00e-02 | 9.40e-02 | **1.76e-11** |
| TokamakMediumCartesian  | cartesian   | 3.92e+00 | 2.25e+00 | **1.15e-05** |

Two separate conclusions:

* **The factor of `R` is spurious.** `ϑ₃` is already the covariant φ-component of the one-form and
  is the canonical momentum on its own; multiplying by `R` turns a quantity conserved to 1e-13 into
  one that drifts by 1e-3, ten orders of magnitude worse.
* **Cartesian coordinates need a different expression entirely.** There the third coordinate is `z`,
  not an angle, and neither candidate is conserved to any useful degree. The conserved quantity is
  the generator of rotation about the z-axis, `x ϑ₂ - y ϑ₁`.

The residual in the cartesian column is discretisation error, not a defect of the quantity — it
converges at fourth order, as befits Gauss(2):

| Δt | rel. var. | ratio |
|---|---|---|
| 1.0   | 1.15e-05 | |
| 0.5   | 7.22e-07 | 16.0 |
| 0.25  | 4.51e-08 | 16.0 |
| 0.125 | 2.82e-09 | 16.0 |

The conservation is now asserted in `test/structure_tests.jl`, but the comparison between
candidates is not, and it is what identifies the right expression should this ever be revisited for
a new coordinate system.


## The gyrokinetic model is the guiding centre model in rescaled time

`scripts/study_gyrokinetic_rescaling.jl`

`GyroKinetics4d` and `GuidingCenter4d` implement the same physics — same equilibrium Hamiltonian,
same generalised vector potential — but the gyrokinetic module clears the phasespace Jacobian
`B*∥` from the denominator of the characteristics, which is what makes its right-hand side
divergence-free and its splitting volume preserving.

The two vector fields are therefore exactly proportional:

| q = (R₁, R₂, R₃, u) | B*∥ | max rel. diff of `v_gk - B*∥ v_gc` |
|---|---|---|
| [6.2, 0.3, 0.0, 0.34]  | 4.354e+04 | 6.64e-16 |
| [5.5, -0.8, 1.1, -0.2] | 2.113e+04 | 1.21e-16 |
| [7.0, 0.5, 2.0, 0.5]   | 9.258e+04 | 1.40e-16 |

and the six subsystems of the splitting sum to the full field exactly (difference 0.00e+00, not
merely small).

**The direction of the reparametrisation is `dt = B*∥ ds`, not the reverse.** `s` is the "slow"
variable: one unit of `s` covers `B*∥ ≈ 2e2` units of physical time for the shipped ITER-like
equilibrium. Integrating the gyrokinetic model over `s` and the guiding centre model over
`t = B*∥(q₀) s` lands in the same place:

| s interval | endpoint difference | ratio |
|---|---|---|
| 4e-03 | 1.73e-08 | |
| 2e-03 | 4.37e-09 | 4.0 |
| 1e-03 | 1.10e-09 | 4.0 |
| 5e-04 | 2.75e-10 | 4.0 |

The residual is second order in the interval and independent of the step, which identifies it as
the error of the *linear* time map rather than of the integrator: the exact relation is
`t = ∫₀ˢ B*∥(q(s')) ds'`, so freezing `B*∥` at `q₀` costs `O(s²)`. Recovering physical time along an
orbit requires integrating `dt/ds = B*∥` alongside.

This is the check that would have caught the bug that was there: the module used to substitute a
coordinate transformation `q̃ = ω₀ q` into the field without its Jacobian, which breaks the
proportionality above.


## Volume preservation of the gyrokinetic splitting

`scripts/study_volume_preservation.jl`

The point of the six-way splitting is that each subsystem freezes two of the four variables and is
symplectic in the remaining two, so a symplectic method applied to each preserves phasespace
volume *exactly*, at any step size — not to the order of the method. That is a claim about the
discrete map, so it is worth measuring rather than assuming.

Determinant of the Jacobian of the one-step map, by central differences (floor ≈ 1e-11):

| Δs | Strang split | ExplicitEuler | Euler / Δs |
|---|---|---|---|
| 1e-03 | 9.150e-12 | 3.152e-08 | 3.152e-05 |
| 3e-03 | 1.458e-12 | 2.837e-07 | 9.456e-05 |
| 1e-02 | 1.654e-11 | 3.156e-06 | 3.156e-04 |
| 3e-02 | 1.836e-13 | 2.851e-05 | 9.503e-04 |
| 1e-01 | 9.070e-12 | 3.208e-04 | 3.208e-03 |

The splitting sits at the measurement floor across two decades of step size, with no trend — its
volume error is zero to the precision available. Explicit Euler, consistent with the same
divergence-free field, has a volume error growing linearly in `Δs`; the last column is constant to
within a few percent, confirming the expected `O(Δs)`.

A tighter measurement would need automatic differentiation of the one-step map rather than central
differences; `ForwardDiff` is not currently a dependency. Note also that a fourth-order method on
the full field is *not* distinguishable from the splitting this way at these step sizes — its
volume error is `O(Δs⁵)`, below the floor — which is why the control is first order.


## Conservation across the model families

`scripts/study_invariant_conservation.jl`

A side-by-side of what each model conserves, which is the end-to-end check that the diagnostics and
the equations of motion agree.

4D guiding centre, `ODEProblem` with Gauss(2):

| equilibrium | \|ΔH/H\| | \|Δp_φ/p_φ\| |
|---|---|---|
| TokamakSmallCylindrical | 1.43e-15 | 1.63e-13 |
| TokamakSmallToroidal    | 4.77e-16 | 3.34e-15 |
| TokamakSmallCartesian   | 1.39e-13 | 1.76e-11 |
| SolovevIterXpoint       | 1.39e-15 | 4.69e-13 |

3D guiding centre, `HODEProblem` with PartitionedGauss(2). The informative column is not the energy
but the constraints: they vanish along the exact flow, so their magnitude is the drift off the
manifold on which `p = ϑ` — the quantity the approximately symplectic methods are designed to keep
bounded, and the reason `compute_constraints` was added to the diagnostics. The constraint the pair
omits is listed separately from the two it retains, because it is the retained pair's drift amplified
by `1/bₘ`.

| equilibrium | pair | \|ΔH/H\| | \|gᵏ\| retained | \|gᵏ\| omitted | \|Δp_φ/p_φ\| |
|---|---|---|---|---|---|
| SolovevIterXpoint        | `:g31` | 5.39e-15 | 1.49e-13 | 9.33e-13 | 0.00e+00 |
| TokamakMediumCartesian   | `:g31` | 3.23e-10 | 1.78e-09 | 1.45e-06 | — |
| TokamakSmallCylindrical  | `:g12` | 6.84e-15 | 1.31e-15 | 1.11e-18 | 0.00e+00 |
| TokamakMediumCylindrical | `:g12` | 2.04e-13 | 7.43e-13 | 2.62e-14 | 0.00e+00 |
| TokamakSmallToroidal     | `:g12` | 1.12e-15 | 1.22e-17 | 1.50e-20 | 0.00e+00 |

The three equilibria at the bottom are ones that could not be started at all before the constraint
pair became selectable. `SolovevSymmetricField`, which also could not, still cannot be run over a
thousand steps: it is the stiffest equilibrium in the package at `‖ϑ‖ ≈ 256` and the drift takes
Newton into a NaN a few hundred steps in. The `:g12` rows show the amplification working the *other*
way — `b₃` is the large component there, so the omitted constraint comes out better conserved than
the retained pair rather than worse. `TokamakMediumCartesian` is the case where it hurts: `b₁` falls
to 3.5e-4 along that orbit.

`p_φ` is conserved *exactly* for the cylindrical and toroidal cases because φ is cyclic and the
partitioned symplectic method conserves the conjugate momentum of a cyclic coordinate to round-off.

Gyrokinetic guiding centre, over the module's default span (~200 units of physical time):

| method | \|ΔH/H\| |
|---|---|
| Gauss(2) on the full field | 1.21e-07 |
| Strang composition of the splitting | 7.66e-04 |

**The splitting is markedly worse on energy than a symplectic method of comparable cost.** That is
not a bug — it is volume preserving, not symplectic, so nothing constrains its energy error — but
it is four orders of magnitude, and it means the trade the scheme makes is a real one rather than
free. Over the very short orbits used in the tests both sit near 1e-14, so this only shows up at
realistic integration lengths.

Whether the energy error is *bounded* or *drifts* over long runs is the question this does not
answer, and is the study most worth doing next. See `TODO.md` in the repository root.


## The residual scale of the ITER-size equilibria

`scripts/study_solver_tolerances.jl`

An *absolute* solver tolerance has to sit above the round-off floor of the residual it is applied
to, and that floor is a property of the equilibrium. The variational residuals carry the momentum
equation `p = ϑ(q)`, so the floor is `‖ϑ‖ eps`, evaluated here at each equilibrium's shipped
`initial_conditions_barely_passing`:

| family | equilibrium | ‖ϑ‖ or ‖p‖ | `‖ϑ‖ eps` | `f_abstol = 1E-15` |
|---|---|---|---|---|
| GC4d    | SolovevIterXpoint      | 14.95   | 3.3e-15 | unreachable |
| GC4d    | SolovevIter            | 14.95   | 3.3e-15 | unreachable |
| GC4d    | SolovevSymmetricField  | 256.3   | 5.7e-14 | unreachable |
| GC4d    | TokamakIterCylindrical | 29.07   | 6.5e-15 | unreachable |
| GC4d    | TokamakMediumCartesian | 1.17    | 2.6e-16 | reachable |
| GC4d    | TokamakSmallCartesian  | 0.0244  | 5.4e-18 | reachable |
| Pauli3d | SolovevIter            | 0.137   | 3.0e-17 | reachable |

`‖A‖ ~ B₀R₀²`, and ITER has `B₀ = 5.3`, `R₀ = 6.2`, which is why exactly the ITER-size 4D guiding
centre equilibria are the unreachable ones. `SolovevSymmetricField` is the worst of them at
`‖ϑ‖ ≈ 256`.

An earlier version of this page added that `SolovevSymmetricField` "never showed up in CI only
because its block runs `ode` rather than `iode`, and `Gauss` converges on these problems in one or
two iterations regardless of tolerance". **That is wrong**, and measurement says so: at the
`f_abstol = 1E-12` the suite uses, `Gauss(2)` on its `odeproblem` averages 29.4 Newton iterations
and reaches the iteration cap on 45 of 100 steps at `Δt = 1E2`, and on 4 % of steps at `Δt = 1E1`.
It goes quiet only at `Δt = 1E0`, which is ten thousand steps for the same span. The capping is a
step-size limit of the equilibrium, not a tolerance effect, and the solves that cap do not converge
at any iteration budget.

This was not caused by the move to SimpleSolvers 0.10 — the same block capped 28 of 100 steps under
GeometricIntegrators 0.16.10 / SimpleSolvers 0.9.2, at 0.524 ms/step against 0.346 now. 0.10 reports
more of the failures and is 1.5x faster per step; what it added is the diagnosis, namely that
`φ'(0)` is inconsistent with the merit the line search measures.

The module declares `DEFAULT_TIMESTEP = 1E0` over `DEFAULT_TIMESPAN = (0.0, 1E3)`, and every test
runs at it. That rate is not a convenient choice for `Gauss(2)` but the only stable one measured:
over the same physical span, the maximum relative energy error is

| Δt | `odeproblem`, `Gauss(2)` | `iodeproblem`, `VPRKGauss(2)` |
|---|---|---|
| 5E3 | 9.2e-01 | 2.0e-01 |
| 1E3 | 1.0e+00 | 8.7e-01 |
| 1E2 | 1.8e+01 | 2.3e+03 |
| 1E1 | 1.8e+11 | 1.5e+18 |
| 1E0 | **1.9e-03** | 6.2e+25 |

### The variational formulation is unstable here, and finer steps make it worse

The right-hand column is the finding. The `odeproblem` becomes well behaved at `Δt = 1E0`; the
`iodeproblem` becomes catastrophically worse, `|R|` reaching 1.3e+5 against a machine radius of 2.5.
The degenerate variational discretisation carries a parasitic mode that plain VPRK does not control
and that grows as the step shrinks — the reverse of the usual convergence expectation, and the
reason the old `Δt = 5E3` block looked healthy: it took two steps.

The solver's iteration cap is a *symptom* of this, not a cause. At `Δt = 1E0` it is reached on 97 %
of steps, but it is being asked to solve for a state that has already left the device, and no solver
setting repairs that. Measured: `DogLeg` throws on a NaN direction for all four orbits, `StrongWolfe`
throws a `SingularException` on one and caps just as hard on the rest, `Static` throws on two, and
`PartitionedGauss(2)` is singular.

**`SymmetricProjection(VPRKGauss(2))` is what fixes it**, because it constrains exactly the drift off
`p = ϑ(q)` that the parasitic mode feeds on. All four orbits then integrate silently and stay bounded
at `|R| ≤ 6.5`, the same excursion `Gauss(2)` gives. `MidpointProjection` is not sufficient — it holds
`barely_trapped` and throws on a NaN direction for the other three.

Since the parasitic mode is a property of the degenerate formulation rather than of this
equilibrium, **every variational orbit in the 4D guiding centre tests is integrated with the
projected method**. The one exception is the `ThetaPinchField` loop, where `p` is an exact invariant,
so there is no drift to project and the projection's inner initial guess would re-introduce a
degenerate `p` history.

This restores *stability*, not accuracy: on `SolovevSymmetricField` the projected energy error is
O(1) against the `odeproblem`'s 1.9e-3 over the same span. On this equilibrium the variational
formulation should not be read as a quantitative result, and the tests assert only that the
integration completes. Elsewhere the projected variational orbits track the `ode` ones closely —
4.7e-6 against 3.1e-6 on `TokamakMediumCartesian`, and to three digits on `SolovevIterXpoint`.

### Larger time steps do not help, and the apparent evidence that they do is an artifact

Read down the `Δt` column of the first table and the error appears to *fall* as the step grows past
`5E3`, which would suggest running these examples at `Δt = 1E4` or beyond. It does not survive a
controlled comparison. That table holds the *span* fixed, so a larger step means proportionally fewer
steps and proportionally less opportunity for the instability to accumulate — at `Δt = 1E5` over
`t ∈ [0, 1E6]` the run is ten steps, which cannot resolve the orbit at all.

Holding the *step count* fixed at 1000 instead, and letting the span follow, gives the ordinary
picture — monotone in `Δt`, in the direction one expects:

| Δt | span | `odeproblem` max\|ΔH/H\| | `iodeproblem` + projection |
|---|---|---|---|
| 1E-1 | 1E2 | 9.4e-08 | 5.7e-07 |
| 1E0  | 1E3 | **1.9e-03** | 9.6e-02 |
| 1E1  | 1E4 | 1.8e+11 | throws |
| 1E2  | 1E5 | 4.9e+22 | throws |
| 1E3  | 1E6 | 9.2e+24 | 8.3e+201 |
| 1E4  | 1E7 | 1.3e+32 | 2.6e+212 |

`Δt = 1E4` is therefore the worst configuration measured, not a viable one. `Δt = 1E0` is the largest
step at which this equilibrium is resolved; `Δt = 1E-1` buys four more orders of accuracy for ten
times the work, which is the trade available if it is ever wanted.

With this and one mis-set time step in the Pauli tests corrected, the suite emits no
nonlinear-solver warnings at all, and asserts as much.

Two knobs are involved, and they fail for different models. Cost per step at the ITER Solov'ev
X-point, with the Newton iteration count behind it:

| setting | GC4d `iode` ms/step | mean iters | Pauli3d `iode` ms/step | mean iters |
|---|---|---|---|---|
| `f_abstol=1E-15, f_reltol=1E-15` | 0.130 | 3.8 | 0.044 | 3.3 |
| `f_abstol=1E-15`                 | 0.126 | 3.7 | 0.039 | 2.9 |
| `f_abstol=1E-14`                 | 0.046 | 1.0 | 0.031 | 2.0 |
| `f_abstol=1E-12`                 | 0.045 | 1.0 | 0.031 | 2.0 |

The 4D guiding centre degrades from `f_abstol = 1E-15` alone and recovers between 1E-15 and 1E-14,
exactly where its floor lies. The Pauli model has a small momentum and is untouched by `f_abstol`;
it degrades only once `f_reltol` is pinned as well, i.e. what hurts it is the loss of the relative
term `f_reltol ‖F(x₀)‖` that lets a large-magnitude solve converge at all. `f_abstol = 1E-12` with
`f_reltol` left at its default is the only setting comfortable for both, and is what the test suite
uses. **The final state is unchanged to ~1e-14 across every row** — the extra iterations buy
nothing.

The unreachable rows are cheap because of `max_stalls`: a solve whose iterate stops moving gives up
after two such steps rather than running to `max_iterations`. What that is worth is visible by
comparing the same workloads under GeometricIntegrators 0.16.10 / SimpleSolvers 0.9.2, which lacks
it, with 0.17 / 0.10.1:

| setting | GC4d `iode` 0.9.2 | 0.10.1 | Pauli3d `iode` 0.9.2 | 0.10.1 |
|---|---|---|---|---|
| `f_abstol=1E-15, f_reltol=1E-15` | 30.8 ms/step, 246.8 iters | 0.130, 3.8 | 6.32 ms/step, 127.8 iters | 0.044, 3.3 |
| `f_abstol=1E-15`                 | 21.7 ms/step, 172.2 iters | 0.126, 3.7 | 0.272 ms/step, 2.9 iters | 0.039, 2.9 |
| `f_abstol=1E-14`                 | 0.047, 1.0                | 0.046, 1.0 | 0.058, 2.0               | 0.031, 2.0 |
| `f_abstol=1E-12`                 | 0.047, 1.0                | 0.045, 1.0 | 0.030, 2.0               | 0.031, 2.0 |

An unsatisfiable tolerance therefore costs about four Newton iterations rather than 170, while the
rows that *do* converge are unchanged to within run-to-run noise. This does not make an `f_abstol`
below the residual floor correct — the solver still cannot deliver what was asked, and says so — but
it removes the pathological cost of asking, which is what turns such a mistake into a ten-minute CI
problem rather than a warning.

Note this also rules out simply dropping the options. `GeometricIntegratorsBase` 0.5.1 scales its
default with the size of the stage system,
`f_abstol = max(8, solversize(method, problem)) * eps(datatype(problem))`, which for the workloads
here is `1.8e-15` (`solversize = 8`, the 4D guiding centre and Pauli variational solves) or
`2.7e-15` (`solversize = 12`, the 3D `hodeproblem` under `PartitionedGauss(2)`). Both are still
below the ITER floor, so the sized default does not rescue these problems and the explicit
`f_abstol` stays.

Once the tolerance is sound the remaining cost has nothing to do with ITER. Two changes bear on it —
this package's right-hand side rewrite, and `ElectromagneticFields` 0.6.3 — so all four combinations
were measured, in interleaved randomised rounds on one machine, minimum per row. **A** is the code
before the rewrite on 0.6.2, **B** after it on 0.6.2, **C** after it on 0.6.3, **D** before it on
0.6.3. **C is what the package does now**; **D** exists so that the two factors can be separated
rather than asserted.

| workload | A | B | C | D | A→C |
|---|---|---|---|---|---|
| GC3d SolovevIterXpoint `hodeproblem` + extrapolation      |  5.87 | 0.875 | 0.302 | 1.758 | 19.5× |
| GC3d SolovevIterXpoint `hodeproblem`                      |  1.59 | 0.230 | 0.082 | 0.579 | 19.5× |
| GC3d SolovevIterXpoint `hodeproblem_compact`              |  2.22 | 0.242 | 0.095 | 0.804 | 23.3× |
| GC3d SolovevIterXpoint `hodeproblem_canonical`            | 19.94 | 1.566 | 0.389 | 5.999 | 51.2× |
| GC3d TokamakMediumCartesian `hodeproblem` + extrapolation |  8.09 | 1.454 | 0.375 | 2.324 | 21.6× |
| GC3d TokamakMediumCartesian `hodeproblem`                 |  2.21 | 0.389 | 0.104 | 0.701 | 21.4× |
| GC3d TokamakMediumCartesian `hodeproblem_canonical`       | 32.33 | 4.348 | 0.594 | 7.887 | 54.4× |
| GC4d SolovevIterXpoint `ode`                              | 0.094 | 0.094 | 0.036 | 0.036 |  2.6× |
| GC4d SolovevIterXpoint `iode`                             | 0.120 | 0.120 | 0.042 | 0.042 |  2.8× |
| Pauli3d SolovevIter `hodeproblem`                         | 0.065 | 0.065 | 0.043 | 0.043 |  1.5× |
| Pauli3d SolovevIter `iode`                                | 0.041 | 0.042 | 0.029 | 0.029 |  1.4× |

ms/step, Newton at 1.0 iterations on every 3D row throughout. The last four rows are the control: the
4D guiding centre and Pauli source trees are byte-identical between the two commits, so A = B and
C = D there, as they come out — but they use the same generated field functions, so they do gain from
0.6.3.

**For `hodeproblem` the two changes are almost exactly multiplicative.** The rewrite is 6.9× measured
on 0.6.2 (A/B) and 7.1× on 0.6.3 (D/C); the release is 2.7× on the old code (A/D) and 2.8× on the new
(B/C). 6.9 × 2.8 = 19.4 against 19.5 measured, so quoting the two factors separately is legitimate
here.

**For `hodeproblem_canonical` they are not**, and the interaction is the interesting part: the rewrite
is 12.7× on 0.6.2 but 15.4× on 0.6.3, and the release is 3.3× on the old code but 4.0× on the new.
That row is the one dominated by second derivatives, and the elimination helps those more than first
derivatives — `d²b₁dx₁dx₁` from 1733 ns to 80 against `db₁dx₁`'s 549 to 57 — so the two changes
reinforce each other rather than merely stacking.

What the test suite ran before is the first row and what it runs now is the second, so the figure it
sees is **72×**. The gap to the Pauli `hodeproblem` that prompted the whole comparison closes from
**24× to 1.9×**.

Two caveats on the absolute figures. They are lower throughout than the ones this section carried
before, because the harness now repeats each integration inside one timed region until it exceeds
50 ms — a single 100-step integration of the Pauli model is 8 ms, where operating-system noise
dominates — and because it uses `@timed` with a per-sample assertion that no compilation entered the
sample. The ratios are what transfer between harnesses, not the milliseconds. And repeatedly
integrating the same hundred steps keeps 0.6.2's very much larger generated code hot in the
instruction cache in a way a cold workload would not, so if anything these figures *understate* what
0.6.3 is worth.

`TokamakMediumCartesian` is *slower* than the ITER Solov'ev, in both columns and for a reason that
has nothing to do with the solver: in a cartesian chart the metric is the identity, so `dgⁱʲdxₖ` and
its second derivatives fold away to literal zeros — but the *field* then carries the whole
`sqrt(x²+y²)` toroidal geometry into every component of `db`, which is where the cost is. In the
curvilinear charts that geometry sits in the metric instead, which is diagonal and cheap.

The 3D `hodeproblem` path cost two orders of magnitude more per step than the Pauli one even though
Newton converges in a single iteration. That was for a long time attributed to "nested `ForwardDiff`: the second derivatives of the 3D
Hamiltonian differentiate field functions that are themselves AD-generated, and the backtracking
line search re-evaluates them". **Every part of that is wrong for this path**, and profiling says so:

* There is no nesting, and no AD in the field layer at all. `ElectromagneticFields` generates its
  field functions **symbolically, with SymEngine, at precompile time** — its dependencies are
  `Combinatorics, Documenter, LaTeXStrings, LinearAlgebra, NaNMath, RecipesBase, SymEngine`, and
  `ForwardDiff` is not among them. The only automatic differentiation anywhere in the stack is the
  Newton Jacobian, which `SimpleSolvers` takes of the residual.
* `hodeproblem` never touches a second derivative of the Hamiltonian. Its right-hand side runs
  `λ₁`/`λ₂` → `bracket_gg`/`bracket_gH` → `dgᵏd…`, `dHdq`, `dHdp`, all first derivatives. The Hessian
  in `guiding_center_3d_canonical.jl` belongs to `hodeproblem_canonical` alone, and is hand-written
  closed form rather than differentiated at runtime.
* The line search is 8 % of wall clock and the Jacobian 13 %. Per step the right-hand side is
  evaluated **93.8 times on `Float64` against 4.0 times on `Dual`** — automatic differentiation sees
  four per cent of the calls, so it cannot be the cost whatever it is doing.

Only the last clause survived. `initial_guess!` was **70 %** of wall clock, and
`MidpointExtrapolation(5)` is not the default for either method — `PartitionedGauss` and `VPRKGauss`
both declare `default_iguess() = HermiteExtrapolation()`. Nothing recorded why the 3D tests asked for
it; the one documented reason is the Pauli theta pinch, whose momentum is constant along the initial
trajectory.

What the cost actually was, in three layers. All three are now addressed; the first two here, the
third upstream.

1. **The extrapolator**, 3.0–3.5× across all thirteen equilibria and all three formulations, for
   bit-identical trajectories and no change in iteration count.
2. **The right-hand side re-entering the generated field code.** On `ElectromagneticFields` 0.6.2,
   `db₁dx₁` for the ITER Solov'ev X-point was 1905 statements containing 108 separate evaluations of
   `log(x₁)` and cost 549 ns, against a `g¹¹` that folds to a constant. One evaluation of the Hamilton-Dirac
   right-hand side called it more than seventy times: each of the twelve Poisson-bracket terms
   re-entered from the top, and `λ₁` and `λ₂` each recomputed the `λₒ` they share. Evaluating each
   injected function once into a `FieldValues` and reading the results back took
   `guiding_center_3d_v` from 24.7 µs to 3.6 µs — 6.9×, which is exactly the end-to-end per-step
   ratio in the table above, so the whole of that factor is the right-hand side and none of it the
   solver. On 0.6.3 the same call is 770 ns, 32× below where it started.
3. **The expression swell itself**, which lives in `ElectromagneticFields` rather than here.
   `convert(Expr, ::Basic)` wrote out SymEngine's expanded tree verbatim and SymEngine shares
   nothing, so a subexpression common to fifty branches of `∂ⱼ(Bᵢ / sqrt(gᵏˡBₖBₗ))` was emitted fifty
   times. **[0.6.3](https://github.com/JuliaPlasma/ElectromagneticFields.jl/pull/10) eliminates
   common subexpressions in the generated bodies**, which takes `db₁dx₁` from 549 ns to 57 and
   `d²b₁dx₁dx₁` — which `hodeproblem_canonical` needs — from 1733 ns to 80. The pass only *names*
   subexpressions, so the values do not change: every one of the 27456 field values this package
   injects, across all 32 equilibrium modules of the three families, is bit-identical between 0.6.2
   and 0.6.3, and so is every trajectory integrated from them. `Project.toml` requires 0.6.3 for this
   reason, and the tables above are measured against it.
