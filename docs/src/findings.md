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

Every number here depends on which `ElectromagneticFields` is resolved, and in two different ways.
Every *timing* depends on it because that package emits the field functions the right-hand sides are
built from, and 0.6.3 made them several times faster than 0.6.2. Every *value* on six of the eleven
equilibria depends on it because **0.7.0 reversed the magnetic field in the four left-handed charts**;
see the next section. Before quoting a number against a modified environment, check which version is
actually loaded — and do not check it with `Pkg.status`, which for a path dependency prints the
version recorded in the manifest rather than the version at the path:

```
julia --project=scripts -e 'println(Base.locate_package(Base.identify_package("ElectromagneticFields")))'
```

That resolves the load path without loading anything. A `scripts/Manifest.toml` carrying a
`path = "../../ElectromagneticFields"` stanza with no matching `[sources]` entry in
`scripts/Project.toml` is how one round of these measurements was silently taken against a working
checkout on a feature branch; deleting the manifest and re-resolving is the fix.


## The reversed field of the left-handed charts

`ElectromagneticFields` 0.7.0
([PR #11](https://github.com/JuliaPlasma/ElectromagneticFields.jl/pull/11)) fixed an orientation
error: `B¹ = ⋆dA` is orientation-dependent and `hodge²¹` was being handed the *unsigned* volume
element `J = |det DF|`. Four of that package's charts are left-handed, `(R, Z, ϕ)` and `(r, θ, ϕ)`
both having `det DF < 0`, and in every one of them `B` came out **antiparallel** to what the
cartesian chart gives at the same physical point. Six of this package's equilibria are built on
them:

| affected — `b` reversed | unaffected — cartesian charts |
|---|---|
| `TokamakSmallCylindrical`, `TokamakMediumCylindrical`, `TokamakIterCylindrical` | `TokamakSmallCartesian`, `TokamakMediumCartesian` |
| `TokamakSmallToroidal` | `SolovevSymmetricField` |
| `SolovevIter`, `SolovevIterXpoint` | `Dipole3d`, `QuadraticPotentials3d`, `ThetaPinchField`, `SymmetricQuadratic` |

**`SolovevSymmetricField` is the one to watch**: despite the name it is a *cartesian* chart, so it is
right-handed and unaffected, unlike every other `Solovev*` equilibrium. The `coords` column of the
conditioning table below said "cylindrical" for it, which was wrong.

### The models are exactly equivariant, so the fix is a sign on `u`

All three guiding-centre-related families are built from

```
ϑ = A + u b ,    H = ½ u² + μ |B| + φ
```

and `A`, `|B|` and `φ` are untouched by the upstream change — only the unit vector `b` reverses. Both
`ϑ` and `H` are therefore invariant under

```
b → -b   together with   u → -u
```

so **the old dynamics at parallel velocity `u` is exactly the new dynamics at `-u`**. That is not an
approximation: the 4D vector field comes out bit-identical component for component (with `v_u` negated),
and integrating fifty steps of every affected equilibrium under 0.6.3 and under 0.7.0-with-`u`-negated
reproduces the final state to `0.000e+00` or a few units in the last place, in all of `GuidingCenter3d`,
`GuidingCenter4d` and `PauliParticle3d`. The cartesian equilibria are bit-identical with `u` left alone,
which is the control.

`GuidingCenter3d` is invariant for a second reason worth noting: its state is `(q, p)` with `p = ϑ`,
which does not change, and although each constraint `gᵏ = bᵢvⱼ - bⱼvᵢ` flips sign, so does each
multiplier, and the products `λ ∂g/∂p` that enter the right-hand side do not. The sign enters only
through the initial `p = ϑ(q, u)`.

### What actually needed changing

Almost nothing, and that is the interesting part. The `initial_conditions_*` of every module were
already written chart-independently — `u = 8.117E-4` in the cartesian, cylindrical *and* toroidal
charts of the small tokamak — so under 0.6.3 they silently denoted **different physical particles in
different charts**, and under 0.7.0 they finally denote the same one. They are left exactly as they
were.

What did need changing was a hand-compensation. The module-level default of the small tokamak read

```
uᵢ = +0.00045135897235326736   # tokamak_small_cartesian.jl
uᵢ = -0.00045135897235326736   # tokamak_small_cylindrical.jl, tokamak_small_toroidal.jl
```

at the same position with the same `μ` — someone had matched the physical particle across charts by
negating `u` where the field came out backwards. That compensation is now stale and is removed, so all
three charts again start one particle. It is the clearest evidence that the reversal was known about
and worked around rather than unnoticed.

### The orientation shows through anywhere a coordinate curl appears

Two derived quantities in this package are built from `∂₂ϑ₃ - ∂₃ϑ₂`-shaped expressions, which are
`det(DF)` times the contravariant curl rather than `|det DF|` times it. Both therefore carry the
chart's handedness, and both come out **negative in the left-handed charts**:

* `λₒ = {g₁, g₂}` and the compact form's `D`, in the 3D guiding centre model. In the table below
  `λₒ(:g23)` and `b₂` now share a sign in the left-handed charts and have opposite signs in the
  right-handed ones, where before the reversal made it look uniform. Nothing in the dynamics depends
  on it — the multipliers divide by `λₒ` and are multiplied by derivatives carrying the same factor.
* `ωabs` in the gyrokinetic model, which is `det(DF) · B*∥` and not `B*∥`. Dividing the factor back
  out recovers a chart-independent `B*∥`: the small tokamak gives `0.9511` in its cartesian,
  cylindrical and toroidal charts alike, where `ωabs` reads `+0.9511`, `-0.9986` and `-0.0499`. This
  one is visible, because the rescaled time is defined by `dt = ωabs ds`, so **`s` runs opposite to
  physical time in six of the eight gyrokinetic equilibria.** The splitting is still volume
  preserving — the coordinate divergence of a coordinate curl vanishes identically whatever its sign
  — and the proportionality `v_gk = ωabs v_gc` still holds exactly. Until 0.7.0 the two sign errors
  cancelled and `ωabs` came out positive everywhere.


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
pair flips `λₒ` and both multipliers together — so only `|λₒ|` is meaningful. On top of that, `λₒ` is
built from a coordinate curl and so carries the chart's orientation as well, which is why in the table
below `λₒ(:g23)` and `b₂` have opposite signs in the right-handed charts and the *same* sign in the
left-handed ones. Neither factor affects the dynamics; see the previous section.


### Which pairs are usable

At the shipped initial condition of each equilibrium. `λₒ = {g₁, g₂}` is the bracket both
multipliers divide by; `D` is the same quantity as the compact form computes it, as a weighted mean
over all three pairs that never divides by a single component of `b`.

| equilibrium | coords | b₁ | b₂ | b₃ | λₒ(:g31) | λₒ(:g12) | λₒ(:g23) | D | default |
|---|---|---|---|---|---|---|---|---|---|
| Dipole3d                 | cartesian   | -4.082e-01 | -8.165e-01 | 4.082e-01 | -3.402e+01 | 3.402e+01  | 6.804e+01  | 8.333e+01  | `:g12` |
| QuadraticPotentials3d    | cartesian   | 2.000e-03  | -3.000e-03 | 1.000e+00 | 2.000e-01  | 9.999e+01  | 3.000e-01  | 9.999e+01  | `:g31` |
| TokamakSmallCartesian    | cartesian   | **0**      | 9.997e-01  | 2.499e-02 | **0**      | 2.379e-02  | -9.516e-01 | 9.519e-01  | `:g23` |
| TokamakMediumCartesian   | cartesian   | **0**      | 9.923e-01  | 1.240e-01 | **0**      | 4.812e-01  | -3.849e+00 | 3.879e+00  | `:g23` |
| TokamakSmallCylindrical  | cylindrical | **0**      | 2.499e-02  | 1.050e+00 | **0**      | -1.049e+00 | 2.498e-02  | -9.995e-01 | `:g12` |
| TokamakMediumCylindrical | cylindrical | **0**      | 1.240e-01  | 2.481e+00 | **0**      | -2.406e+01 | 1.203e+00  | -9.698e+00 | `:g12` |
| TokamakIterCylindrical   | cylindrical | **0**      | -3.888e-01 | 2.303e+00 | **0**      | -8.149e+01 | -1.375e+01 | -3.538e+01 | `:g12` |
| SolovevIter              | cylindrical | **0**      | -2.387e-02 | 2.500e+00 | **0**      | -5.093e+02 | -4.862e+00 | -2.037e+02 | `:g12` |
| SolovevIterXpoint        | cylindrical | -5.926e-03 | -2.626e-02 | 2.500e+00 | 1.207e+00  | -5.093e+02 | -5.349e+00 | -2.037e+02 | `:g31` |
| SolovevSymmetricField    | cartesian   | **0**      | 1.000e+00  | **0**      | **0**     | **0**      | -2.278e+02 | 2.278e+02  | `:g23` |
| TokamakSmallToroidal     | toroidal    | **0**      | 1.250e-03  | 1.050e+00 | **0**      | -5.246e-02 | 6.245e-05  | -4.997e-02 | `:g12` |

Every `b` in the six left-handed rows is the negative of what this table carried before
`ElectromagneticFields` 0.7.0, and `D` with it; the five cartesian rows are bit-identical. `D` being
negative where the chart is left-handed is not a defect: see "The orientation shows through anywhere a
coordinate curl appears" above.

Three rows also moved for an unrelated reason. `TokamakMediumCylindrical`, `TokamakSmallToroidal` and
`TokamakMediumCartesian` had their initial conditions moved from a cartesian `y = 0.1` to `y = 0` to
match their 4D and Pauli counterparts — in the toroidal chart that changes `r` from 0.0548 to 0.05 and
so the flux surface, and in the cartesian one it puts the condition on the `b₁ = b_x = 0` midplane, so
`TokamakMediumCartesian` moves from `:g31` to `:g23`. That last is the better pair anyway: `b₂ = 0.9923`
is the largest component of `b` there, giving `λₒ = -3.85` against `:g31`'s `-0.154` at `y = 0.1`, and
`hodeproblem` conserves the constraints about three times better on it. What it costs is the three-pair
comparison, which moves to `SolovevIterXpoint`.

**Eight of eleven equilibria have `b₁ = 0` exactly at the initial condition the package ships**, so
`λₒ` vanishes for `(g³, g¹)` and its multipliers are infinite. Those problems cannot be started at
all — the Newton solver hits a NaN in its direction vector on the first step — which is why the
package used to integrate on four equilibria of eleven. `SolovevSymmetricField` is the case that
made all three pairs necessary rather than two: `b = e₂` there, so `b₁` and `b₃` vanish together and
only `(g², g³)` survives.

Only three have all three pairs regular at once: `Dipole3d`, `QuadraticPotentials3d` and
`SolovevIterXpoint`. `TokamakMediumCartesian` was a fourth until its initial condition moved to `y = 0`.

Three things about this are worth recording:

* **It is not a curvilinear-coordinates problem.** Both cartesian tokamaks fail for `(g³, g¹)`: their
  initial conditions sit at `y = z = 0`, where `B_x = 0`. What decides the outcome is where the
  initial condition lies relative to the field, not the coordinate system.
* **A cartesian chart is the *harder* case for this pair, not the easier one.** In a cylindrical or
  toroidal chart `b₁ = b_R` vanishes on the midplane but `b₃ = b_φ` is the large component, so `(g¹, g²)`
  is well conditioned there. In a cartesian chart the toroidal direction is split across `b₁` and `b₂`,
  so the pair that survives depends on where on the flux surface the particle sits — which is why both
  cartesian tokamaks default to `(g², g³)` and every curvilinear one to `(g¹, g²)`.
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
| SolovevIterXpoint      | hode `:g31`         | 7.909e-16 | 1.738e-15 | 4.085e-14 |
|                        | hode `:g12`         | 5.877e-17 | 8.690e-16 | 2.927e-14 |
|                        | hode `:g23`         | 1.472e-17 | 8.690e-16 | 2.925e-14 |
|                        | canonical `:g31`    | 5.957e-14 | 9.784e-14 | 4.624e-12 |
|                        | compact `:g31`      | 2.226e-16 | 8.690e-16 | 2.008e-14 |
|                        | compact `:parallel` | —         | 8.690e-16 | 2.927e-14 |

`SolovevIterXpoint` is the sharpest equilibrium for this comparison: no component of `b` vanishes at
its initial condition, so all three pairs are regular there at once, and being curvilinear it is far
less stiff than the cartesian charts. The three agree to **8e-16** — the accuracy of the nonlinear
solve at this tolerance — which is the substantive check that the three pairs parameterise the same
constrained system and that the compact form is equivalent to the Hamilton-Dirac one it was derived
from. The same holds on every other equilibrium for the pairs that are regular there; the full table is
in the script's output.

This table was `TokamakMediumCartesian` until its initial condition moved to `y = 0`, where `b₁ = b_x`
vanishes identically and `:g31` becomes singular. That chart also gave a three-pair agreement of only
6e-9 — seven orders looser than the ITER Solov'ev — so the check is stronger where it now lives, not
weaker. `Dipole3d` and `QuadraticPotentials3d` are the other two equilibria with all three pairs
regular; the dipole is the subject of the next section and the quadratic potentials agree to 4e-15.

Two entries in that table are not agreements. `Dipole3d` disagrees by 2e-4 between its
well-conditioned pair and its two badly conditioned ones, which is the subject of the next section.
And `SolovevSymmetricField` cannot usefully run `hodeproblem_canonical`: its pair is well conditioned,
`λₒ ≈ -228` and never below -173 along the orbit, but the `∂λ/∂q`, `∂λ/∂p` terms the canonicalised
form adds amplify the constraint drift without bound. The run does not stop — it reaches `|ΔH/H|` of
1.5e+228 and `max|gᵏ|` of 2.6e+114 over the hundred steps, at 44.2 Newton iterations per solve against
1.1 for the other formulations and the cap reached on most of them — which is a destroyed trajectory
reported as a number rather than an error, and worth knowing when reading that row. The same orbit runs
to completion under `hodeproblem` and `hodeproblem_compact` at 4.5e-12 and 2.4e-11. It is the only
equilibrium where this happens, and it is a property of that formulation rather than of the constraint
pair. This equilibrium is a cartesian chart and so untouched by the 0.7.0 field reversal.


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

`Dipole3d` is therefore one of two equilibria whose `default_constraints` is chosen on conditioning
rather than on bare regularity, and it is the only one chosen on conditioning *along the orbit*: all
three of its pairs are regular at the initial condition and it is the trajectory that takes `b₁` and
`b₂` through zero. `TokamakMediumCartesian` is the other, and its reason is simpler — `b₂ = 0.9923` is
the largest component of `b` at its initial condition, giving `λₒ = -3.85` where `:g12` gives `0.48`,
and `:g31` is not merely badly conditioned there but singular. The remaining nine were left where
regularity put them; see `TODO.md`.

The related effect is that the constraint the pair omits is conserved `1/bₘ` times worse than the two
it retains, since `b₁g² - b₂g¹ + b₃g³ = 0` determines it from them pointwise. Over a thousand steps of
`TokamakMediumCartesian` with `:g23` the retained pair holds to 2.3e-9 while the omitted `g¹` reaches
9.9e-7 — a factor of 430, `b₁` being the small component that pair divides out. That is why
`compute_constraints` reports all three. It read 2e-9 against 1.4e-6 on the `:g31` pair this module
used to default to, so the effect survives the change of pair rather than being an artefact of it.


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

* `Dipole3d` runs at 2.0–2.1 Newton iterations per stage across its Hamilton-Dirac, compact and
  `:parallel` variants, and 4.1 for the canonicalised one, where every other equilibrium here converges
  in 1.0. Its ratios are therefore comparisons of solver work as much as of expressions:
  `compact :parallel` comes out at 1.13 and the canonicalised form at 12.2. That last figure read 29.6
  before the sign error in `∂λ₂/∂q` and `∂λ₂/∂p` was fixed, which had Newton at 3.8; the right-hand
  side itself costs the same either way.
* `SolovevSymmetricField`'s canonicalised row runs at **44.2** Newton iterations per stage against
  1.1 for its other variants, and produces the destroyed trajectory described above, so its 34x
  cost ratio measures the solver failing rather than the expression. `TODO.md` records it as open.

Every other row runs at exactly 1.0 iterations in every column, so the rest of the table is a clean
comparison of expressions. The 0.7.0 field reversal moves nothing in this table: it changes the sign of
a term in the generated code and not the amount of arithmetic, and the ranges re-measure inside their
stated bounds (1.02–1.10, 1.11–1.37 and 4.8–7.2 excluding the two rows above).

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
| TokamakSmallCylindrical | cylindrical | **1.97e-14** | 4.38e-03 | 8.45e-02 |
| TokamakSmallToroidal    | toroidal    | **8.81e-16** | 4.38e-03 | 1.59e+00 |
| SolovevIterXpoint       | cylindrical | **9.75e-14** | 1.05e-02 | 1.07e+00 |
| TokamakSmallCartesian   | cartesian   | 9.00e-02 | 9.40e-02 | **1.76e-11** |
| TokamakMediumCartesian  | cartesian   | 3.92e+00 | 2.25e+00 | **1.15e-05** |

The three curvilinear rows changed with the 0.7.0 field reversal — they are the affected charts — and
the two cartesian rows are bit-identical. Both conclusions below are unchanged by it, and the margins
they rest on are if anything wider.

Two separate conclusions:

* **The factor of `R` is spurious.** `ϑ₃` is already the covariant φ-component of the one-form and
  is the canonical momentum on its own; multiplying by `R` turns a quantity conserved to 1e-14 into
  one that drifts by 4e-3, eleven orders of magnitude worse.
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

The two vector fields are therefore exactly proportional, with `ωabs` the factor:

| q = (R₁, R₂, R₃, u) | `ωabs` | max rel. diff of `v_gk - ωabs · v_gc` |
|---|---|---|
| [6.2, 0.3, 0.0, 0.34]  | -4.354e+04 | 0.00e+00 |
| [5.5, -0.8, 1.1, -0.2] | -2.113e+04 | 1.21e-16 |
| [7.0, 0.5, 2.0, 0.5]   | -9.258e+04 | 1.18e-16 |

and the six subsystems of the splitting sum to the full field exactly (difference 0.00e+00, not
merely small).

**The factor is negative here, and that is the chart and not a defect.** `ωabs` is `det(DF) · B*∥`,
and the Solov'ev chart is left-handed, so it carries the opposite sign to the physical
`B*∥ > 0`; see "The orientation shows through anywhere a coordinate curl appears" above. Six of the
eight gyrokinetic equilibria are in this position and the two cartesian ones are not. The consequence
is that **`s` runs opposite to physical time in those six**: the relation is `dt = ωabs ds` with
`ωabs < 0`. Everything else survives it — the proportionality above is exact, the splitting is still
volume preserving, and the comparison below still recovers the same orbit — because a coordinate curl
is divergence-free whatever its sign. Until `ElectromagneticFields` 0.7.0 the reversed `b` cancelled
this and `ωabs` came out positive in every chart.

**`s` is the "slow" variable either way.** One unit of `s` covers `|ωabs| ≈ 2e2` units of physical
time at the shipped ITER-like initial condition. Integrating the gyrokinetic model over `s` and the
guiding centre model over `t = ωabs(q₀) s` — a *negative* interval, so backwards — lands in the same
place:

| s interval | endpoint difference | ratio |
|---|---|---|
| 4e-03 | 1.69e-08 | |
| 2e-03 | 4.20e-09 | 4.0 |
| 1e-03 | 1.04e-09 | 4.0 |
| 5e-04 | 2.61e-10 | 4.0 |

The residual is second order in the interval and independent of the step, which identifies it as
the error of the *linear* time map rather than of the integrator: the exact relation is
`t = ∫₀ˢ ωabs(q(s')) ds'`, so freezing `ωabs` at `q₀` costs `O(s²)`. Recovering physical time along an
orbit requires integrating `dt/ds = ωabs` alongside.

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
| 1e-03 | 1.425e-12 | 3.768e-08 | 3.768e-05 |
| 3e-03 | 1.022e-11 | 3.393e-07 | 1.131e-04 |
| 1e-02 | 7.738e-12 | 3.774e-06 | 3.774e-04 |
| 3e-02 | 6.045e-12 | 3.406e-05 | 1.135e-03 |
| 1e-01 | 8.374e-12 | 3.822e-04 | 3.822e-03 |

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
| TokamakSmallCylindrical | 1.11e-15 | 1.61e-14 |
| TokamakSmallToroidal    | 4.77e-16 | 8.81e-16 |
| TokamakSmallCartesian   | 1.39e-13 | 1.76e-11 |
| SolovevIterXpoint       | 2.09e-15 | 5.31e-14 |

3D guiding centre, `HODEProblem` with PartitionedGauss(2). The informative column is not the energy
but the constraints: they vanish along the exact flow, so their magnitude is the drift off the
manifold on which `p = ϑ` — the quantity the approximately symplectic methods are designed to keep
bounded, and the reason `compute_constraints` was added to the diagnostics. The constraint the pair
omits is listed separately from the two it retains, because it is the retained pair's drift amplified
by `1/bₘ`.

| equilibrium | pair | \|ΔH/H\| | \|gᵏ\| retained | \|gᵏ\| omitted | \|Δp_φ/p_φ\| |
|---|---|---|---|---|---|
| SolovevIterXpoint        | `:g31` | 2.26e-15 | 6.54e-14 | 4.69e-11 | 0.00e+00 |
| TokamakMediumCartesian   | `:g23` | 3.21e-10 | 2.29e-09 | 9.94e-07 | — |
| TokamakSmallCylindrical  | `:g12` | 1.13e-14 | 1.14e-15 | 8.64e-19 | 0.00e+00 |
| TokamakMediumCylindrical | `:g12` | 3.78e-12 | 4.26e-12 | 8.63e-13 | 0.00e+00 |
| TokamakSmallToroidal     | `:g12` | 7.96e-16 | 6.37e-18 | 3.29e-21 | 0.00e+00 |

The three equilibria at the bottom are ones that could not be started at all before the constraint
pair became selectable. `SolovevSymmetricField`, which also could not, still cannot be run over a
thousand steps: it is the stiffest equilibrium in the package at `‖ϑ‖ ≈ 256` and the drift takes
Newton into a NaN a few hundred steps in. The `:g12` rows show the amplification working the *other*
way — `b₃` is the large component there, so the omitted constraint comes out better conserved than
the retained pair rather than worse. `TokamakMediumCartesian` is the case where it hurts, by a factor
of 430: `(g², g³)` divides by `b₂ = 0.9923` and so holds its retained pair to 2.3e-9, but the omitted
`g¹` carries `b₁`, which is zero at the initial condition and small along the whole orbit, and reaches
9.9e-7. Which way the amplification runs is decided by the *ratio* of the largest to the smallest
component of `b`, not by which pair was chosen.

`p_φ` is conserved *exactly* for the cylindrical and toroidal cases because φ is cyclic and the
partitioned symplectic method conserves the conjugate momentum of a cyclic coordinate to round-off.

Gyrokinetic guiding centre, over the module's default span (~200 units of physical time):

| method | \|ΔH/H\| |
|---|---|
| Gauss(2) on the full field | 2.67e-10 |
| Strang composition of the splitting | 3.72e-05 |

**The splitting is markedly worse on energy than a symplectic method of comparable cost.** That is
not a bug — it is volume preserving, not symplectic, so nothing constrains its energy error — but
it is five orders of magnitude, and it means the trade the scheme makes is a real one rather than
free. Over the very short orbits used in the tests both sit near 1e-14, so this only shows up at
realistic integration lengths.

Whether the energy error is *bounded* or *drifts* over long runs is the question this does not
answer, and is the study most worth doing next. See `TODO.md` in the repository root.


## The three families on one initial condition

`scripts/study_model_agreement.jl`

Three model families describe the guiding centre, and until the initial conditions were aligned nothing
compared them. `PauliParticle3d` declared `initial_conditions_*` for three of its eight equilibria, was
missing `barely_passing` on one of those three, and `GuidingCenter3d` carried a cartesian `y = 0.1`
offset on three modules that its 4D counterpart did not. All three families now declare the same
`(x, u, μ)` for every equilibrium they share, and `test/model_agreement_tests.jl` asserts it.

The comparison runs with `periodic = false`. The guiding centre problems wrap their angular coordinate
into the chart's range and the Pauli problems do not, so with wrapping left on a passing orbit's `ϕ` is
compared against the same `ϕ` plus a multiple of 2π — which reads as complete disagreement, and is how
two of these rows first looked.

### `GuidingCenter3d` and `GuidingCenter4d` are the same model

The 4D form carries the state `(x, u)`; the 3D form carries the position and its conjugate momentum
`p = ϑ(x, u)`, with the parallel velocity eliminated by the constraints. A disagreement would be a bug,
so this is an equivalence check rather than a comparison. Relative difference of the final `(x, u)` over
a thousand steps, `PartitionedGauss(2)` against `Gauss(2)`:

| equilibrium | chart | Δx | Δu |
|---|---|---|---|
| TokamakSmallCartesian   | cartesian   | 1.3e-10 – 3.7e-08 | 3.5e-10 – 3.3e-09 |
| TokamakSmallCylindrical | cylindrical | 9.3e-14 – 1.3e-11 | 7.6e-13 – 5.4e-12 |
| TokamakSmallToroidal    | toroidal    | 7.1e-15 – 1.8e-13 | 2.9e-14 – 1.0e-13 |
| TokamakIterCylindrical  | cylindrical | 6.5e-15 – 2.6e-14 | 1.0e-14 – 8.0e-14 |
| SolovevIter             | cylindrical | 2.9e-16 – 4.6e-15 | 0.0e+00 – 4.2e-15 |
| SolovevIterXpoint       | cylindrical | 4.1e-17 – 2.2e-16 | 8.3e-16 – 2.5e-15 |

Ranges over the four shipped conditions. The curvilinear charts agree at the accuracy of the nonlinear
solves. The cartesian chart is three to five orders looser, and that is a property of the chart rather
than of the equivalence — there the toroidal transit appears as a rotation of `(x, y)` that has to be
resolved, where in the other charts it is a coordinate that simply winds. The same asymmetry is why the
cartesian chart of this tokamak carries `Δt = 10` where its siblings carry 500.

### The Pauli particle differs by its initial condition, not its dynamics

`PauliParticle3d` is a *different* model — a regular Lagrangian on the full six-dimensional phasespace,
differing from the charged particle only by `μ|B|` — whose slow manifold is the guiding centre. So it
agrees only to the order at which the initial condition is placed on that manifold, and the interesting
question is which of the two that residual measures.

`v = u b⃗` is the manifold to lowest order. Against the 4D guiding centre, with `ρ/L` the gyroradius
over the field's gradient scale length:

| equilibrium | ρ/L | Δx at `v = u b⃗` | Δx at `v = v_gc` | gain |
|---|---|---|---|---|
| TokamakSmallCylindrical | 2.1e-03 | 1.7e-06 – 1.0e-05 | 1.2e-06 – 2.8e-05 | 0.1 – 8.7 |
| TokamakSmallToroidal    | 2.1e-03 | 5.6e-06 – 1.5e-05 | 4.2e-06 – 2.7e-05 | 0.2 – 3.5 |
| TokamakSmallCartesian   | 2.1e-03 | 4.5e-06 – 1.6e-04 | 1.8e-05 – 1.6e-04 | 0.1 – 8.7 |
| TokamakIterCylindrical  | 1.7e-02 | 2.1e-04 – 2.6e-04 | 6.8e-06 – 1.7e-05 | 14 – 37 |
| SolovevIter             | 2.4e-01 | 2.1e-03 – 1.8e-02 | 2.1e-06 – 2.8e-05 | 650 – 997 |
| SolovevIterXpoint       | 2.4e-01 | 2.1e-03 – 1.8e-02 | 1.9e-06 – 2.2e-05 | 811 – 1080 |

In the third column the agreement is ordered by `ρ/L`, as the expansion says it should be. But **within
one equilibrium it is not a power of either small parameter**: holding position and span fixed on
`SolovevIterXpoint` and taking `μ` down by four decades leaves it at 1.7e-2, and taking `u` down by
four decades leaves it at 1.5e-3. Both saturate rather than converge, so something that is not a
gyroradius effect sets the floor.

**It is the initial condition.** The guiding centre does not move along `b`; it moves along `b` plus the
∇B and curvature drifts. So `u b⃗` and the guiding centre velocity differ by the drift velocity, the
Pauli particle starts that far off the manifold, and it carries the difference as a gyration for the
whole orbit — an offset that does not shrink when `μ` or `u` does, because the drift shrinks with them.
Supplying the guiding centre velocity itself as the Pauli initial velocity is the fourth column, and
**the spread across equilibria collapses from four decades to two**: everything lands between 1.9e-6 and
1.8e-4, and the ITER Solov'ev, three orders worse than the small tokamak before, is no longer
distinguishable from it. That is the evidence that the two models really do share a slow manifold.

The correction is only a partial one, and the small tokamak rows say so: there it is 2–10× *worse*, and
the sign is consistent across charts and conditions rather than noise. Supplying `v_gc` fixes the
`O(v_drift)` mismatch in the initial *velocity*, but the Pauli state is the *particle* position, which
differs from the guiding centre position by the gyroradius vector, and that `O(ρ)` mismatch is left in
place. Where the velocity term dominates the correction wins by three orders; where `ρ/L` is small
enough that the residual is already down at the position term, correcting one without the other moves
the total the wrong way. Both together would be the first-order slow manifold proper, and
`PauliParticle3d.initial_conditions(x, u, μ)` does neither — it would have to depend on
`GuidingCenter4d` for a drift formula, and the lowest-order condition is the self-contained and
conventional one. See `TODO.md`.


## The residual scale of the ITER-size equilibria

`scripts/study_solver_tolerances.jl`

An *absolute* solver tolerance has to sit above the round-off floor of the residual it is applied
to, and that floor is a property of the equilibrium. The variational residuals carry the momentum
equation `p = ϑ(q)`, so the floor is `‖ϑ‖ eps`, evaluated here at each equilibrium's shipped
`initial_conditions_barely_passing`:

| family | equilibrium | ‖ϑ‖ or ‖p‖ | `‖ϑ‖ eps` | `f_abstol = 1E-15` |
|---|---|---|---|---|
| GC4d    | SolovevIterXpoint      | 14.94   | 3.3e-15 | unreachable |
| GC4d    | SolovevIter            | 14.94   | 3.3e-15 | unreachable |
| GC4d    | SolovevSymmetricField  | 256.3   | 5.7e-14 | unreachable |
| GC4d    | TokamakIterCylindrical | 30.3    | 6.7e-15 | unreachable |
| GC4d    | TokamakMediumCartesian | 1.17    | 2.6e-16 | reachable |
| GC4d    | TokamakSmallCartesian  | 0.0244  | 5.4e-18 | reachable |
| Pauli3d | SolovevIter            | 0.137   | 3.0e-17 | reachable |
| Pauli3d | SolovevIterXpoint      | 0.137   | 3.0e-17 | reachable |
| Pauli3d | TokamakIterCylindrical | 0.1835  | 4.1e-17 | reachable |

The last two Pauli rows are new: those modules had no `initial_conditions_barely_passing` to evaluate
at until the three families' initial conditions were aligned. `‖ϑ‖` moved only in the third decimal
for the affected equilibria, `ϑ = A + u b` being dominated by `‖A‖`, so the floors and the conclusion
below are the same as before the field reversal.

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
| `f_abstol=1E-15, f_reltol=1E-15` | 0.128 | 3.7 | 0.101 | 3.3 |
| `f_abstol=1E-15`                 | 0.117 | 3.4 | 0.040 | 3.0 |
| `f_abstol=1E-14`                 | 0.046 | 1.0 | 0.031 | 2.0 |
| `f_abstol=1E-12`                 | 0.046 | 1.0 | 0.031 | 2.0 |

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
   and 0.6.3, and so is every trajectory integrated from them.

`Project.toml` now requires 0.7.0, and the A/B/C/D table above is the one set of measurements on this
page that cannot be re-taken against it: A, B and D are states of two source trees that no longer
exist together. It is left as measured. What can be said is that 0.7.0 changes the *sign* of terms in
the generated code and not the amount of arithmetic, and that re-running the C column against it
reproduces every row to within the run-to-run spread of the harness — 0.086 against 0.082 for the
GC3d ITER Solov'ev `hodeproblem`, 0.106 against 0.104 for the medium cartesian one, 0.047 against
0.042 for the GC4d `iode`. Unlike 0.6.3, though, 0.7.0 is **not** value-preserving: it reverses `b` on
six of the eleven equilibria, so trajectories do change there. See the first section of this page.
