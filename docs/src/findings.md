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


## Conditioning of the 3D guiding centre constraint formulations

`scripts/study_guiding_center_3d_conditioning.jl`

The constrained canonical formulation of Li, Zhang & Liu replaces the parallel-velocity constraint
by two of the three independent components of `v × b = 0`,

```
g¹ = b₁ v₃ - b₃ v₁ ,    g² = b₂ v₃ - b₃ v₂ ,    g³ = b₁ v₂ - b₂ v₁ ,
```

and the Poisson bracket in the denominator of both Lagrange multipliers carries a factor of the
`b`-component belonging to the constraint the pair leaves out: `b₁` for `(g³, g¹)` (`constraints =
:g31`), `b₃` for `(g¹, g²)` (`:g12`), `b₂` for `(g², g³)` (`:g23`). The package once implemented
`(g³, g¹)` alone, so it was singular wherever `b₁ = 0`.


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
|                        | canonical `:g31`    | 7.480e-09 | 7.658e-09 | 1.329e-08 |
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
| canonical `:g31` | 3.27e-05 | 8.11e-07 | NaN in the Newton direction |
| canonical `:g12` | 5.62e-06 | 1.57e-08 | 1.57e-08 |

Four orders of magnitude at `t = 3`, and past `t = 30` the badly conditioned pair loses the orbit
altogether — `|ΔH/H| = 2e+03` is not an error estimate, it is a destroyed trajectory. With `:g12`
both the Hamilton-Dirac and the canonicalised form hold their levels out to `t = 300`.

`Dipole3d` is therefore the one equilibrium whose `default_constraints` is chosen on the conditioning
of the pair *along the orbit* rather than at the initial condition. The others were left where
regularity put them; see `TODO.md`.

The related effect is that the constraint the pair omits is conserved `1/bₘ` times worse than the two
it retains, since `b₁g² - b₂g¹ + b₃g³ = 0` determines it from them pointwise. Over a thousand steps
of `TokamakMediumCartesian` with `:g31` the retained pair holds to 2e-9 while the omitted `g²`
reaches 1.4e-6, the orbit having passed through `b₁ = 3.5e-4`. That is why `compute_constraints`
reports all three.


### What the formulations cost

Per step, taking the Hamilton-Dirac form as unity. Measured across all eleven equilibria; the spread
between them is a few percent.

| formulation | relative cost | why |
|---|---|---|
| compact, literal pair | 0.86 – 0.95 | replaces both multipliers by a single bracket ratio; no second derivatives |
| Hamilton-Dirac | 1 | two multipliers over six bracket terms |
| compact, `:parallel` | 1.21 – 1.44 | the regularised `D` averages all three brackets, so nine terms rather than three |
| canonicalised | 8.9 – 29.6 | `∂λ/∂q` and `∂λ/∂p` need every second derivative of the field |

The three pairs cost the same as each other to within noise — they are the same expressions with
permuted indices — so the pair should be chosen on conditioning alone. The order of magnitude on the
canonicalised form is the number worth remembering: it buys exact symplecticity rather than
approximate, and gives up two to three orders of magnitude of energy and constraint conservation for
it at these step sizes. Its spread is wider than the others because part of the cost is the harder
nonlinear solve, not only the right-hand side: `Dipole3d`, the worst case at 29.6, is also the only
one where Newton needs more than one iteration per stage on average.

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
`‖ϑ‖ ≈ 256`; it never showed up in CI only because its block runs `ode` rather than `iode`.

Two knobs are involved, and they fail for different models. Cost per step at the ITER Solov'ev
X-point, with the Newton iteration count behind it:

| setting | GC4d `iode` ms/step | mean iters | Pauli3d `iode` ms/step | mean iters |
|---|---|---|---|---|
| `f_abstol=1E-15, f_reltol=1E-15` | 11.89 | 4.0 | 0.74 | 3.9 |
| `f_abstol=1E-15`                 | 11.91 | 4.0 | 0.74 | 3.9 |
| `f_abstol=1E-14`                 |  0.33 | 2.0 | 0.12 | 3.0 |
| `f_abstol=1E-12`                 |  0.33 | 2.0 | 0.11 | 2.9 |

The 4D guiding centre degrades from `f_abstol = 1E-15` alone and recovers between 1E-15 and 1E-14,
exactly where its floor lies. The Pauli model has a small momentum and is untouched by `f_abstol`;
it degrades only once `f_reltol` is pinned as well, i.e. what hurts it is the loss of the relative
term `f_reltol ‖F(x₀)‖` that lets a large-magnitude solve converge at all. `f_abstol = 1E-12` with
`f_reltol` left at its default is the only setting comfortable for both, and is what the test suite
uses. **The final state is unchanged to ~1e-14 across every row** — the extra iterations buy
nothing.

Note this also rules out simply dropping the options: the library default is `f_abstol = 8 eps() =
1.8e-15`, itself below the ITER floor.

Once the tolerance is sound the remaining cost has nothing to do with ITER:

| workload | ms/step | mean iters |
|---|---|---|
| GC3d SolovevIterXpoint `hodeproblem` + extrapolation      | 17.9 | 1.0 |
| GC3d TokamakMediumCartesian `hodeproblem` + extrapolation | 23.8 | 1.0 |
| GC4d SolovevIterXpoint `iode`                      |  0.3 | 2.0 |
| Pauli3d SolovevIter `iode`                         |  0.1 | 2.9 |

`TokamakMediumCartesian` is *slower* than the ITER Solov'ev, and the 3D `hodeproblem` path costs two
orders of magnitude more per step than the Pauli one even though Newton converges in a single
iteration. The cost is nested `ForwardDiff`: the second derivatives of the 3D Hamiltonian
differentiate field functions that are themselves AD-generated, and `MidpointExtrapolation(5)`
triples it. Making those derivatives cheaper is the real fix; cutting the two 3D blocks from 1000
steps to 100 was the pragmatic one.
