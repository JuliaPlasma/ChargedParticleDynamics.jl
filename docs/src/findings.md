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


## Conditioning of the 3D guiding centre constraint pair

`scripts/study_guiding_center_3d_conditioning.jl`

The constrained canonical formulation of Li, Zhang & Liu replaces the parallel-velocity constraint
by two of the three independent components of `v × b = 0`,

```
g¹ = b₁ v₃ - b₃ v₁ ,    g² = b₂ v₃ - b₃ v₂ ,    g³ = b₁ v₂ - b₂ v₁ ,
```

and the Poisson bracket in the denominator of both Lagrange multipliers is
`b₃ [B + (p-A)·(∇×b)]` for the pair `(g¹, g²)` and `b₁ [ … ]` for `(g¹, g³)`. This package uses
`(g³, g¹)`, so it is singular wherever `b₁ = 0`.

At the shipped initial condition of each equilibrium:

| equilibrium | coords | b₁ | λₒ | λ₁ | λ₂ |
|---|---|---|---|---|---|
| Dipole3d                 | cartesian   | -4.082e-01 | -3.402e+01 | 7.779e-03 | 5.874e-15 |
| QuadraticPotentials3d    | cartesian   | 2.000e-03  | 2.000e-01  | 1.514e+00 | 6.500e-03 |
| TokamakSmallCartesian    | cartesian   | -0.000e+00 | 0.000e+00  | **-Inf**  | **Inf**   |
| TokamakMediumCartesian   | cartesian   | -3.966e-02 | -1.537e-01 | 5.232e-02 | -4.172e-01 |
| TokamakSmallCylindrical  | cylindrical | 0.000e+00  | 0.000e+00  | **Inf**   | **-Inf**  |
| TokamakMediumCylindrical | cylindrical | 0.000e+00  | 0.000e+00  | **Inf**   | **-Inf**  |
| TokamakIterCylindrical   | cylindrical | 0.000e+00  | 0.000e+00  | **Inf**   | **Inf**   |
| SolovevIter              | cylindrical | 0.000e+00  | 0.000e+00  | **Inf**   | **Inf**   |
| SolovevIterXpoint        | cylindrical | 5.926e-03  | 1.207e+00  | 7.113e-01 | 7.470e-03 |
| SolovevSymmetricField    | cylindrical | -0.000e+00 | 0.000e+00  | **NaN**   | **-Inf**  |
| TokamakSmallToroidal     | toroidal    | 0.000e+00  | 0.000e+00  | **Inf**   | **-Inf**  |

**Seven of eleven equilibria have `b₁ = 0` exactly at the initial condition the package ships**, so
`λₒ` vanishes and the multipliers are infinite. `hode` cannot be started there at all — the Newton
solver hits a NaN in the direction vector on the first step.

Two things about this are worth recording:

* **It is not a curvilinear-coordinates problem.** `TokamakSmallCartesian` fails too: its initial
  condition sits at `y = z = 0`, where `B_x = 0`. What decides the outcome is where the initial
  condition lies relative to the field, not the coordinate system. The four survivors are simply
  the ones whose initial condition is off the symmetry plane.
* **The survivors are exactly the equilibria the test suite exercises**, plus the two toy problems
  that have their own scripts. That is not a coincidence — the other 3D test blocks were commented
  out, and `scripts/guiding_center_3d.jl` displaces its initial condition off the midplane by
  `sqrt(eps())`, both of which are workarounds for this.

Switching to `(g¹, g²)` puts `b₃ = b_φ` in the denominator — the largest component of a tokamak
field, and non-zero on the midplane — which would make most of these equilibria usable. The choice
is hard-coded in `guiding_center_3d_equations.jl`, with the alternative commented out beside it.
See `TODO.md` in the repository root.


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
same generalised vector potential — but the gyrokinetic module clears the phase-space Jacobian
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
symplectic in the remaining two, so a symplectic method applied to each preserves phase-space
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
but the two constraints: they vanish along the exact flow, so their magnitude is the drift off the
manifold on which `p = ϑ` — the quantity the approximately symplectic methods are designed to keep
bounded, and the reason `compute_constraints` was added to the diagnostics.

| equilibrium | \|ΔH/H\| | \|g₁\| | \|g₂\| | \|Δp_φ/p_φ\| |
|---|---|---|---|---|
| SolovevIterXpoint      | 5.74e-15 | 1.13e-15 | 1.47e-13 | 0.00e+00 |
| TokamakMediumCartesian | 3.23e-10 | 7.07e-10 | 1.78e-09 | — |

Only the two equilibria that can be started at all appear; see the conditioning section above.
`p_φ` is conserved *exactly* for the cylindrical case because φ is cyclic and the partitioned
symplectic method conserves the conjugate momentum of a cyclic coordinate to round-off.

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
