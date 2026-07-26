
# Guiding Center Dynamics in 3D

The [four-dimensional guiding centre model](guiding_center_4d.md) evolves the guiding centre position ``r = (x,y,z)`` together with the parallel velocity ``u``.
The three-dimensional model instead treats the parallel velocity as a dependent variable and evolves only the position together with its canonically conjugate momentum, so that the guiding centre dynamics takes the form of a *constrained canonical* Hamiltonian system.
The implementation follows Xinjie Li, Ruili Zhang and Jian Liu, *Approximately symplectic
Runge-Kutta methods for the guiding center dynamics*, Journal of Computational Physics **563**,
115069 (2026), [doi:10.1016/j.jcp.2026.115069](https://doi.org/10.1016/j.jcp.2026.115069), and its
companion Ruili Zhang, Jian Liu, Tong Liu, Wenxiang Li and Xiaogang Wang, *Canonical Hamiltonian
guiding center theory and classical intrinsic magnetic moment*, Frontiers of Physics **21**(2),
026200 (2026).

Equation numbers below refer to the first of those two papers.


## Where the constraints come from

The starting point is the guiding centre one-form ``\vartheta(r,u) = A (r) + u \, b (r)``, where ``A`` is the magnetic vector potential, ``B = \nabla \times A`` and ``b = B / \vert B \vert``.
Identifying the momentum with the one-form,
```math
p = A (r) + u \, b (r) ,
```
the parallel velocity is recovered from ``p`` and ``r`` as
```math
v (r,p) = p - A (r) ,
\qquad
u (r,p) = b (r) \cdot v (r,p) ,
```
and the Hamiltonian becomes
```math
H (r,p) = \tfrac{1}{2} \, \vert v (r,p) \vert^{2} + \mu \, \vert B (r) \vert + \varphi (r) ,
```
with the magnetic moment ``\mu`` and the electrostatic potential ``\varphi``.

The Legendre transform is degenerate: ``\partial L / \partial \dot{r} = M(r) \dot{r} + A(r)`` with ``M = b \otimes b`` of rank one, so ``p`` has only one independent component along ``b`` and the three momenta are tied together by
```math
\frac{p_{1} - A_{1}}{b_{1}} = \frac{p_{2} - A_{2}}{b_{2}} = \frac{p_{3} - A_{3}}{b_{3}} ,
```
that is by ``v \times b = 0``. Two of the three components of that are independent, and Li, Zhang and Liu label the three as (Eqs. 4 and 5)
```math
\begin{aligned}
g^{1} (r,p) &= b_{1} \, v_{3} - b_{3} \, v_{1} , \\
g^{2} (r,p) &= b_{2} \, v_{3} - b_{3} \, v_{2} , \\
g^{3} (r,p) &= b_{1} \, v_{2} - b_{2} \, v_{1} .
\end{aligned}
```
They are the components of ``b \times v`` up to labelling and sign, ``c_{k} = \varepsilon_{kab} \, b_{a} v_{b}``, with
```math
c_{1} = g^{2} , \qquad c_{2} = - g^{1} , \qquad c_{3} = g^{3} ,
```
which is what makes them one indexed family rather than three separate expressions. In the source they are `gᵏ(Val(k), t, q, p)`, with `constraint_indices` supplying the pair of field indices each is built from, and every first and second derivative of them is written once against that index pair in `src/guiding_center_3d/guiding_center_3d_constraints.jl`.


## Choosing the constraint pair

Any two of the three ``g^{k}`` confine the solution to the same manifold, and the choice is not free in practice. Requiring ``\dot{g}_{1} = \dot{g}_{2} = 0`` determines the two Lagrange multipliers as (Eq. 10)
```math
\lambda_{1} = \frac{\{ g_{2} , H \}}{\{ g_{1} , g_{2} \}} ,
\qquad
\lambda_{2} = - \frac{\{ g_{1} , H \}}{\{ g_{1} , g_{2} \}} ,
```
and the Poisson bracket in the denominator evaluates to (Eq. 22)
```math
\{ g_{1} , g_{2} \} = b_{m} \left[ \vert B \vert + (p - A) \cdot (\nabla \times b) \right] ,
```
where ``m`` is the index of the constraint the pair leaves out. So each pair is singular on the surface where one particular component of ``b`` vanishes:

| `constraints` | pair | singular where |
|---|---|---|
| `:g31` | ``(g^{3}, g^{1})`` | ``b_{1} = 0`` |
| `:g12` | ``(g^{1}, g^{2})`` | ``b_{3} = 0`` |
| `:g23` | ``(g^{2}, g^{3})`` | ``b_{2} = 0`` |

This is not a fine point. In cylindrical coordinates ``b_{1} = b_{R}`` vanishes on the midplane of an axisymmetric equilibrium, which is exactly where most of the initial conditions in this package sit: seven of the eleven equilibria that ship an initial condition have ``b_{1} = 0`` there exactly, and `SolovevSymmetricField` has ``b_{1} = b_{3} = 0`` at once, leaving `:g23` as its only usable pair. The package once implemented ``(g^{3}, g^{1})`` alone and could not be started on any of them.

Every problem constructor therefore takes a `constraints` keyword, and every equilibrium declares its own `default_constraints` that is regular at its own initial conditions:

```julia
hodeproblem(initial_conditions_barely_passing())                      # the equilibrium's own pair
hodeproblem(initial_conditions_barely_passing(); constraints = :g12)  # ⟨g¹, g²⟩, singular at b₃ = 0
```

The symbol is resolved once, in the constructor, into the `Val`-based pair `constraint_pair` returns, so the right-hand side stays specialised on it. There is deliberately no run-time switching between pairs along an orbit. `scripts/study_guiding_center_3d_conditioning.jl` tabulates the conditioning of all three pairs for every equilibrium.

``b_{m} \ne 0`` at the initial condition is only the entry requirement. Along the orbit, how large ``b_{m}`` is decides how well the *omitted* constraint is conserved. The identity
```math
b_{1} g^{2} - b_{2} g^{1} + b_{3} g^{3} = b \cdot (b \times v) = 0
```
holds pointwise, so the omitted constraint is the retained pair divided by ``b_{m}``:
```math
\vert g_{\mathrm{omitted}} \vert \lesssim \frac{\max \vert g_{\mathrm{retained}} \vert}{\min \vert b_{m} \vert} ,
```
with the minimum over the orbit. On `TokamakMediumCartesian` with `:g31` the retained pair holds to ``2 \times 10^{-9}`` over a thousand steps while the omitted ``g^{2}`` reaches ``1.4 \times 10^{-6}``, because the orbit passes through ``b_{1} = 3.5 \times 10^{-4}``. Nothing is wrong there — the numerical flow enforces the pair it was given — but it is the reason `compute_constraints` reports all three components.

Where ``b_{m}`` gets small enough, the whole solution goes, not just the omitted constraint. `Dipole3d` is the case in this package: all three pairs are regular at its initial condition, but its orbit takes ``b_{1}`` and ``b_{2}`` through zero, and the two pairs that divide by them hold the constraints to ``10^{-3}`` rather than ``10^{-8}`` over a hundred steps and lose the orbit altogether by ``t = 30``. Its `default_constraints` is therefore chosen on the conditioning of the pair along the orbit; the other ten are chosen on regularity at the initial condition, which is a weaker criterion. `scripts/study_guiding_center_3d_conditioning.jl` §2 measures both.


## The three formulations

Dirac's theory gives three equivalent systems, and this package builds all three.

### Hamilton-Dirac — `hodeproblem`, Eq. (11)

With the multipliers substituted, the constraints are imposed only through the initial condition and the dynamics is an ordinary differential equation,
```math
\begin{aligned}
\dot{r} (t) &= + \dfrac{\partial H}{\partial p} + \lambda_{1} \dfrac{\partial g_{1}}{\partial p} + \lambda_{2} \dfrac{\partial g_{2}}{\partial p} , \\
\dot{p} (t) &= - \dfrac{\partial H}{\partial r} - \lambda_{1} \dfrac{\partial g_{1}}{\partial r} - \lambda_{2} \dfrac{\partial g_{2}}{\partial r} ,
\end{aligned}
```
which is a Hamiltonian system in the sense of a [`HODEProblem`](https://juliagni.github.io/GeometricEquations.jl/stable/). Its flow preserves ``H`` and all three ``g^{k}`` exactly, and the restriction to the constraint manifold is symplectic — which is why an ordinary symplectic Runge-Kutta method applied to it is *approximately* symplectic, with the defect bounded by the constraint drift.

### Canonicalised — `hodeproblem_canonical`, Eq. (13)

The multipliers can alternatively be absorbed into the Hamiltonian,
```math
\tilde{H} (r,p) = H (r,p) + \lambda_{1} (r,p) \, g_{1} (r,p) + \lambda_{2} (r,p) \, g_{2} (r,p) ,
```
which agrees with ``H`` on the constraint manifold. Differentiating ``\tilde{H}`` adds terms proportional to the constraints,
```math
\begin{aligned}
\dot{r} (t) &= \dfrac{\partial H}{\partial p} + \lambda_{1} \dfrac{\partial g_{1}}{\partial p} + \lambda_{2} \dfrac{\partial g_{2}}{\partial p} + \dfrac{\partial \lambda_{1}}{\partial p} g_{1} + \dfrac{\partial \lambda_{2}}{\partial p} g_{2} , \\
\dot{p} (t) &= - \dfrac{\partial H}{\partial r} - \lambda_{1} \dfrac{\partial g_{1}}{\partial r} - \lambda_{2} \dfrac{\partial g_{2}}{\partial r} - \dfrac{\partial \lambda_{1}}{\partial r} g_{1} - \dfrac{\partial \lambda_{2}}{\partial r} g_{2} ,
\end{aligned}
```
which vanish for exact initial data but not for a numerical solution that has drifted off the manifold. This system is genuinely canonical, so a symplectic method applied to it is symplectic exactly. The price is the second derivatives of the field, needed for ``\partial \lambda / \partial r`` and ``\partial \lambda / \partial p``, and a right-hand side that amplifies rather than ignores the constraint drift. It is by a wide margin the most expensive of the three — nine to thirty times the Hamilton-Dirac form per step, the wide spread coming from the harder nonlinear solve as much as from the right-hand side — and the least robust: at the step sizes the test suite uses it conserves energy and the constraints two to three orders of magnitude *worse* than the form it was derived from, and on `SolovevSymmetricField` it is the one formulation that cannot complete an orbit at all. Section 3 of `scripts/study_guiding_center_3d_conditioning.jl` has the numbers.

### Compact — `hodeproblem_compact`, Eq. (29)

Rearranged, the Hamilton-Dirac right-hand side splits into a part that is regular where ``b_{m} = 0`` and a part proportional to the constraints ``g^{j}``, which vanish along the solution. Dropping the latter removes most of the singularity, and is the second of the paper's two remedies for it. Not all of it, as it turns out: what Eq. (29) is left with still divides by ``b_{m}`` once, in the factor ``(p_{m} - A_{m}) / b_{m}``. Removing that too is what makes the implementation here depart from the paper; see below.

Eq. (29) as printed cannot be ported here. It is written with cartesian vector identities — ``\nabla \times b``, ``b \times \nabla B`` — and with ``H = \tfrac{1}{2} \sum_{i} (p_{i} - A_{i})^{2}``, whereas eight of the thirteen equilibria in this package are curvilinear and its Hamiltonian carries ``g^{11}, g^{22}, g^{33}``. The implementation therefore derives the same system from the model's own objects, coordinate-generally. The derivation is recorded in full at the top of `src/guiding_center_3d/guiding_center_3d_compact.jl`; in outline:

* Written in the antisymmetric labelling ``c = b \times v``, the momentum derivative of a constraint is
  ```math
  \frac{\partial c_{k}}{\partial p_{l}} = \varepsilon_{kal} \, b_{a} ,
  ```
  with no metric in it.
* For the pair retaining ``(c_{i}, c_{j})``, the cyclic complement of ``m``, the multiplier contribution to ``\dot{r}`` collapses to
  ```math
  C_{l} = \sum_{k} \lambda^{k} \frac{\partial c_{k}}{\partial p_{l}}
        = - \frac{b_{m} \, \{ c_{l} , H \}}{\{ c_{i} , c_{j} \}}
  ```
  uniformly in ``l``. For ``l \in \{i, j\}`` that is exact algebra; for ``l = m`` it needs ``\sum_{l} b_{l} c_{l} = b \cdot (b \times v) = 0``, which holds identically and hence gives ``\sum_{l} b_{l} \{ c_{l}, H \} = 0`` on the constraint manifold. That identity is the one thing Eq. (29) uses the constraints for.
* Since ``\{ c_{i}, c_{j} \} = b_{m} D`` with the same ``D`` for all three ``m``, the singular factor ``b_{m} / \{ c_{i}, c_{j} \}`` may be replaced by ``1/D`` with
  ```math
  D = \frac{\sum_{m} b_{m} \, \{ c_{i} , c_{j} \}(m)}{\sum_{m} b_{m} b_{m}} ,
  ```
  whose denominator cannot vanish because ``\vert b \vert = 1``. This is `compact_denominator`, and in a cartesian chart it reduces to ``\vert B \vert + (p-A) \cdot (\nabla \times b)`` as it should.
* The momentum equation follows the same way. Splitting ``\partial c_{k} / \partial q_{l}`` into its ``\partial b / \partial x`` and ``\partial A / \partial x`` parts, the second contracts straight back into ``C``, and the first does too once ``v_{a} = u \, b_{a}`` is used on the manifold. What is left is
  ```math
  \dot{r}_{l} = \frac{\partial H}{\partial p_{l}} + C_{l} ,
  \qquad
  \dot{p}_{l} = - \frac{\partial H}{\partial q_{l}} + \sum_{a} \left( \frac{\partial A_{a}}{\partial x_{l}} + u \, \frac{\partial b_{a}}{\partial x_{l}} \right) C_{a} .
  ```

Two things about this are worth stating plainly, because they are departures from the paper rather than restatements of it.

First, **the compact form is independent of the constraint pair.** Eq. (29) depends on the pair only through the factor ``(p_{m} - A_{m}) / b_{m}``, which is the parallel velocity in three different spellings, all equal on the constraint manifold. Using ``u = b \cdot v`` — the `u(t, q, p)` the model already has, regular everywhere — makes the dependence disappear. This is what `constraints = :parallel`, the default of `hodeproblem_compact`, selects.

Second, **the literal Eq. (29) is still singular.** ``(p_{m} - A_{m}) / b_{m}`` is ``0/0`` exactly where the pair was singular to begin with, so the ``\nu_{j}`` and ``\omega_{j}`` terms are not the only obstruction. Passing `:g31`, `:g12` or `:g23` to `hodeproblem_compact` reproduces the paper's expressions for that pair, singularity included; section 1 of `scripts/study_guiding_center_3d_conditioning.jl` shows both side by side.

The equivalence is pinned by a test: `test/structure_tests.jl` spells Eq. (29) out directly from the injected field functions for the three cartesian equilibria — including `QuadraticPotentials3d`, the one with a non-zero electrostatic potential — and checks it against `guiding_center_3d_compact_v`/`_f` to round-off.

On cost, the literal form is the cheapest of the three formulations, 0.86 to 0.95 of the Hamilton-Dirac form per step across the eleven equilibria: it replaces the two multipliers by a single bracket ratio and needs no second derivatives at all. The regularisation costs about half as much again, because ``D`` averages the bracket over all three ``m`` and so evaluates nine bracket terms where the literal form evaluates three. That leaves `:parallel` at 1.2 to 1.45 times the Hamilton-Dirac form — a cheap price for never dividing by a component of ``b``, and an order of magnitude below the canonicalised system.


## Which symbol is which

| source | equation |
|---|---|
| `gᵏ(Val(k), t, q, p)` | ``g^{k}``, Eqs. (4) and (5); `g₁`, `g₂`, `g₃` are the three fixed values |
| `cₗ(Val(l), t, q, p)` | ``c_{l} = (b \times v)_{l}``, the same three in antisymmetric labelling |
| `dgᵏdqₗ`, `dgᵏdpₗ`, `d²gᵏdqₗdqₘ`, `d²gᵏdqₗdpₘ` | the derivatives of ``g^{k}``; ``\partial^{2} g^{k} / \partial p \partial p`` vanishes identically |
| `bracket_gg`, `bracket_gH` | ``\{ g^{i}, g^{j} \}`` and ``\{ g^{k}, H \}``, Eqs. (22) and (23) |
| `λₒ`, `λ₁`, `λ₂` | ``\{ g_{1}, g_{2} \}`` and the two multipliers, Eq. (10) |
| `dλ₁dqⱼ`, `dλ₁dpⱼ`, `dλ₂dqⱼ`, `dλ₂dpⱼ` | the multiplier gradients of Eq. (12) |
| `bracket_cH`, `bracket_cc` | ``\{ c_{l}, H \}`` and ``\{ c_{i}, c_{j} \}``, the compact form's ingredients |
| `compact_denominator` | ``D``, Eq. (24) |
| `compact_multipliers` | ``C_{l}``, the multiplier contribution to ``\dot{r}``, all three at once |
| `constraint_pair`, `compact_index` | the `constraints` keyword resolved to indices |


## Usage

Each equilibrium is wrapped in its own module, which injects the magnetic field code and provides
`initial_conditions` for converting ``(r,u)`` into ``(r,p)`` as well as `hodeproblem`,
`hodeproblem_canonical` and `hodeproblem_compact`:

```julia
using GeometricIntegrators
using ChargedParticleDynamics.GuidingCenter3d.Dipole3d

problem = hodeproblem(initial_conditions_dipole())
solution = integrate(problem, PartitionedGauss(1))
```

`compute_constraints(solution)` returns all three ``g^{k}`` along the orbit, whichever pair the
problem retained, since all three vanish along the exact flow.


## Modules

```@docs
ChargedParticleDynamics.GuidingCenter3d
```

Each equilibrium is its own module, and every one of them `include`s the same
`guiding_center_3d_equations.jl`, `guiding_center_3d_canonical.jl` and
`guiding_center_3d_diagnostics.jl` — the first of which pulls in
`guiding_center_3d_constraints.jl` and the second `guiding_center_3d_compact.jl`. The model's
functions are therefore documented once, below, under `TokamakSmallCylindrical`, and hold verbatim
for all thirteen. What differs between the modules is the chart, the equilibrium parameters, the
initial conditions and the constraint pair they default to, which is what these docstrings record:

```@docs
ChargedParticleDynamics.GuidingCenter3d.Dipole3d
ChargedParticleDynamics.GuidingCenter3d.QuadraticPotentials3d
ChargedParticleDynamics.GuidingCenter3d.SymmetricField
ChargedParticleDynamics.GuidingCenter3d.ThetaPinchField
ChargedParticleDynamics.GuidingCenter3d.TokamakSmallCartesian
ChargedParticleDynamics.GuidingCenter3d.TokamakSmallToroidal
ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCartesian
ChargedParticleDynamics.GuidingCenter3d.TokamakMediumCylindrical
ChargedParticleDynamics.GuidingCenter3d.TokamakIterCylindrical
ChargedParticleDynamics.GuidingCenter3d.SolovevSymmetricField
ChargedParticleDynamics.GuidingCenter3d.SolovevIter
ChargedParticleDynamics.GuidingCenter3d.SolovevIterXpoint
```

```@autodocs
Modules = [ChargedParticleDynamics.GuidingCenter3d.TokamakSmallCylindrical]
```
