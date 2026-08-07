
# Gyrokinetic Guiding Centre Dynamics in 4D

`GyroKinetics4d` integrates the characteristics of the gyrokinetic Vlasov equation with a
*volume-preserving* splitting method, rather than with the symplectic or variational methods used
by the other models in this package.

It follows the notes *A volume preserving particle pusher for gyrokinetic PIC codes*.


## The characteristics

With gyrocentre phasespace variables ``Z = (R, V_{\parallel}, \mu, \theta)``, the single-particle
Lagrangian of a species of mass ``m`` and charge ``q`` is
```math
L (Z, \dot{Z}) = q \, A^{\star} \cdot \dot{R} + \frac{m}{q} \, \mu \, \dot{\theta} - H ,
```
with the generalised vector potential and its curl
```math
A^{\star} (R, V_{\parallel}) = A_{0} + \frac{m V_{\parallel}}{q} \, b_{0} ,
\qquad
B^{\star} = \nabla \times A^{\star} ,
\qquad
B^{\star}_{\parallel} = b_{0} \cdot B^{\star} .
```
``\theta`` is ignorable and ``\mu`` invariant, so the dynamics reduces to the four variables
``(R, V_{\parallel})``. The Euler-Lagrange equations give
```math
\begin{aligned}
\dfrac{dR}{dt} &= \dfrac{1}{B^{\star}_{\parallel}} \bigg[ \dfrac{1}{m} \dfrac{\partial H}{\partial V_{\parallel}} B^{\star} - \dfrac{1}{q} \nabla H \times b_{0} \bigg] , \\
\dfrac{dV_{\parallel}}{dt} &= - \dfrac{1}{B^{\star}_{\parallel}} \, \dfrac{1}{m} \, \nabla H \cdot B^{\star} .
\end{aligned}
```
This is the same system the [4D guiding centre model](guiding_center_4d.md) integrates.


## Rescaled time

The phasespace volume element is ``dZ = B^{\star}_{\parallel} \, dR \, dV_{\parallel} \, d\mu \, d\theta``,
so it is ``B^{\star}_{\parallel} f`` rather than ``f`` that is transported in the flat measure.
Introducing ``\tilde{f} = B^{\star}_{\parallel} f`` clears the denominator above and leaves
```math
\begin{aligned}
\dfrac{dR}{ds} &= \dfrac{1}{m} \dfrac{\partial H}{\partial V_{\parallel}} B^{\star} - \dfrac{1}{q} \nabla H \times b_{0} , &
\dfrac{dV_{\parallel}}{ds} &= - \dfrac{1}{m} \nabla H \cdot B^{\star} ,
\end{aligned}
```
which is what this module integrates.

!!! warning "The independent variable is not the physical time"
    The right-hand side is the guiding centre vector field multiplied by
    ``B^{\star}_{\parallel}``, so the independent variable ``s`` is related to the physical time by
    ```math
    dt = B^{\star}_{\parallel} \, ds .
    ```
    For the supplied ITER-like equilibrium ``B^{\star}_{\parallel} \approx 2 \times 10^{2}``, so one
    unit of ``s`` is some two hundred units of physical time and the default time step is scaled
    to match. A step chosen for the 4D guiding centre model would be a factor of
    ``B^{\star}_{\parallel}`` too large here.

    Because ``B^{\star}_{\parallel}`` depends on ``(R, V_{\parallel})``, the relation between
    ``s`` and ``t`` is not linear along an orbit; recovering the physical time requires
    integrating ``dt/ds = B^{\star}_{\parallel}`` alongside. The *orbits* are unaffected, which is
    all a PIC pusher needs.

!!! note "In a curvilinear chart the factor is `ωabs = J B*∥`"
    The expressions above are written in the flat measure of a cartesian chart, where the phasespace
    Jacobian is ``B^{\star}_{\parallel}``. In a general chart it is the Liouville density
    ``\sqrt{\det \Omega} = J \, B^{\star}_{\parallel}``, which is what
    [`ωabs`](@ref ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dTokamakSmallCylindrical.ωabs)
    computes and what the implemented ``dt = \omega_{abs} \, ds`` uses. Its magnitude is therefore
    chart-dependent — the three charts of the small tokamak give `0.9511`, `0.9986` and `0.0499` for
    a common `B*∥ = 0.9511` — while its *sign* is not: it is a measure density and positive in every
    chart. Getting that sign wrong reverses the orbit, which is what happened before each module
    declared its `ORIENTATION`; see [Findings](@ref) and `ωabs`.

The point of the rescaling is that the new right-hand side is divergence-free,
```math
\nabla \cdot \bigg( B^{\star}_{\parallel} \dfrac{dR}{dt} \bigg)
+ \dfrac{\partial}{\partial V_{\parallel}} \bigg( B^{\star}_{\parallel} \dfrac{dV_{\parallel}}{dt} \bigg) = 0 ,
```
so its flow preserves phasespace volume, and a suitable discretisation can preserve it exactly.


## Potentials and splitting

Defining the two vector potentials
```math
\beta = \frac{1}{m} \frac{\partial H}{\partial V_{\parallel}} A^{\star} ,
\qquad
\gamma = \frac{1}{m} A^{\star} \times \nabla H ,
```
the equations of motion take the manifestly divergence-free form
```math
\begin{aligned}
\dfrac{dR}{ds} &= \dfrac{\partial \gamma}{\partial V_{\parallel}} + \nabla \times \beta , &
\dfrac{dV_{\parallel}}{ds} &= - \nabla \cdot \gamma .
\end{aligned}
```
For the equilibrium Hamiltonian ``H_{0} = \tfrac{1}{2} m V_{\parallel}^{2} + \mu B_{0}``, which is
what is implemented here, these reduce to
```math
\beta_{0} = V_{\parallel} A^{\star} ,
\qquad
\gamma_{0} = \frac{\mu}{m} A^{\star} \times \nabla B_{0} .
```

This form splits into six subsystems,
```math
\begin{aligned}
\dot{R}_{1} &= + \partial_{2} \beta_{3} , & \dot{R}_{1} &= - \partial_{3} \beta_{2} , & \dot{R}_{2} &= + \partial_{3} \beta_{1} , \\
\dot{R}_{2} &= - \partial_{1} \beta_{3} , & \dot{R}_{3} &= + \partial_{1} \beta_{2} , & \dot{R}_{3} &= - \partial_{2} \beta_{1} ,
\end{aligned}
```
```math
\begin{aligned}
\dot{R}_{i} &= + \partial_{V_{\parallel}} \gamma_{i} , &
\dot{V}_{\parallel} &= - \partial_{i} \gamma_{i} , &
i &= 1, 2, 3 ,
\end{aligned}
```
each freezing two of the four variables and symplectic in the remaining two. Applying a symplectic
integrator to each subsystem preserves the volume of that subsystem exactly, so any composition of
them is volume preserving; a symmetric composition of second-order methods gives a second-order
volume-preserving scheme.

The full field is `odeproblem`, the splitting `sodeproblem`.


## Noncanonical Hamiltonian structure

The same equations can be written as ``\dot{z} = \mathcal{J}(z) \nabla H`` with the Poisson matrix
```math
\mathcal{J} = \begin{pmatrix}
0 & - b_{3}/q & + b_{2}/q & + B^{\star}_{1}/m \\
+ b_{3}/q & 0 & - b_{1}/q & + B^{\star}_{2}/m \\
- b_{2}/q & + b_{1}/q & 0 & + B^{\star}_{3}/m \\
- B^{\star}_{1}/m & - B^{\star}_{2}/m & - B^{\star}_{3}/m & 0
\end{pmatrix} ,
```
where ``\mathcal{J} / B^{\star}_{\parallel}`` is the inverse of the symplectic matrix of the
guiding centre model. This structure is not used by the integrators here — it is the
volume-preserving splitting above that the module implements — but `ω` provides the corresponding
two-form.


## Usage

```julia
using GeometricIntegrators
using ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dSolovevIterXpoint

# volume-preserving Strang composition of the six subsystems
problem  = sodeproblem(initial_conditions_deeply_passing())
solution = integrate(problem, Composition(Tuple(Gauss(1) for _ in 1:6), Strang()))

# or the full vector field
problem  = odeproblem(initial_conditions_deeply_passing())
solution = integrate(problem, Gauss(2))
```


## Status

The vector field is verified against the [4D guiding centre model](guiding_center_4d.md): the two
agree to machine precision after multiplication by ``B^{\star}_{\parallel}``, the six subsystems
sum exactly to the full field, and the Strang composition holds the determinant of the one-step
map at ``1`` for every step size tested, while a first-order method's volume error grows linearly
with the step.

Two limitations are worth recording:

* Only the equilibrium Hamiltonian ``H_{0}`` is implemented. The notes also give the zero
  Larmor radius Hamiltonian, with the ``\phi``, ``A_{\parallel}`` and
  ``\langle \tilde{A}_{\parallel} \rangle^{2}`` terms; those are absent here, as they are from the
  4D guiding centre model (see the [Model Audit](@ref)).
* Only the ITER-like Solov'ev equilibrium with X-point is provided, where the other model families
  cover a dozen equilibria each.

The coordinate transformation of `coordinate_transformations.jl` is not applied by the problem
constructors, and is retained only as a utility: it is a preconditioner for the nonlinear solver
rather than a change of variables in the model, so it belongs inside an integrator. No such
integrator exists yet; see `TODO.md` in the repository root for the design.


## Equilibria

```@docs
ChargedParticleDynamics.GyroKinetics4d
```

Each equilibrium is its own module, and each one `include`s the same `gc_common.jl`,
`gc_equations.jl` and `coordinate_transformations.jl`, so the model's functions below are documented
once, under `GuidingCenter4dTokamakSmallCylindrical`, and hold verbatim for all eight. What differs
between them is the chart, the equilibrium parameters and the rescaled default time step, which is
what these docstrings record:

```@docs
ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dTokamakSmallCartesian
ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dTokamakSmallToroidal
ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dTokamakMediumCartesian
ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dTokamakMediumCylindrical
ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dTokamakIterCylindrical
ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dSolovevIter
ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dSolovevIterXpoint
```


## Module

```@autodocs
Modules = [ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dTokamakSmallCylindrical]
```
