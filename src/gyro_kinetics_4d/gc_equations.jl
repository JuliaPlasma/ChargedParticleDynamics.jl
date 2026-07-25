
import ElectromagneticFields

using GeometricEquations: ODEProblem, SODEProblem

export guiding_center_4d_ode, guiding_center_4d_sode


# `GeometricEquations` recognises periodicity only as an `(xmin, xmax)` tuple of arrays; the range
# of the position coordinates comes from the coordinate system, via the `rangemin`/`rangemax`
# injected with the field code. The parallel velocity is not periodic.
function guiding_center_4d_periodicity(::Type{T}, periodic=true) where {T}
    xmin = -T(Inf) * ones(T, 4)
    xmax = +T(Inf) * ones(T, 4)

    if periodic
        xmin[1:3] .= rangemin(xmin[1:3])
        xmax[1:3] .= rangemax(xmax[1:3])
    end

    return (xmin, xmax)
end

guiding_center_4d_periodicity(::AbstractVector{<:AbstractArray{T}}, periodic=true) where {T<:Number} = guiding_center_4d_periodicity(T, periodic)
guiding_center_4d_periodicity(::AbstractArray{T}, periodic=true) where {T<:Number} = guiding_center_4d_periodicity(T, periodic)


@doc raw"""
    guiding_center_4d_ode(qᵢ, parameters; kwargs...)

The gyrokinetic guiding centre characteristics in the rescaled time of the volume-preserving
formulation,

```math
\begin{aligned}
\dfrac{d R}{d s} &= \dfrac{\partial \gamma}{\partial u} + \nabla \times \beta , &
\dfrac{d u}{d s} &= - \nabla \cdot \gamma ,
\end{aligned}
```

with the vector potentials ``\beta`` and ``\gamma`` of [`β`](@ref) and [`γ`](@ref).

The independent variable is **not** the physical time: the right-hand side is the guiding centre
vector field multiplied by ``B^{\star}_{\parallel}``, so ``dt = B^{\star}_{\parallel} \, ds``. One
unit of ``s`` therefore covers ``B^{\star}_{\parallel} \approx 200`` units of physical time for the
supplied ITER-like equilibrium, and the default time step is scaled accordingly.

That factor is the phase-space Jacobian, absorbed into the distribution function as
``\tilde{f} = B^{\star}_{\parallel} f``. What it buys is that the right-hand side is then
divergence-free, which is what makes the splitting of [`guiding_center_4d_sode`](@ref) volume
preserving.
"""
function guiding_center_4d_ode(qᵢ=qᵢ, parameters=parameters; tspan=tspan, tstep=Δt, periodic=true)
    ODEProblem(
        v,
        tspan, tstep, qᵢ;
        parameters=parameters,
        invariants=(h=hamiltonian,),
        periodicity=guiding_center_4d_periodicity(qᵢ, periodic)
    )
end


@doc raw"""
    guiding_center_4d_sode(qᵢ, parameters; kwargs...)

The same characteristics as [`guiding_center_4d_ode`](@ref), split into the six subsystems

```math
\begin{aligned}
\dot{R}_{1} &= + \partial_{2} \beta_{3} , & \dot{R}_{1} &= - \partial_{3} \beta_{2} , & \dot{R}_{2} &= + \partial_{3} \beta_{1} , \\
\dot{R}_{2} &= - \partial_{1} \beta_{3} , & \dot{R}_{3} &= + \partial_{1} \beta_{2} , & \dot{R}_{3} &= - \partial_{2} \beta_{1} ,
\end{aligned}
```

```math
\begin{aligned}
\dot{R}_{i} &= + \partial_{u} \gamma_{i} , &
\dot{u} &= - \partial_{i} \gamma_{i} , &
i &= 1,2,3 ,
\end{aligned}
```

each of which freezes two of the four variables and is symplectic in the remaining two. Composing
symplectic integrators for the subsystems therefore preserves phase-space volume exactly; a
symmetric composition of second-order methods gives a second-order volume-preserving scheme.
"""
function guiding_center_4d_sode(qᵢ=qᵢ, parameters=parameters; tspan=tspan, tstep=Δt, periodic=true)
    SODEProblem(
        (v₁, v₂, v₃, v₄, v₅, v₆),
        tspan, tstep, qᵢ;
        parameters=parameters,
        invariants=(h=hamiltonian,),
        periodicity=guiding_center_4d_periodicity(qᵢ, periodic)
    )
end
