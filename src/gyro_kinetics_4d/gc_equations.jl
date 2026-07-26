
import ElectromagneticFields

using GeometricEquations: ODEProblem, SODEProblem

export odeproblem, sodeproblem


# `GeometricEquations` recognises periodicity only as an `(xmin, xmax)` tuple of arrays; the range
# of the position coordinates comes from the coordinate system, via the `rangemin`/`rangemax`
# injected with the field code. The parallel velocity is not periodic.
#
# Those two take an evaluation point only for uniformity with the other generated field functions.
# The range of a coordinate is a property of the chart, and `minx¹`…`maxx³` are baked in as
# literals, so the argument is discarded; pass the origin rather than `±Inf`.
function guiding_center_4d_periodicity(::Type{T}, periodic=true) where {T}
    xmin = -T(Inf) * ones(T, 4)
    xmax = +T(Inf) * ones(T, 4)

    if periodic
        xmin[1:3] .= rangemin(zeros(T, 3))
        xmax[1:3] .= rangemax(zeros(T, 3))
    end

    return (xmin, xmax)
end

guiding_center_4d_periodicity(::AbstractVector{<:AbstractArray{T}}, periodic=true) where {T<:Number} = guiding_center_4d_periodicity(T, periodic)
guiding_center_4d_periodicity(::AbstractArray{T}, periodic=true) where {T<:Number} = guiding_center_4d_periodicity(T, periodic)


@doc raw"""
    odeproblem(qᵢ; kwargs...)
    odeproblem(ics::NamedTuple; kwargs...)

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

That factor is the phasespace Jacobian, absorbed into the distribution function as
``\tilde{f} = B^{\star}_{\parallel} f``. What it buys is that the right-hand side is then
divergence-free, which is what makes the splitting of [`sodeproblem`](@ref) volume
preserving.

The second form takes the named tuple that every `initial_conditions_*` of this module returns,
whose `params` carry that condition's own magnetic moment `μ`.
"""
function odeproblem(qᵢ=qᵢ; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters(), periodic=true)
    ODEProblem(
        v,
        timespan, timestep, qᵢ;
        parameters=parameters,
        invariants=(h=hamiltonian,),
        periodicity=guiding_center_4d_periodicity(qᵢ, periodic)
    )
end


@doc raw"""
    sodeproblem(qᵢ; kwargs...)
    sodeproblem(ics::NamedTuple; kwargs...)

The same characteristics as [`odeproblem`](@ref), split into the six subsystems

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
symplectic integrators for the subsystems therefore preserves phasespace volume exactly; a
symmetric composition of second-order methods gives a second-order volume-preserving scheme.
"""
function sodeproblem(qᵢ=qᵢ; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters(), periodic=true)
    SODEProblem(
        (v₁, v₂, v₃, v₄, v₅, v₆),
        timespan, timestep, qᵢ;
        parameters=parameters,
        invariants=(h=hamiltonian,),
        periodicity=guiding_center_4d_periodicity(qᵢ, periodic)
    )
end


# `initial_conditions_*` returns `(q = …, params = …)`; take it directly.
for problem in (:odeproblem, :sodeproblem)
    @eval $problem(ics::NamedTuple; kwargs...) = $problem(ics.q; parameters = ics.params, kwargs...)
end
