
import ElectromagneticFields

using GeometricEquations: ODEProblem, IODEProblem, LODEProblem
using GeometricSolutions: GeometricSolution, DataSeries, TimeSeries


export odeproblem, iodeproblem, iodeproblem_λ,
    lodeproblem,
    iodeproblem_dg, lodeproblem_formal_lagrangian



# `rangemin`/`rangemax` take an evaluation point only for uniformity with the other generated field
# functions; the range of a coordinate is a property of the chart, and `ElectromagneticFields` bakes
# `minx¹`…`maxx³` in as literals, so the argument is discarded. Pass the origin rather than the
# `±Inf` that `xmin`/`xmax` are initialised to. The parallel velocity is not periodic.
function guiding_center_4d_periodicity(::Type{T}, periodic=true) where {T}
    xmin = -Inf * ones(T, 4)
    xmax = +Inf * ones(T, 4)

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

The 4D guiding centre dynamics as an explicit `ODEProblem` in the state ``q = (X, u)``, the
guiding centre position and the parallel velocity. The vector field is obtained from
``\Omega \dot{q} = - \nabla H``.

The second form takes the named tuple that every `initial_conditions_*` of this module returns,
whose `params` carry that condition's own magnetic moment `μ`.
"""
function odeproblem(qᵢ=qᵢ; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters(), periodic=true)
    ODEProblem(
        guiding_center_4d_v,
        timespan, timestep, qᵢ;
        parameters=parameters,
        invariants=(h=hamiltonian,),
        periodicity=guiding_center_4d_periodicity(qᵢ, periodic)
    )
end


@doc raw"""
    iodeproblem(qᵢ; kwargs...)
    iodeproblem(ics::NamedTuple; kwargs...)

The same dynamics as [`odeproblem`](@ref) in implicit form, as an `IODEProblem` built from the
one-form ``\vartheta = A + u b``, so that the variational and projection integrators can be applied
to it. The initial momentum is ``\vartheta(q_{i})``.
"""
function iodeproblem(qᵢ=qᵢ; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters(), periodic=true)
    IODEProblem(
        guiding_center_4d_ϑ,
        guiding_center_4d_f,
        guiding_center_4d_g,
        timespan, timestep, qᵢ, guiding_center_4d_pᵢ(timespan[begin], qᵢ);
        parameters=parameters,
        invariants=(h=hamiltonian,),
        periodicity=guiding_center_4d_periodicity(qᵢ, periodic),
        v̄=guiding_center_4d_v
    )
end

"""
    iodeproblem_λ(qᵢ; kwargs...)

[`iodeproblem`](@ref) with the Lagrange multiplier initialised explicitly, for the integrators that
need a starting value for it rather than taking zero.
"""
function iodeproblem_λ(qᵢ=qᵢ; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters(), periodic=true)
    IODEProblem(
        guiding_center_4d_ϑ,
        guiding_center_4d_f,
        guiding_center_4d_g,
        timespan, timestep, qᵢ, guiding_center_4d_pᵢ(timespan[begin], qᵢ), guiding_center_4d_λᵢ(timespan[begin], qᵢ, parameters);
        parameters=parameters,
        invariants=(h=hamiltonian,),
        periodicity=guiding_center_4d_periodicity(qᵢ, periodic),
        v̄=guiding_center_4d_v
    )
end

@doc raw"""
    lodeproblem(qᵢ; kwargs...)
    lodeproblem(ics::NamedTuple; kwargs...)

The 4D guiding centre dynamics as an `LODEProblem`, from the degenerate phasespace Lagrangian
``L = \vartheta \cdot \dot{q} - H`` together with the two-form ``\omega``. Same equations as
[`iodeproblem`](@ref); the Lagrangian and the two-form are carried in addition, for the variational
integrators that want them.
"""
function lodeproblem(qᵢ=qᵢ; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters(), periodic=true)
    LODEProblem(
        guiding_center_4d_ϑ,
        guiding_center_4d_f,
        guiding_center_4d_g,
        guiding_center_4d_ω, lagrangian,
        timespan, timestep, qᵢ, guiding_center_4d_pᵢ(timespan[begin], qᵢ);
        parameters=parameters,
        invariants=(h=hamiltonian,),
        periodicity=guiding_center_4d_periodicity(qᵢ, periodic),
        v̄=guiding_center_4d_v
    )
end

function iodeproblem_dg(qᵢ=qᵢ; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters(), periodic=true, κ=0.0)
    # The output array comes first in every GeometricEquations callback, and `g` is called as
    # `g(g, t, q, v, λ, params)`; the κ-form of `g` ignores `v`.
    guiding_center_4d_ϑ_κ(θ, t, q, v, params) = guiding_center_4d_ϑ(θ, t, q, v, params, κ)
    guiding_center_4d_f_κ(f, t, q, v, params) = guiding_center_4d_f(f, t, q, v, params, κ)
    guiding_center_4d_g_κ(g, t, q, v, λ, params) = guiding_center_4d_g(g, t, q, λ, params, κ)

    pᵢ = zero(qᵢ)
    guiding_center_4d_ϑ_κ(pᵢ, timespan[begin], qᵢ, zero(qᵢ), parameters)

    IODEProblem(
        guiding_center_4d_ϑ_κ,
        guiding_center_4d_f_κ,
        guiding_center_4d_g_κ,
        timespan, timestep, qᵢ, pᵢ;
        parameters=parameters,
        invariants=(h=hamiltonian,),
        periodicity=guiding_center_4d_periodicity(qᵢ, periodic),
        v̄=guiding_center_4d_v
    )
end

@doc raw"""
    lodeproblem_formal_lagrangian(qᵢ; kwargs...)

!!! warning "Not yet the formal Lagrangian"
    This is currently a verbatim copy of [`lodeproblem`](@ref) — the same `LODEProblem` from the
    same `ϑ`, `f`, `g`, `ω` and `lagrangian`. It should instead be the formal Lagrangian of the
    equations of motion in noncanonical Hamiltonian form: writing them as
    ``F_{i}(z, \dot{z}) = \Omega_{ij}(z) \dot{z}^{j} + \partial_{i} H(z) = 0``, the formal
    Lagrangian is ``L(z, y, \dot{z}) = y^{i} F_{i}(z, \dot{z})`` on the doubled state ``(z, y)``,
    whose Euler-Lagrange equations return the system on variation of ``y`` and its adjoint on
    variation of ``z``. That is an eight-dimensional problem, not the four-dimensional one built
    here. See `TODO.md`.
"""
function lodeproblem_formal_lagrangian(qᵢ=qᵢ; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters(), periodic=true)
    LODEProblem(
        guiding_center_4d_ϑ,
        guiding_center_4d_f,
        guiding_center_4d_g,
        guiding_center_4d_ω, lagrangian,
        timespan, timestep, qᵢ, guiding_center_4d_pᵢ(timespan[begin], qᵢ);
        parameters=parameters,
        invariants=(h=hamiltonian,),
        periodicity=guiding_center_4d_periodicity(qᵢ, periodic),
        v̄=guiding_center_4d_v
    )
end


# Every `initial_conditions_*` returns `(q = …, params = …)`, so each constructor takes that named
# tuple directly: `odeproblem(initial_conditions_deeply_passing())`. The parameters travel with the
# initial condition because `μ` differs between them.
for problem in (:odeproblem, :iodeproblem, :iodeproblem_λ, :lodeproblem,
                :iodeproblem_dg, :lodeproblem_formal_lagrangian)
    @eval $problem(ics::NamedTuple; kwargs...) = $problem(ics.q; parameters = ics.params, kwargs...)
end
