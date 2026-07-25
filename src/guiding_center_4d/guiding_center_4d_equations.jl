
import ElectromagneticFields

using GeometricEquations: ODEProblem, IODEProblem, LODEProblem
using GeometricSolutions: GeometricSolution, DataSeries, TimeSeries


export guiding_center_4d_ode, guiding_center_4d_iode, guiding_center_4d_iode_λ,
    guiding_center_4d_lode,
    guiding_center_4d_dg, guiding_center_4d_formal_lagrangian



function guiding_center_4d_periodicity(::Type{T}, periodic=true) where {T}
    xmin = -Inf * ones(T, 4)
    xmax = +Inf * ones(T, 4)

    if periodic
        xmin[1:3] .= rangemin(xmin[1:3])
        xmax[1:3] .= rangemax(xmax[1:3])
    end

    return (xmin, xmax)
end

guiding_center_4d_periodicity(::AbstractVector{<:AbstractArray{T}}, periodic=true) where {T<:Number} = guiding_center_4d_periodicity(T, periodic)
guiding_center_4d_periodicity(::AbstractArray{T}, periodic=true) where {T<:Number} = guiding_center_4d_periodicity(T, periodic)



function guiding_center_4d_ode(qᵢ=qᵢ, parameters=parameters; tspan=tspan, tstep=Δt, periodic=true)
    ODEProblem(
        guiding_center_4d_v,
        tspan, tstep, qᵢ;
        parameters=parameters,
        invariants=(h=hamiltonian,),
        periodicity=guiding_center_4d_periodicity(qᵢ, periodic)
    )
end


function guiding_center_4d_iode(qᵢ=qᵢ, parameters=parameters; tspan=tspan, tstep=Δt, periodic=true)
    IODEProblem(
        guiding_center_4d_ϑ,
        guiding_center_4d_f,
        guiding_center_4d_g,
        tspan, tstep, qᵢ, guiding_center_4d_pᵢ(tspan[begin], qᵢ);
        parameters=parameters,
        invariants=(h=hamiltonian,),
        periodicity=guiding_center_4d_periodicity(qᵢ, periodic),
        v̄=guiding_center_4d_v
    )
end

function guiding_center_4d_iode_λ(qᵢ=qᵢ, parameters=parameters; tspan=tspan, tstep=Δt, periodic=true)
    IODEProblem(
        guiding_center_4d_ϑ,
        guiding_center_4d_f,
        guiding_center_4d_g,
        tspan, tstep, qᵢ, guiding_center_4d_pᵢ(tspan[begin], qᵢ), guiding_center_4d_λᵢ(tspan[begin], qᵢ, parameters);
        parameters=parameters,
        invariants=(h=hamiltonian,),
        periodicity=guiding_center_4d_periodicity(qᵢ, periodic),
        v̄=guiding_center_4d_v
    )
end

function guiding_center_4d_lode(qᵢ=qᵢ, parameters=parameters; tspan=tspan, tstep=Δt, periodic=true)
    LODEProblem(
        guiding_center_4d_ϑ,
        guiding_center_4d_f,
        guiding_center_4d_g,
        guiding_center_4d_ω, lagrangian,
        tspan, tstep, qᵢ, guiding_center_4d_pᵢ(tspan[begin], qᵢ);
        parameters=parameters,
        invariants=(h=hamiltonian,),
        periodicity=guiding_center_4d_periodicity(qᵢ, periodic),
        v̄=guiding_center_4d_v
    )
end

function guiding_center_4d_dg(qᵢ=qᵢ, parameters=parameters; tspan=tspan, tstep=Δt, periodic=true, κ=0.0)
    # The output array comes first in every GeometricEquations callback, and `g` is called as
    # `g(g, t, q, v, λ, params)`; the κ-form of `g` ignores `v`.
    guiding_center_4d_ϑ_κ(θ, t, q, v, params) = guiding_center_4d_ϑ(θ, t, q, v, params, κ)
    guiding_center_4d_f_κ(f, t, q, v, params) = guiding_center_4d_f(f, t, q, v, params, κ)
    guiding_center_4d_g_κ(g, t, q, v, λ, params) = guiding_center_4d_g(g, t, q, λ, params, κ)

    pᵢ = zero(qᵢ)
    guiding_center_4d_ϑ_κ(pᵢ, tspan[begin], qᵢ, zero(qᵢ), parameters)

    IODEProblem(
        guiding_center_4d_ϑ_κ,
        guiding_center_4d_f_κ,
        guiding_center_4d_g_κ,
        tspan, tstep, qᵢ, pᵢ;
        parameters=parameters,
        invariants=(h=hamiltonian,),
        periodicity=guiding_center_4d_periodicity(qᵢ, periodic),
        v̄=guiding_center_4d_v
    )
end

function guiding_center_4d_formal_lagrangian(qᵢ=qᵢ, parameters=parameters; tspan=tspan, tstep=Δt, periodic=true)
    LODEProblem(
        guiding_center_4d_ϑ,
        guiding_center_4d_f,
        guiding_center_4d_g,
        guiding_center_4d_ω, lagrangian,
        tspan, tstep, qᵢ, guiding_center_4d_pᵢ(tspan[begin], qᵢ);
        parameters=parameters,
        invariants=(h=hamiltonian,),
        periodicity=guiding_center_4d_periodicity(qᵢ, periodic),
        v̄=guiding_center_4d_v
    )
end
