
import ElectromagneticFields

using GeometricEquations: ODEProblem, IODEProblem, LODEProblem
using GeometricSolutions: GeometricSolution, DataSeries, TimeSeries


export guiding_center_4d_ode, guiding_center_4d_iode, guiding_center_4d_iode_λ,
       guiding_center_4d_iode_dec128, guiding_center_4d_lode,
       guiding_center_4d_dg, guiding_center_4d_formal_lagrangian



function guiding_center_4d_periodicity(::Type{T}, periodic=true) where {T}
    period = zeros(T,4)

    if periodic
        period[1:3] .= periodicity(period[1:3])
    end

    return period
end

guiding_center_4d_periodicity(::AbstractVector{<:AbstractArray{T}}, periodic=true) where {T <: Number} = guiding_center_4d_periodicity(T, periodic)
guiding_center_4d_periodicity(::AbstractArray{T}, periodic=true) where {T <: Number} = guiding_center_4d_periodicity(T, periodic)


function guiding_center_4d_ode(qᵢ = qᵢ, parameters = parameters; tspan = tspan, tstep = Δt, periodic = true)
    ODEProblem(
        guiding_center_4d_v,
        tspan, tstep, qᵢ;
        parameters = parameters,
        invariants = (h = hamiltonian,),
        periodicity = guiding_center_4d_periodicity(qᵢ, periodic)
    )
end


function guiding_center_4d_iode(qᵢ = qᵢ, parameters = parameters; tspan = tspan, tstep = Δt, periodic = true)
    IODEProblem(
        guiding_center_4d_ϑ,
        guiding_center_4d_f,
        guiding_center_4d_g,
        tspan, tstep, qᵢ, guiding_center_4d_pᵢ(tspan[begin], qᵢ);
        parameters = parameters,
        invariants = (h = hamiltonian,),
        periodicity = guiding_center_4d_periodicity(qᵢ, periodic),
        v̄ = guiding_center_4d_v
    )
end

function guiding_center_4d_iode_dec128(qᵢ = qᵢ, parameters = parameters; tspan = tspan, tstep = Δt, periodic = true)
    IODEProblem(
        guiding_center_4d_ϑ,
        guiding_center_4d_f,
        guiding_center_4d_g,
        tspan, tstep, Dec128.(qᵢ), guiding_center_4d_pᵢ(tspan[begin], Dec128.(qᵢ));
        parameters = parameters,
        invariants = (h = hamiltonian,),
        periodicity = guiding_center_4d_periodicity(qᵢ, periodic),
        v̄ = guiding_center_4d_v
    )
end

function guiding_center_4d_iode_λ(qᵢ = qᵢ, parameters = parameters; tspan = tspan, tstep = Δt, periodic = true)
    IODEProblem(
        guiding_center_4d_ϑ,
        guiding_center_4d_f,
        guiding_center_4d_g,
        tspan, tstep, qᵢ, guiding_center_4d_pᵢ(tspan[begin], qᵢ), guiding_center_4d_λ₀(qᵢ);
        parameters = parameters,
        invariants = (h = hamiltonian,),
        periodicity = guiding_center_4d_periodicity(qᵢ, periodic),
        v̄ = guiding_center_4d_v
    )
end

function guiding_center_4d_lode(qᵢ = qᵢ, parameters = parameters; tspan = tspan, tstep = Δt, periodic = true)
    LODEProblem(
        guiding_center_4d_ϑ,
        guiding_center_4d_f,
        guiding_center_4d_g,
        guiding_center_4d_ω, lagrangian,
        tspan, tstep, qᵢ, guiding_center_4d_pᵢ(tspan[begin], qᵢ);
        parameters = parameters,
        invariants = (h = hamiltonian,),
        periodicity = guiding_center_4d_periodicity(qᵢ, periodic),
        v̄ = guiding_center_4d_v
    )
end

function guiding_center_4d_dg(qᵢ = qᵢ, parameters = parameters; tspan = tspan, tstep = Δt, periodic = true, κ = 0.0)
    guiding_center_4d_ϑ_κ(t, q, v, p, params) = guiding_center_4d_ϑ(t, q, v, p, params, κ)
    guiding_center_4d_f_κ(t, q, v, f, params) = guiding_center_4d_f(t, q, v, f, params, κ)
    guiding_center_4d_g_κ(t, q, λ, g, params) = guiding_center_4d_g(t, q, λ, g, params, κ)

    IODEProblem(
        guiding_center_4d_ϑ_κ,
        guiding_center_4d_f_κ,
        guiding_center_4d_g_κ,
        tspan, tstep, qᵢ, qᵢ;
        parameters = parameters,
        invariants = (h = hamiltonian,),
        periodicity = guiding_center_4d_periodicity(qᵢ, periodic),
        v̄ = guiding_center_4d_v
    )
end

function guiding_center_4d_formal_lagrangian(qᵢ = qᵢ, parameters = parameters; tspan = tspan, tstep = Δt, periodic = true)
    LODEProblem(
        guiding_center_4d_ϑ,
        guiding_center_4d_f,
        guiding_center_4d_g,
        guiding_center_4d_ω, lagrangian,
        tspan, tstep, qᵢ, guiding_center_4d_pᵢ(tspan[begin], qᵢ);
        parameters = parameters,
        invariants = (h = hamiltonian,),
        periodicity = guiding_center_4d_periodicity(qᵢ, periodic),
        v̄ = guiding_center_4d_v
    )
end
