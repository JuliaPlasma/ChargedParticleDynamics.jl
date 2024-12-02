
import ElectromagneticFields

using GeometricEquations: ODEProblem, SODEProblem

export guiding_center_4d_ode, guiding_center_4d_sode


function guiding_center_4d_periodicity(::Type{T}, periodic=true) where {T}
    period = zeros(T,4)

    if periodic
        period[1:3] .= periodicity(period[1:3])
    end

    return period
end

guiding_center_4d_periodicity(::AbstractVector{<:AbstractArray{T}}, periodic=true) where {T <: Number} = guiding_center_4d_periodicity(T, periodic)
guiding_center_4d_periodicity(::AbstractArray{T}, periodic=true) where {T <: Number} = guiding_center_4d_periodicity(T, periodic)


transform_v(V) = (q̇,t,q,parameters) -> V(q̇, t, transform_q̃_to_q(t, q, parameters), parameters)

function guiding_center_4d_ode(qᵢ = qᵢ, parameters = parameters; tspan = tspan, tstep = Δt, periodic = true)
    ODEProblem(
        transform_v(v),
        tspan, tstep, qᵢ; 
        parameters = parameters,
        invariants = (h = hamiltonian,),
        periodicity = guiding_center_4d_periodicity(qᵢ, periodic)
    )
end

function guiding_center_4d_sode(qᵢ = qᵢ, parameters = parameters; tspan = tspan, tstep = Δt, periodic = true)
    SODEProblem(
        Tuple(transform_v(V) for V in (v₁, v₂, v₃, v₄, v₅, v₆)),
        tspan, tstep, qᵢ;
        parameters = parameters,
        invariants = (h = hamiltonian,),
        periodicity = guiding_center_4d_periodicity(qᵢ, periodic)
    )
end
