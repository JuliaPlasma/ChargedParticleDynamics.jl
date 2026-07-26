using LinearAlgebra
using Parameters

import GeometricEquations: HODEProblem
import GeometricSolutions: GeometricSolution, DataSeries, TimeSeries

export hamiltonian, hamiltonian_canonical
export hodeproblem, hodeproblem_canonical


@views ϑ₁(t, Q) = A₁(t, Q[1:3]) + Q[4] * b₁(t, Q[1:3])
@views ϑ₂(t, Q) = A₂(t, Q[1:3]) + Q[4] * b₂(t, Q[1:3])
@views ϑ₃(t, Q) = A₃(t, Q[1:3]) + Q[4] * b₃(t, Q[1:3])

v₁(t, q, p) = p[1] - A₁(t, q)
v₂(t, q, p) = p[2] - A₂(t, q)
v₃(t, q, p) = p[3] - A₃(t, q)

dv₁dq₁(t, q, p) = -dA₁dx₁(t, q)
dv₁dq₂(t, q, p) = -dA₁dx₂(t, q)
dv₁dq₃(t, q, p) = -dA₁dx₃(t, q)

dv₂dq₁(t, q, p) = -dA₂dx₁(t, q)
dv₂dq₂(t, q, p) = -dA₂dx₂(t, q)
dv₂dq₃(t, q, p) = -dA₂dx₃(t, q)

dv₃dq₁(t, q, p) = -dA₃dx₁(t, q)
dv₃dq₂(t, q, p) = -dA₃dx₂(t, q)
dv₃dq₃(t, q, p) = -dA₃dx₃(t, q)

dv₁dp₁(t, q, p) = one(eltype(p))
dv₁dp₂(t, q, p) = zero(eltype(p))
dv₁dp₃(t, q, p) = zero(eltype(p))

dv₂dp₁(t, q, p) = zero(eltype(p))
dv₂dp₂(t, q, p) = one(eltype(p))
dv₂dp₃(t, q, p) = zero(eltype(p))

dv₃dp₁(t, q, p) = zero(eltype(p))
dv₃dp₂(t, q, p) = zero(eltype(p))
dv₃dp₃(t, q, p) = one(eltype(p))


u(t, q, p) = v₁(t, q, p) * g¹¹(t, q) * b₁(t, q) + v₂(t, q, p) * g²²(t, q) * b₂(t, q) + v₃(t, q, p) * g³³(t, q) * b₃(t, q)

# dudq₁(t, q, p) = v₁(t, q, p) * db₁dx₁(t, q) + v₂(t, q, p) * db₂dx₁(t, q) + v₃(t, q, p) * db₃dx₁(t, q) +
#                  dv₁dq₁(t, q, p) * b₁(t, q) + dv₂dq₁(t, q, p) * b₂(t, q) + dv₃dq₁(t, q, p) * b₃(t, q)
# dudq₂(t, q, p) = v₁(t, q, p) * db₁dx₂(t, q) + v₂(t, q, p) * db₂dx₂(t, q) + v₃(t, q, p) * db₃dx₂(t, q) +
#                  dv₁dq₂(t, q, p) * b₁(t, q) + dv₂dq₂(t, q, p) * b₂(t, q) + dv₃dq₂(t, q, p) * b₃(t, q)
# dudq₃(t, q, p) = v₁(t, q, p) * db₁dx₃(t, q) + v₂(t, q, p) * db₂dx₃(t, q) + v₃(t, q, p) * db₃dx₃(t, q) +
#                  dv₁dq₃(t, q, p) * b₁(t, q) + dv₂dq₃(t, q, p) * b₂(t, q) + dv₃dq₃(t, q, p) * b₃(t, q)

# dudp₁(t, q, p) = b₁(t, q)
# dudp₂(t, q, p) = b₂(t, q)
# dudp₃(t, q, p) = b₃(t, q)


function initial_momentum(tᵢ, Qᵢ::AbstractArray{T}) where {T<:Number}
    pᵢ = zeros(T, 3)
    pᵢ[1] = ϑ₁(tᵢ, Qᵢ)
    pᵢ[2] = ϑ₂(tᵢ, Qᵢ)
    pᵢ[3] = ϑ₃(tᵢ, Qᵢ)
    return pᵢ
end

function fix_initial_momentum(tᵢ, qᵢ::AbstractArray{T}, pᵢ::AbstractArray{T}) where {T<:Number}
    # `u` raises the index with the inverse metric; spelling the contraction out here without it
    # made this disagree with `u(t, q, p)` above in every curvilinear coordinate system.
    initial_momentum(tᵢ, [qᵢ..., u(tᵢ, qᵢ, pᵢ)])
end

function initial_conditions(tᵢ, Qᵢ::AbstractArray{T}) where {T<:Number}
    qᵢ = Qᵢ[1:3]
    pᵢ = initial_momentum(tᵢ, Qᵢ)
    (q=qᵢ, p=pᵢ)
end


# `rangemin`/`rangemax` take an evaluation point only for uniformity with the other generated field
# functions (`B(t,x)`, `A₁(t,x)`, …). The range of a coordinate is a property of the chart, not a
# field sampled at a point: `ElectromagneticFields` bakes `minx¹`…`maxx³` in as literals when it
# injects the field code, so the argument is discarded. Pass the origin rather than the `±Inf` that
# `xmin`/`xmax` are initialised to, which read as if the range were being asked for at the point at
# infinity.
function guiding_center_3d_periodicity(::Type{T}, periodic=true) where {T}
    xmin = -Inf * ones(T, 3)
    xmax = +Inf * ones(T, 3)

    if periodic
        xmin[1:3] .= rangemin(zeros(T, 3))
        xmax[1:3] .= rangemax(zeros(T, 3))
    end

    return (xmin, xmax)
end

guiding_center_3d_periodicity(::AbstractVector{<:AbstractArray{T}}, periodic=true) where {T<:Number} = guiding_center_3d_periodicity(T, periodic)
guiding_center_3d_periodicity(::AbstractArray{T}, periodic=true) where {T<:Number} = guiding_center_3d_periodicity(T, periodic)


hamiltonian(t, q, p, params) = g¹¹(t, q) * v₁(t, q, p)^2 / 2 + g²²(t, q) * v₂(t, q, p)^2 / 2 + g³³(t, q) * v₃(t, q, p)^2 / 2 + params.μ * B(t, q) + φ(t, q)
# The same Hamiltonian with the parallel velocity in place of the full kinetic energy. The two
# agree wherever the constraint `v × b = 0` holds, since `v` is then purely parallel and `|b| = 1`,
# so their difference along a trajectory is a measure of the drift off the constraint manifold —
# which is what `scripts/guiding_center_3d_*.jl` plot it for. It omitted `φ`, unlike `hamiltonian`
# beside it, which made the two incomparable for any equilibrium with a potential.
hamiltonian_u(t, q, p, params) = u(t, q, p)^2 / 2 + params.μ * B(t, q) + φ(t, q)
hamiltonian_canonical(t, q, p, params) = hamiltonian(t, q, p, params) + λ₁(t, q, p, params) * g₁(t, q, p) +
                                         λ₂(t, q, p, params) * g₂(t, q, p)


dHdq₁(t, q, p, params) = v₁(t, q, p) * g¹¹(t, q) * dv₁dq₁(t, q, p) +
                         v₂(t, q, p) * g²²(t, q) * dv₂dq₁(t, q, p) +
                         v₃(t, q, p) * g³³(t, q) * dv₃dq₁(t, q, p) +
                         v₁(t, q, p) * dg¹¹dx₁(t, q) * v₁(t, q, p) / 2 +
                         v₂(t, q, p) * dg²²dx₁(t, q) * v₂(t, q, p) / 2 +
                         v₃(t, q, p) * dg³³dx₁(t, q) * v₃(t, q, p) / 2 +
                         params.μ * dBdx₁(t, q) -
                         E₁(t, q)

dHdq₂(t, q, p, params) = v₁(t, q, p) * g¹¹(t, q) * dv₁dq₂(t, q, p) +
                         v₂(t, q, p) * g²²(t, q) * dv₂dq₂(t, q, p) +
                         v₃(t, q, p) * g³³(t, q) * dv₃dq₂(t, q, p) +
                         v₁(t, q, p) * dg¹¹dx₂(t, q) * v₁(t, q, p) / 2 +
                         v₂(t, q, p) * dg²²dx₂(t, q) * v₂(t, q, p) / 2 +
                         v₃(t, q, p) * dg³³dx₂(t, q) * v₃(t, q, p) / 2 +
                         params.μ * dBdx₂(t, q) -
                         E₂(t, q)

dHdq₃(t, q, p, params) = v₁(t, q, p) * g¹¹(t, q) * dv₁dq₃(t, q, p) +
                         v₂(t, q, p) * g²²(t, q) * dv₂dq₃(t, q, p) +
                         v₃(t, q, p) * g³³(t, q) * dv₃dq₃(t, q, p) +
                         v₁(t, q, p) * dg¹¹dx₃(t, q) * v₁(t, q, p) / 2 +
                         v₂(t, q, p) * dg²²dx₃(t, q) * v₂(t, q, p) / 2 +
                         v₃(t, q, p) * dg³³dx₃(t, q) * v₃(t, q, p) / 2 +
                         params.μ * dBdx₃(t, q) -
                         E₃(t, q)

dHdp₁(t, q, p, params) = v₁(t, q, p) * g¹¹(t, q) * dv₁dp₁(t, q, p) +
                         v₂(t, q, p) * g²²(t, q) * dv₂dp₁(t, q, p) +
                         v₃(t, q, p) * g³³(t, q) * dv₃dp₁(t, q, p)

dHdp₂(t, q, p, params) = v₁(t, q, p) * g¹¹(t, q) * dv₁dp₂(t, q, p) +
                         v₂(t, q, p) * g²²(t, q) * dv₂dp₂(t, q, p) +
                         v₃(t, q, p) * g³³(t, q) * dv₃dp₂(t, q, p)

dHdp₃(t, q, p, params) = v₁(t, q, p) * g¹¹(t, q) * dv₁dp₃(t, q, p) +
                         v₂(t, q, p) * g²²(t, q) * dv₂dp₃(t, q, p) +
                         v₃(t, q, p) * g³³(t, q) * dv₃dp₃(t, q, p)


# function dHdq(dH, t, q, p, params)
#     dH[1] = dHdq₁(t, q, p, params)
#     dH[2] = dHdq₂(t, q, p, params)
#     dH[3] = dHdq₃(t, q, p, params)
#     nothing
# end

# function dHdp(dH, t, q, p, params)
#     dH[1] = dHdp₁(t, q, p, params)
#     dH[2] = dHdp₂(t, q, p, params)
#     dH[3] = dHdp₃(t, q, p, params)
#     nothing
# end

g₁(t, q, p) = b₁(t, q) * v₂(t, q, p) - b₂(t, q) * v₁(t, q, p)
g₂(t, q, p) = b₁(t, q) * v₃(t, q, p) - b₃(t, q) * v₁(t, q, p)
# g₁(t, q, p) = b₁(t, q) * v₃(t, q, p) - b₃(t, q) * v₁(t, q, p)
# g₂(t, q, p) = b₂(t, q) * v₃(t, q, p) - b₃(t, q) * v₂(t, q, p)

g₁(t, q, p, params) = g₁(t, q, p)
g₂(t, q, p, params) = g₂(t, q, p)

dg₁dq₁(t, q, p) = (db₁dx₁(t, q) * v₂(t, q, p) + b₁(t, q) * dv₂dq₁(t, q, p)) -
                  (db₂dx₁(t, q) * v₁(t, q, p) + b₂(t, q) * dv₁dq₁(t, q, p))
dg₁dq₂(t, q, p) = (db₁dx₂(t, q) * v₂(t, q, p) + b₁(t, q) * dv₂dq₂(t, q, p)) -
                  (db₂dx₂(t, q) * v₁(t, q, p) + b₂(t, q) * dv₁dq₂(t, q, p))
dg₁dq₃(t, q, p) = (db₁dx₃(t, q) * v₂(t, q, p) + b₁(t, q) * dv₂dq₃(t, q, p)) -
                  (db₂dx₃(t, q) * v₁(t, q, p) + b₂(t, q) * dv₁dq₃(t, q, p))

dg₂dq₁(t, q, p) = (db₁dx₁(t, q) * v₃(t, q, p) + b₁(t, q) * dv₃dq₁(t, q, p)) -
                  (db₃dx₁(t, q) * v₁(t, q, p) + b₃(t, q) * dv₁dq₁(t, q, p))
dg₂dq₂(t, q, p) = (db₁dx₂(t, q) * v₃(t, q, p) + b₁(t, q) * dv₃dq₂(t, q, p)) -
                  (db₃dx₂(t, q) * v₁(t, q, p) + b₃(t, q) * dv₁dq₂(t, q, p))
dg₂dq₃(t, q, p) = (db₁dx₃(t, q) * v₃(t, q, p) + b₁(t, q) * dv₃dq₃(t, q, p)) -
                  (db₃dx₃(t, q) * v₁(t, q, p) + b₃(t, q) * dv₁dq₃(t, q, p))

dg₁dp₁(t, q, p) = b₁(t, q) * dv₂dp₁(t, q, p) - b₂(t, q) * dv₁dp₁(t, q, p)
dg₁dp₂(t, q, p) = b₁(t, q) * dv₂dp₂(t, q, p) - b₂(t, q) * dv₁dp₂(t, q, p)
dg₁dp₃(t, q, p) = b₁(t, q) * dv₂dp₃(t, q, p) - b₂(t, q) * dv₁dp₃(t, q, p)

dg₂dp₁(t, q, p) = b₁(t, q) * dv₃dp₁(t, q, p) - b₃(t, q) * dv₁dp₁(t, q, p)
dg₂dp₂(t, q, p) = b₁(t, q) * dv₃dp₂(t, q, p) - b₃(t, q) * dv₁dp₂(t, q, p)
dg₂dp₃(t, q, p) = b₁(t, q) * dv₃dp₃(t, q, p) - b₃(t, q) * dv₁dp₃(t, q, p)


# dg₁dq₁(t, q, p) = (db₁dx₁(t, q) * v₃(t, q, p) + b₁(t, q) * dv₃dq₁(t, q, p)) -
#                   (db₃dx₁(t, q) * v₁(t, q, p) + b₃(t, q) * dv₁dq₁(t, q, p))
# dg₁dq₂(t, q, p) = (db₁dx₂(t, q) * v₃(t, q, p) + b₁(t, q) * dv₃dq₂(t, q, p)) -
#                   (db₃dx₂(t, q) * v₁(t, q, p) + b₃(t, q) * dv₁dq₂(t, q, p))
# dg₁dq₃(t, q, p) = (db₁dx₃(t, q) * v₃(t, q, p) + b₁(t, q) * dv₃dq₃(t, q, p)) -
#                   (db₃dx₃(t, q) * v₁(t, q, p) + b₃(t, q) * dv₁dq₃(t, q, p))

# dg₂dq₁(t, q, p) = (db₂dx₁(t, q) * v₃(t, q, p) + b₂(t, q) * dv₃dq₁(t, q, p)) -
#                   (db₃dx₁(t, q) * v₂(t, q, p) + b₃(t, q) * dv₂dq₁(t, q, p))
# dg₂dq₂(t, q, p) = (db₂dx₂(t, q) * v₃(t, q, p) + b₂(t, q) * dv₃dq₂(t, q, p)) -
#                   (db₃dx₂(t, q) * v₂(t, q, p) + b₃(t, q) * dv₂dq₂(t, q, p))
# dg₂dq₃(t, q, p) = (db₂dx₃(t, q) * v₃(t, q, p) + b₂(t, q) * dv₃dq₃(t, q, p)) -
#                   (db₃dx₃(t, q) * v₂(t, q, p) + b₃(t, q) * dv₂dq₃(t, q, p))

# dg₁dp₁(t, q, p) = b₁(t, q) * dv₃dp₁(t, q, p) - b₃(t, q) * dv₁dp₁(t, q, p)
# dg₁dp₂(t, q, p) = b₁(t, q) * dv₃dp₂(t, q, p) - b₃(t, q) * dv₁dp₂(t, q, p)
# dg₁dp₃(t, q, p) = b₁(t, q) * dv₃dp₃(t, q, p) - b₃(t, q) * dv₁dp₃(t, q, p)

# dg₂dp₁(t, q, p) = b₂(t, q) * dv₃dp₁(t, q, p) - b₃(t, q) * dv₂dp₁(t, q, p)
# dg₂dp₂(t, q, p) = b₂(t, q) * dv₃dp₂(t, q, p) - b₃(t, q) * dv₂dp₂(t, q, p)
# dg₂dp₃(t, q, p) = b₂(t, q) * dv₃dp₃(t, q, p) - b₃(t, q) * dv₂dp₃(t, q, p)


λₒ(t, q, p) = (
    dg₁dq₁(t, q, p) * dg₂dp₁(t, q, p) +
    dg₁dq₂(t, q, p) * dg₂dp₂(t, q, p) +
    dg₁dq₃(t, q, p) * dg₂dp₃(t, q, p) -
    dg₁dp₁(t, q, p) * dg₂dq₁(t, q, p) -
    dg₁dp₂(t, q, p) * dg₂dq₂(t, q, p) -
    dg₁dp₃(t, q, p) * dg₂dq₃(t, q, p)
)

λ₁(t, q, p, params) = +(
    dg₂dq₁(t, q, p) * dHdp₁(t, q, p, params) +
    dg₂dq₂(t, q, p) * dHdp₂(t, q, p, params) +
    dg₂dq₃(t, q, p) * dHdp₃(t, q, p, params) -
    dg₂dp₁(t, q, p) * dHdq₁(t, q, p, params) -
    dg₂dp₂(t, q, p) * dHdq₂(t, q, p, params) -
    dg₂dp₃(t, q, p) * dHdq₃(t, q, p, params)
) / λₒ(t, q, p)

λ₂(t, q, p, params) = -(
    dg₁dq₁(t, q, p) * dHdp₁(t, q, p, params) +
    dg₁dq₂(t, q, p) * dHdp₂(t, q, p, params) +
    dg₁dq₃(t, q, p) * dHdp₃(t, q, p, params) -
    dg₁dp₁(t, q, p) * dHdq₁(t, q, p, params) -
    dg₁dp₂(t, q, p) * dHdq₂(t, q, p, params) -
    dg₁dp₃(t, q, p) * dHdq₃(t, q, p, params)
) / λₒ(t, q, p)



# function guiding_center_3d_ode_v(v, t, x, params)
#     q = @views x[1:3]
#     p = @views x[4:6]

#     v[1] = dHdp₁(t, q, p, params) + λ₁(t, q, p, params) * dg₁dp₁(t, q, p) + λ₂(t, q, p, params) * dg₂dp₁(t, q, p)
#     v[2] = dHdp₂(t, q, p, params) + λ₁(t, q, p, params) * dg₁dp₂(t, q, p) + λ₂(t, q, p, params) * dg₂dp₂(t, q, p)
#     v[3] = dHdp₃(t, q, p, params) + λ₁(t, q, p, params) * dg₁dp₃(t, q, p) + λ₂(t, q, p, params) * dg₂dp₃(t, q, p)
#     v[4] = -dHdq₁(t, q, p, params) - λ₁(t, q, p, params) * dg₁dq₁(t, q, p) - λ₂(t, q, p, params) * dg₂dq₁(t, q, p)
#     v[5] = -dHdq₂(t, q, p, params) - λ₁(t, q, p, params) * dg₁dq₂(t, q, p) - λ₂(t, q, p, params) * dg₂dq₂(t, q, p)
#     v[6] = -dHdq₃(t, q, p, params) - λ₁(t, q, p, params) * dg₁dq₃(t, q, p) - λ₂(t, q, p, params) * dg₂dq₃(t, q, p)

#     nothing
# end

# function ode(x₀, parameters; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, periodic=false)
#     # println("3D Guiding Center model initial constraints g₁ = $(g₁(t₀, q₀, p₀)) and g₂ = $(g₂(t₀, q₀, p₀))")
#     ODEProblem(
#         guiding_center_3d_ode_v,
#         timespan, timestep, x₀;
#         parameters=parameters,
#         periodicity=guiding_center_3d_periodicity(q₀, periodic))
# end


function guiding_center_3d_v(v, t, q, p, params)
    v[1] = dHdp₁(t, q, p, params) + λ₁(t, q, p, params) * dg₁dp₁(t, q, p) + λ₂(t, q, p, params) * dg₂dp₁(t, q, p)
    v[2] = dHdp₂(t, q, p, params) + λ₁(t, q, p, params) * dg₁dp₂(t, q, p) + λ₂(t, q, p, params) * dg₂dp₂(t, q, p)
    v[3] = dHdp₃(t, q, p, params) + λ₁(t, q, p, params) * dg₁dp₃(t, q, p) + λ₂(t, q, p, params) * dg₂dp₃(t, q, p)
    nothing
end

function guiding_center_3d_f(f, t, q, p, params)
    f[1] = -dHdq₁(t, q, p, params) - λ₁(t, q, p, params) * dg₁dq₁(t, q, p) - λ₂(t, q, p, params) * dg₂dq₁(t, q, p)
    f[2] = -dHdq₂(t, q, p, params) - λ₁(t, q, p, params) * dg₁dq₂(t, q, p) - λ₂(t, q, p, params) * dg₂dq₂(t, q, p)
    f[3] = -dHdq₃(t, q, p, params) - λ₁(t, q, p, params) * dg₁dq₃(t, q, p) - λ₂(t, q, p, params) * dg₂dq₃(t, q, p)
    nothing
end

"""
    hodeproblem(q₀, p₀; kwargs...)
    hodeproblem(x₀ = qᵢ; kwargs...)
    hodeproblem(ics::NamedTuple; kwargs...)

The constrained canonical guiding centre system as an `HODEProblem` in the position and its
conjugate momentum.

The first form takes the position and momentum directly. The second takes the four-component state
``(x, u)`` and recovers the momentum from it through `initial_conditions`; this is the form
the module constant `qᵢ` is in. The third takes the named tuple that every `initial_conditions_*`
returns, so `hodeproblem(initial_conditions_barely_passing())` carries that condition's own `μ`.
"""
function hodeproblem(q₀::AbstractVector, p₀::AbstractVector; timespan = DEFAULT_TIMESPAN,
                     timestep = DEFAULT_TIMESTEP, parameters = default_parameters(), periodic = true)
    HODEProblem(
        guiding_center_3d_v,
        guiding_center_3d_f,
        hamiltonian,
        timespan, timestep, q₀, p₀;
        parameters=parameters,
        periodicity=guiding_center_3d_periodicity(q₀, periodic))
end

function hodeproblem(x₀::AbstractVector = qᵢ; timespan = DEFAULT_TIMESPAN, kwargs...)
    ics = initial_conditions(timespan[begin], x₀)
    hodeproblem(ics.q, ics.p; timespan = timespan, kwargs...)
end

hodeproblem(ics::NamedTuple; kwargs...) = hodeproblem(ics.q, ics.p; parameters = ics.params, kwargs...)
