using LinearAlgebra
using Parameters

import GeometricEquations: HODEProblem
import GeometricSolutions: GeometricSolution, DataSeries, TimeSeries

export hamiltonian
export hode


ϑ₁(t, Q) = A₁(t, Q[1:3]) + Q[4] * b₁(t, Q[1:3])
ϑ₂(t, Q) = A₂(t, Q[1:3]) + Q[4] * b₂(t, Q[1:3])
ϑ₃(t, Q) = A₃(t, Q[1:3]) + Q[4] * b₃(t, Q[1:3])

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

u(t, q, p) = v₁(t, q, p) * b₁(t, q) + v₂(t, q, p) * b₂(t, q) + v₃(t, q, p) * b₃(t, q)

dudq₁(t, q, p) = v₁(t, q, p) * db₁dx₁(t, q) + v₂(t, q, p) * db₂dx₁(t, q) + v₃(t, q, p) * db₃dx₁(t, q) +
                 dv₁dq₁(t, q, p) * b₁(t, q) + dv₂dq₁(t, q, p) * b₂(t, q) + dv₃dq₁(t, q, p) * b₃(t, q)
dudq₂(t, q, p) = v₁(t, q, p) * db₁dx₂(t, q) + v₂(t, q, p) * db₂dx₂(t, q) + v₃(t, q, p) * db₃dx₂(t, q) +
                 dv₁dq₂(t, q, p) * b₁(t, q) + dv₂dq₂(t, q, p) * b₂(t, q) + dv₃dq₂(t, q, p) * b₃(t, q)
dudq₃(t, q, p) = v₁(t, q, p) * db₁dx₃(t, q) + v₂(t, q, p) * db₂dx₃(t, q) + v₃(t, q, p) * db₃dx₃(t, q) +
                 dv₁dq₃(t, q, p) * b₁(t, q) + dv₂dq₃(t, q, p) * b₂(t, q) + dv₃dq₃(t, q, p) * b₃(t, q)

dudp₁(t, q, p) = b₁(t, q)
dudp₂(t, q, p) = b₂(t, q)
dudp₃(t, q, p) = b₃(t, q)


function initial_momentum(tᵢ, qᵢ::AbstractArray{T}) where {T<:Number}
    pᵢ = zeros(T, 3)
    pᵢ[1] = ϑ₁(tᵢ, qᵢ)
    pᵢ[2] = ϑ₂(tᵢ, qᵢ)
    pᵢ[3] = ϑ₃(tᵢ, qᵢ)
    return pᵢ
end

guiding_center_3d_periodicity(::Type{T}, periodic=true) where {T} = periodicity(zeros(T, 3))
guiding_center_3d_periodicity(::AbstractVector{<:AbstractArray{T}}, periodic=true) where {T<:Number} = guiding_center_3d_periodicity(T, periodic)
guiding_center_3d_periodicity(::AbstractArray{T}, periodic=true) where {T<:Number} = guiding_center_3d_periodicity(T, periodic)


hamiltonian(t, q, p, params) = u(t, q, p)^2 / 2 + params.μ * B(t, q)

dHdq₁(t, q, p, params) = u(t, q, p) * dudq₁(t, q, p) + params.μ * dBdx₁(t, q)
dHdq₂(t, q, p, params) = u(t, q, p) * dudq₂(t, q, p) + params.μ * dBdx₂(t, q)
dHdq₃(t, q, p, params) = u(t, q, p) * dudq₃(t, q, p) + params.μ * dBdx₃(t, q)

dHdp₁(t, q, p, params) = u(t, q, p) * dudp₁(t, q, p)
dHdp₂(t, q, p, params) = u(t, q, p) * dudp₂(t, q, p)
dHdp₃(t, q, p, params) = u(t, q, p) * dudp₃(t, q, p)

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

dg₁dq₁(t, q, p) = db₁dx₁(t, q) * v₂(t, q, p) + b₁(t, q) * dv₂dq₁(t, q, p) - db₂dx₁(t, q) * v₁(t, q, p) - b₂(t, q) * dv₁dq₁(t, q, p)
dg₁dq₂(t, q, p) = db₁dx₂(t, q) * v₂(t, q, p) + b₁(t, q) * dv₂dq₂(t, q, p) - db₂dx₂(t, q) * v₁(t, q, p) - b₂(t, q) * dv₁dq₂(t, q, p)
dg₁dq₃(t, q, p) = db₁dx₃(t, q) * v₂(t, q, p) + b₁(t, q) * dv₂dq₃(t, q, p) - db₂dx₃(t, q) * v₁(t, q, p) - b₂(t, q) * dv₁dq₃(t, q, p)

dg₂dq₁(t, q, p) = db₁dx₁(t, q) * v₃(t, q, p) + b₁(t, q) * dv₃dq₁(t, q, p) - db₃dx₁(t, q) * v₁(t, q, p) - b₃(t, q) * dv₁dq₁(t, q, p)
dg₂dq₂(t, q, p) = db₁dx₂(t, q) * v₃(t, q, p) + b₁(t, q) * dv₃dq₂(t, q, p) - db₃dx₂(t, q) * v₁(t, q, p) - b₃(t, q) * dv₁dq₂(t, q, p)
dg₂dq₃(t, q, p) = db₁dx₃(t, q) * v₃(t, q, p) + b₁(t, q) * dv₃dq₃(t, q, p) - db₃dx₃(t, q) * v₁(t, q, p) - b₃(t, q) * dv₁dq₃(t, q, p)

dg₁dp₁(t, q, p) = -b₂(t, q)
dg₁dp₂(t, q, p) = +b₁(t, q)
dg₁dp₃(t, q, p) = zero(eltype(p))

dg₂dp₁(t, q, p) = -b₃(t, q)
dg₂dp₂(t, q, p) = zero(eltype(p))
dg₂dp₃(t, q, p) = +b₁(t, q)

λₒ(t, q, p) = (dg₁dq₁(t, q, p) * dg₂dp₁(t, q, p)
               + dg₁dq₂(t, q, p) * dg₂dp₂(t, q, p)
               + dg₁dq₃(t, q, p) * dg₂dp₃(t, q, p)
               -
               dg₁dp₁(t, q, p) * dg₂dq₁(t, q, p)
               -
               dg₁dp₂(t, q, p) * dg₂dq₂(t, q, p)
               -
               dg₁dp₃(t, q, p) * dg₁dq₃(t, q, p)
)

λ₁(t, q, p, params) = (dg₂dq₁(t, q, p) * dHdp₁(t, q, p, params)
                       + dg₂dq₂(t, q, p) * dHdp₂(t, q, p, params)
                       + dg₂dq₃(t, q, p) * dHdp₃(t, q, p, params)
                       -
                       dg₂dp₁(t, q, p) * dHdq₁(t, q, p, params)
                       -
                       dg₂dp₂(t, q, p) * dHdq₂(t, q, p, params)
                       -
                       dg₂dp₃(t, q, p) * dHdq₃(t, q, p, params)
) / λₒ(t, q, p)

λ₂(t, q, p, params) = -(dg₁dq₁(t, q, p) * dHdp₁(t, q, p, params)
                        + dg₁dq₂(t, q, p) * dHdp₂(t, q, p, params)
                        + dg₁dq₃(t, q, p) * dHdp₃(t, q, p, params)
                        -
                        dg₁dp₁(t, q, p) * dHdq₁(t, q, p, params)
                        -
                        dg₁dp₂(t, q, p) * dHdq₂(t, q, p, params)
                        -
                        dg₁dp₃(t, q, p) * dHdq₃(t, q, p, params)
) / λₒ(t, q, p)


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

function hode(qᵢ=qᵢ, parameters=parameters; tspan=tspan, tstep=Δt, periodic=true)
    t₀ = tspan[begin]
    q₀ = qᵢ[1:3]
    p₀ = initial_momentum(t₀, qᵢ)

    # println("3D Guiding Center model initial constraints g₁ = $(g₁(t₀, q₀, p₀)) and g₂ = $(g₂(t₀, q₀, p₀))")

    HODEProblem(
        guiding_center_3d_v,
        guiding_center_3d_f,
        hamiltonian,
        tspan, tstep, q₀, p₀;
        parameters=parameters,
        periodicity=guiding_center_3d_periodicity(qᵢ, periodic)
    )
end
