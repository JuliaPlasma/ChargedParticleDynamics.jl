# Common functions for (r,θ,ϕ) tokamak coordinates.

using GeometricIntegrators.Equations
using GeometricIntegrators.Utils


function R(t, q)
    R₀ + q[1] * cos(q[2])
end

function Z(t, q)
    q[1] * sin(q[2])
end

function ϕ(t, q)
    q[3]
end

function u(t, q)
    q[4]
end


function dRdx₁(t,q)
    cos(q[2])
end

function dRdx₂(t,q)
    - q[1] * sin(q[2])
end

function dRdx₃(t,q)
    zero(eltype(q))
end


function α1(t, q)
    A₁(t,q) + u(t,q) * b₁(t,q)
end

function α2(t, q)
    A₂(t,q) + u(t,q) * b₂(t,q)
end

function α3(t, q)
    R(t,q) * ( A₃(t,q) + u(t,q) * b₃(t,q) )
end

function α4(t, q)
    zero(eltype(q))
end


function dα1d1(t, q)
    dA₁dx₁(t,q) + u(t,q) * db₁dx₁(t,q)
end

function dα1d2(t, q)
    dA₁dx₂(t,q) + u(t,q) * db₁dx₂(t,q)
end

function dα1d3(t, q)
    dA₁dx₃(t,q) + u(t,q) * db₁dx₃(t,q)
end

function dα1d4(t, q)
    b₁(t,q)
end

function dα2d1(t, q)
    dA₂dx₁(t,q) + u(t,q) * db₂dx₁(t,q)
end

function dα2d2(t, q)
    dA₂dx₂(t,q) + u(t,q) * db₂dx₂(t,q)
end

function dα2d3(t, q)
    dA₂dx₃(t,q) + u(t,q) * db₂dx₃(t,q)
end

function dα2d4(t, q)
    b₂(t,q)
end

function dα3d1(t, q)
    R(t,q) * ( dA₃dx₁(t,q) + u(t,q) * db₃dx₁(t,q) ) + dRdx₁(t,q) * ( A₃(t,q) + u(t,q) * b₃(t,q) )
end

function dα3d2(t, q)
    R(t,q) * ( dA₃dx₂(t,q) + u(t,q) * db₃dx₂(t,q) ) + dRdx₂(t,q) * ( A₃(t,q) + u(t,q) * b₃(t,q) )
end

function dα3d3(t, q)
    R(t,q) * ( dA₃dx₃(t,q) + u(t,q) * db₃dx₃(t,q) ) + dRdx₃(t,q) * ( A₃(t,q) + u(t,q) * b₃(t,q) )
end

function dα3d4(t, q)
    R(t,q) * b₃(t,q)
end


# function compute_parallel_velocity(t,x,pᵤ)
#     - (pᵤ + 0.5 * ( x[2]^2 + (R₀ - x[1])^2 ) * B₀ / q) * B(t,x) / (B₀ * R₀)
# end


include("guiding_center_4d_common.jl")
