
using GeometricIntegrators.Equations
using GeometricIntegrators.Utils


function u(t, q)
    q[4]
end


function α1(t, q)
    A₁(t,q) + u(t,q) * b₁(t,q)
end

function α2(t, q)
    A₂(t,q) + u(t,q) * b₂(t,q)
end

function α3(t, q)
    A₃(t,q) + u(t,q) * b₃(t,q)
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
    dA₃dx₁(t,q) + u(t,q) * db₃dx₁(t,q)
end

function dα3d2(t, q)
    dA₃dx₂(t,q) + u(t,q) * db₃dx₂(t,q)
end

function dα3d3(t, q)
    dA₃dx₃(t,q) + u(t,q) * db₃dx₃(t,q)
end

function dα3d4(t, q)
    b₃(t,q)
end


include("guiding_center_4d_common.jl")
