"""
Uniform magnetic field of the form
``B(x,y,z) = B_0 e_z``.
"""

const A₀ = 1

A₁(t,q) = - 0.5 * A₀ * q[2]
A₂(t,q) = + 0.5 * A₀ * q[1]
A₃(t,q) = zero(eltype(q))

dA₁dx₁(t,q) = zero(eltype(q))
dA₁dx₂(t,q) = - 0.5 * A₀
dA₁dx₃(t,q) = zero(eltype(q))

dA₂dx₁(t,q) = + 0.5 * A₀
dA₂dx₂(t,q) = zero(eltype(q))
dA₂dx₃(t,q) = zero(eltype(q))

dA₃dx₁(t,q) = zero(eltype(q))
dA₃dx₂(t,q) = zero(eltype(q))
dA₃dx₃(t,q) = zero(eltype(q))

B(t,q) = eltype(q)(A₀)

dBdx₁(t,q) = zero(eltype(q))
dBdx₂(t,q) = zero(eltype(q))
dBdx₃(t,q) = zero(eltype(q))

B₁(t,q) = zero(eltype(q))
B₂(t,q) = zero(eltype(q))
B₃(t,q) = B(t,q)

b₁(t,q) = zero(eltype(q))
b₂(t,q) = zero(eltype(q))
b₃(t,q) = one(eltype(q))

b¹(t,q) = zero(eltype(q))
b²(t,q) = zero(eltype(q))
b³(t,q) = one(eltype(q))

db₁dx₁(t,q) = zero(eltype(q))
db₁dx₂(t,q) = zero(eltype(q))
db₁dx₃(t,q) = zero(eltype(q))

db₂dx₁(t,q) = zero(eltype(q))
db₂dx₂(t,q) = zero(eltype(q))
db₂dx₃(t,q) = zero(eltype(q))

db₃dx₁(t,q) = zero(eltype(q))
db₃dx₂(t,q) = zero(eltype(q))
db₃dx₃(t,q) = zero(eltype(q))

R(t,x) = one(eltype(x))

dRdx₁(t,x) = zero(eltype(x))
dRdx₂(t,x) = zero(eltype(x))
dRdx₃(t,x) = zero(eltype(x))

g₁₁(t,x) = one(eltype(x))
g₂₂(t,x) = one(eltype(x))
g₃₃(t,x) = one(eltype(x))

g¹¹(t,x) = one(eltype(x))
g²²(t,x) = one(eltype(x))
g³³(t,x) = one(eltype(x))
