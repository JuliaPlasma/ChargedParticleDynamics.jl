"""
Axisymmetric magnetic field of the form
``B(x,y,z) = (1 + x^2 + y^2) e_z``.
"""

A₁(t,x) = - x[2] * (2 + x[1]^2 + x[2]^2) / 4
A₂(t,x) = + x[1] * (2 + x[1]^2 + x[2]^2) / 4
A₃(t,x) = zero(eltype(x))

dA₁dx₁(t,x) = - x[1] * x[2] / 2
dA₁dx₂(t,x) = - (2 + x[1]^2 + 3 * x[2]^2) / 4
dA₁dx₃(t,x) = zero(eltype(x))

dA₂dx₁(t,x) = + (2 + 3 * x[1]^2 + x[2]^2) / 4
dA₂dx₂(t,x) = + x[1] * x[2] / 2
dA₂dx₃(t,x) = zero(eltype(x))

dA₃dx₁(t,x) = zero(eltype(x))
dA₃dx₂(t,x) = zero(eltype(x))
dA₃dx₃(t,x) =zero(eltype(x))

B₁(t,x) = zero(eltype(x))
B₂(t,x) = zero(eltype(x))
B₃(t,x) = x[1]^2 + x[2]^2 + 1

B(t,x) = x[1]^2 + x[2]^2 + 1

b₁(t,x) = zero(eltype(x))
b₂(t,x) = zero(eltype(x))
b₃(t,x) = one(eltype(x))

b(t,x) = [b₁(t,x), b₂(t,x), b₃(t,x)]

R(t,x) = one(eltype(x))
Z(t,x) = zero(eltype(x))
phi(t,x) = zero(eltype(x))

dRdx₁(t,x) = zero(eltype(x))
dRdx₂(t,x) = zero(eltype(x))
dRdx₃(t,x) = zero(eltype(x))

db₁dx₁(t,x) = zero(eltype(x))
db₁dx₂(t,x) = zero(eltype(x))
db₁dx₃(t,x) = zero(eltype(x))
db₂dx₁(t,x) = zero(eltype(x))
db₂dx₂(t,x) = zero(eltype(x))
db₂dx₃(t,x) = zero(eltype(x))
db₃dx₁(t,x) = zero(eltype(x))
db₃dx₂(t,x) = zero(eltype(x))
db₃dx₃(t,x) = zero(eltype(x))

dB₁dx₁(t,x) = zero(eltype(x))
dB₁dx₂(t,x) = zero(eltype(x))
dB₁dx₃(t,x) = zero(eltype(x))
dB₂dx₁(t,x) = zero(eltype(x))
dB₂dx₂(t,x) = zero(eltype(x))
dB₂dx₃(t,x) = zero(eltype(x))
dB₃dx₁(t,x) = 2*x[1]
dB₃dx₂(t,x) = 2*x[2]
dB₃dx₃(t,x) = zero(eltype(x))

dBdx₁(t,x) = 2*x[1]
dBdx₂(t,x) = 2*x[2]
dBdx₃(t,x) = zero(eltype(x))
