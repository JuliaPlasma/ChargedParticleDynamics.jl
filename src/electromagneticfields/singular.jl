"""
Singular magnetic field of the form
``B(x,y,z) = (x^2 + y^2)^{-3/2} e_z``.
"""

function denominator(x,y)
   sqrt(x^2 + y^2)
   # sqrt(sqrt(eps()) + x^2 + y^2)
   # sqrt(1 + x^2 + y^2)
end


A₁(t,q) = + q[2] / denominator(q[1], q[2])^3
A₂(t,q) = - q[1] / denominator(q[1], q[2])^3
A₃(t,q) = zero(eltype(q))

dA₁dx₁(t,q) = - 3 * q[1] * q[2] / denominator(q[1], q[2])^5
dA₁dx₂(t,q) = - 3 * q[2] * q[2] / denominator(q[1], q[2])^5 + 1 / denominator(q[1], q[2])^3
dA₁dx₃(t,q) = zero(eltype(q))

dA₂dx₁(t,q) = + 3 * q[1] * q[1] / denominator(q[1], q[2])^5 - 1 / denominator(q[1], q[2])^3
dA₂dx₂(t,q) = + 3 * q[1] * q[2] / denominator(q[1], q[2])^5
dA₂dx₃(t,q) = zero(eltype(q))

dA₃dx₁(t,q) = zero(eltype(q))
dA₃dx₂(t,q) = zero(eltype(q))
dA₃dx₃(t,q) = zero(eltype(q))

B₀(x,y) = 1 / denominator(x, y)^3

B₁(t,q) = zero(eltype(q))
B₂(t,q) = zero(eltype(q))
B₃(t,q) = B₀(q[1], q[2])

B(t,q) = B₀(q[1], q[2])

b₁(t,q) = zero(eltype(q))
b₂(t,q) = zero(eltype(q))
b₃(t,q) = one(eltype(q))
