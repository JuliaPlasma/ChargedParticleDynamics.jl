#
# The canonicalised form of the 3D guiding centre model: `H̃ = H + Σₖ λᵏ gᵏ`, Eq. (13) of Li, Zhang &
# Liu. Because the multipliers are now differentiated, this needs the second derivatives of both the
# Hamiltonian and the constraints. The former are written out per index below; the latter come from
# `guiding_center_3d_constraints.jl` in one indexed definition each.
#

# The named second derivatives of `v`, kept for the Hamiltonian block below, as aliases of the
# generic accessor. Only the upper triangle exists — `∂²v/∂qⱼ∂qₖ` is symmetric.
for i in 1:3, j in 1:3, k in j:3
    @eval $(Symbol("d²v", INDEX_SUBSCRIPTS[i], "dq", INDEX_SUBSCRIPTS[j], "dq", INDEX_SUBSCRIPTS[k]))(t, q, p) =
        d²vᵢdqⱼdqₖ($(Val(i)), $(Val(j)), $(Val(k)), t, q, p)
end


d²Hdq₁dq₁(t, q, p, params) = (
    dv₁dq₁(t, q, p) * g¹¹(t, q) * dv₁dq₁(t, q, p) +
    dv₂dq₁(t, q, p) * g²²(t, q) * dv₂dq₁(t, q, p) +
    dv₃dq₁(t, q, p) * g³³(t, q) * dv₃dq₁(t, q, p) +
    v₁(t, q, p) * g¹¹(t, q) * d²v₁dq₁dq₁(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * d²v₂dq₁dq₁(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * d²v₃dq₁dq₁(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₁(t, q) * dv₁dq₁(t, q, p) * 2 +
    v₂(t, q, p) * dg²²dx₁(t, q) * dv₂dq₁(t, q, p) * 2 +
    v₃(t, q, p) * dg³³dx₁(t, q) * dv₃dq₁(t, q, p) * 2 +
    v₁(t, q, p) * d²g¹¹dx₁dx₁(t, q) * v₁(t, q, p) / 2 +
    v₂(t, q, p) * d²g²²dx₁dx₁(t, q) * v₂(t, q, p) / 2 +
    v₃(t, q, p) * d²g³³dx₁dx₁(t, q) * v₃(t, q, p) / 2 +
    params.μ * d²Bdx₁dx₁(t, q) -
    dE₁dx₁(t, q)
)

# The metric term of ∂²H/∂qᵢ∂qⱼ is vₖ (∂gᵏᵏ/∂qᵢ ∂vₖ/∂qⱼ + ∂gᵏᵏ/∂qⱼ ∂vₖ/∂qᵢ). For i = j the two
# summands coincide and it collapses to a factor of two; for i ≠ j they do not, and writing the
# factor of two there instead made the Hessian asymmetric (by 65 % on the ITER Solov'ev X-point
# equilibrium) wherever the metric is not constant, i.e. for every curvilinear coordinate system.
d²Hdq₁dq₂(t, q, p, params) = (
    dv₁dq₂(t, q, p) * g¹¹(t, q) * dv₁dq₁(t, q, p) +
    dv₂dq₂(t, q, p) * g²²(t, q) * dv₂dq₁(t, q, p) +
    dv₃dq₂(t, q, p) * g³³(t, q) * dv₃dq₁(t, q, p) +
    v₁(t, q, p) * g¹¹(t, q) * d²v₁dq₁dq₂(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * d²v₂dq₁dq₂(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * d²v₃dq₁dq₂(t, q, p) +
    v₁(t, q, p) * (dg¹¹dx₁(t, q) * dv₁dq₂(t, q, p) + dg¹¹dx₂(t, q) * dv₁dq₁(t, q, p)) +
    v₂(t, q, p) * (dg²²dx₁(t, q) * dv₂dq₂(t, q, p) + dg²²dx₂(t, q) * dv₂dq₁(t, q, p)) +
    v₃(t, q, p) * (dg³³dx₁(t, q) * dv₃dq₂(t, q, p) + dg³³dx₂(t, q) * dv₃dq₁(t, q, p)) +
    v₁(t, q, p) * d²g¹¹dx₁dx₂(t, q) * v₁(t, q, p) / 2 +
    v₂(t, q, p) * d²g²²dx₁dx₂(t, q) * v₂(t, q, p) / 2 +
    v₃(t, q, p) * d²g³³dx₁dx₂(t, q) * v₃(t, q, p) / 2 +
    params.μ * d²Bdx₁dx₂(t, q) -
    dE₁dx₂(t, q)
)

d²Hdq₁dq₃(t, q, p, params) = (
    dv₁dq₃(t, q, p) * g¹¹(t, q) * dv₁dq₁(t, q, p) +
    dv₂dq₃(t, q, p) * g²²(t, q) * dv₂dq₁(t, q, p) +
    dv₃dq₃(t, q, p) * g³³(t, q) * dv₃dq₁(t, q, p) +
    v₁(t, q, p) * g¹¹(t, q) * d²v₁dq₁dq₃(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * d²v₂dq₁dq₃(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * d²v₃dq₁dq₃(t, q, p) +
    v₁(t, q, p) * (dg¹¹dx₁(t, q) * dv₁dq₃(t, q, p) + dg¹¹dx₃(t, q) * dv₁dq₁(t, q, p)) +
    v₂(t, q, p) * (dg²²dx₁(t, q) * dv₂dq₃(t, q, p) + dg²²dx₃(t, q) * dv₂dq₁(t, q, p)) +
    v₃(t, q, p) * (dg³³dx₁(t, q) * dv₃dq₃(t, q, p) + dg³³dx₃(t, q) * dv₃dq₁(t, q, p)) +
    v₁(t, q, p) * d²g¹¹dx₁dx₃(t, q) * v₁(t, q, p) / 2 +
    v₂(t, q, p) * d²g²²dx₁dx₃(t, q) * v₂(t, q, p) / 2 +
    v₃(t, q, p) * d²g³³dx₁dx₃(t, q) * v₃(t, q, p) / 2 +
    params.μ * d²Bdx₁dx₃(t, q) -
    dE₁dx₃(t, q)
)

d²Hdq₁dp₁(t, q, p, params) = (
    dv₁dp₁(t, q, p) * g¹¹(t, q) * dv₁dq₁(t, q, p) +
    dv₂dp₁(t, q, p) * g²²(t, q) * dv₂dq₁(t, q, p) +
    dv₃dp₁(t, q, p) * g³³(t, q) * dv₃dq₁(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₁(t, q) * dv₁dp₁(t, q, p) +
    v₂(t, q, p) * dg²²dx₁(t, q) * dv₂dp₁(t, q, p) +
    v₃(t, q, p) * dg³³dx₁(t, q) * dv₃dp₁(t, q, p)
)

d²Hdq₁dp₂(t, q, p, params) = (
    dv₁dp₂(t, q, p) * g¹¹(t, q) * dv₁dq₁(t, q, p) +
    dv₂dp₂(t, q, p) * g²²(t, q) * dv₂dq₁(t, q, p) +
    dv₃dp₂(t, q, p) * g³³(t, q) * dv₃dq₁(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₁(t, q) * dv₁dp₂(t, q, p) +
    v₂(t, q, p) * dg²²dx₁(t, q) * dv₂dp₂(t, q, p) +
    v₃(t, q, p) * dg³³dx₁(t, q) * dv₃dp₂(t, q, p)
)

d²Hdq₁dp₃(t, q, p, params) = (
    dv₁dp₃(t, q, p) * g¹¹(t, q) * dv₁dq₁(t, q, p) +
    dv₂dp₃(t, q, p) * g²²(t, q) * dv₂dq₁(t, q, p) +
    dv₃dp₃(t, q, p) * g³³(t, q) * dv₃dq₁(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₁(t, q) * dv₁dp₃(t, q, p) +
    v₂(t, q, p) * dg²²dx₁(t, q) * dv₂dp₃(t, q, p) +
    v₃(t, q, p) * dg³³dx₁(t, q) * dv₃dp₃(t, q, p)
)


d²Hdq₂dq₂(t, q, p, params) = (
    dv₁dq₂(t, q, p) * g¹¹(t, q) * dv₁dq₂(t, q, p) +
    dv₂dq₂(t, q, p) * g²²(t, q) * dv₂dq₂(t, q, p) +
    dv₃dq₂(t, q, p) * g³³(t, q) * dv₃dq₂(t, q, p) +
    v₁(t, q, p) * g¹¹(t, q) * d²v₁dq₂dq₂(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * d²v₂dq₂dq₂(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * d²v₃dq₂dq₂(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₂(t, q) * dv₁dq₂(t, q, p) * 2 +
    v₂(t, q, p) * dg²²dx₂(t, q) * dv₂dq₂(t, q, p) * 2 +
    v₃(t, q, p) * dg³³dx₂(t, q) * dv₃dq₂(t, q, p) * 2 +
    v₁(t, q, p) * d²g¹¹dx₂dx₂(t, q) * v₁(t, q, p) / 2 +
    v₂(t, q, p) * d²g²²dx₂dx₂(t, q) * v₂(t, q, p) / 2 +
    v₃(t, q, p) * d²g³³dx₂dx₂(t, q) * v₃(t, q, p) / 2 +
    params.μ * d²Bdx₂dx₂(t, q) -
    dE₂dx₂(t, q)
)

d²Hdq₂dq₃(t, q, p, params) = (
    dv₁dq₃(t, q, p) * g¹¹(t, q) * dv₁dq₂(t, q, p) +
    dv₂dq₃(t, q, p) * g²²(t, q) * dv₂dq₂(t, q, p) +
    dv₃dq₃(t, q, p) * g³³(t, q) * dv₃dq₂(t, q, p) +
    v₁(t, q, p) * g¹¹(t, q) * d²v₁dq₂dq₃(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * d²v₂dq₂dq₃(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * d²v₃dq₂dq₃(t, q, p) +
    v₁(t, q, p) * (dg¹¹dx₂(t, q) * dv₁dq₃(t, q, p) + dg¹¹dx₃(t, q) * dv₁dq₂(t, q, p)) +
    v₂(t, q, p) * (dg²²dx₂(t, q) * dv₂dq₃(t, q, p) + dg²²dx₃(t, q) * dv₂dq₂(t, q, p)) +
    v₃(t, q, p) * (dg³³dx₂(t, q) * dv₃dq₃(t, q, p) + dg³³dx₃(t, q) * dv₃dq₂(t, q, p)) +
    v₁(t, q, p) * d²g¹¹dx₂dx₃(t, q) * v₁(t, q, p) / 2 +
    v₂(t, q, p) * d²g²²dx₂dx₃(t, q) * v₂(t, q, p) / 2 +
    v₃(t, q, p) * d²g³³dx₂dx₃(t, q) * v₃(t, q, p) / 2 +
    params.μ * d²Bdx₂dx₃(t, q) -
    dE₂dx₃(t, q)
)

d²Hdq₂dp₁(t, q, p, params) = (
    dv₁dp₁(t, q, p) * g¹¹(t, q) * dv₁dq₂(t, q, p) +
    dv₂dp₁(t, q, p) * g²²(t, q) * dv₂dq₂(t, q, p) +
    dv₃dp₁(t, q, p) * g³³(t, q) * dv₃dq₂(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₂(t, q) * dv₁dp₁(t, q, p) +
    v₂(t, q, p) * dg²²dx₂(t, q) * dv₂dp₁(t, q, p) +
    v₃(t, q, p) * dg³³dx₂(t, q) * dv₃dp₁(t, q, p)
)

d²Hdq₂dp₂(t, q, p, params) = (
    dv₁dp₂(t, q, p) * g¹¹(t, q) * dv₁dq₂(t, q, p) +
    dv₂dp₂(t, q, p) * g²²(t, q) * dv₂dq₂(t, q, p) +
    dv₃dp₂(t, q, p) * g³³(t, q) * dv₃dq₂(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₂(t, q) * dv₁dp₂(t, q, p) +
    v₂(t, q, p) * dg²²dx₂(t, q) * dv₂dp₂(t, q, p) +
    v₃(t, q, p) * dg³³dx₂(t, q) * dv₃dp₂(t, q, p)
)

d²Hdq₂dp₃(t, q, p, params) = (
    dv₁dp₃(t, q, p) * g¹¹(t, q) * dv₁dq₂(t, q, p) +
    dv₂dp₃(t, q, p) * g²²(t, q) * dv₂dq₂(t, q, p) +
    dv₃dp₃(t, q, p) * g³³(t, q) * dv₃dq₂(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₂(t, q) * dv₁dp₃(t, q, p) +
    v₂(t, q, p) * dg²²dx₂(t, q) * dv₂dp₃(t, q, p) +
    v₃(t, q, p) * dg³³dx₂(t, q) * dv₃dp₃(t, q, p)
)


d²Hdq₃dq₃(t, q, p, params) = (
    dv₁dq₃(t, q, p) * g¹¹(t, q) * dv₁dq₃(t, q, p) +
    dv₂dq₃(t, q, p) * g²²(t, q) * dv₂dq₃(t, q, p) +
    dv₃dq₃(t, q, p) * g³³(t, q) * dv₃dq₃(t, q, p) +
    v₁(t, q, p) * g¹¹(t, q) * d²v₁dq₃dq₃(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * d²v₂dq₃dq₃(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * d²v₃dq₃dq₃(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₃(t, q) * dv₁dq₃(t, q, p) * 2 +
    v₂(t, q, p) * dg²²dx₃(t, q) * dv₂dq₃(t, q, p) * 2 +
    v₃(t, q, p) * dg³³dx₃(t, q) * dv₃dq₃(t, q, p) * 2 +
    v₁(t, q, p) * d²g¹¹dx₃dx₃(t, q) * v₁(t, q, p) / 2 +
    v₂(t, q, p) * d²g²²dx₃dx₃(t, q) * v₂(t, q, p) / 2 +
    v₃(t, q, p) * d²g³³dx₃dx₃(t, q) * v₃(t, q, p) / 2 +
    params.μ * d²Bdx₃dx₃(t, q) -
    dE₃dx₃(t, q)
)

d²Hdq₃dp₁(t, q, p, params) = (
    dv₁dp₁(t, q, p) * g¹¹(t, q) * dv₁dq₃(t, q, p) +
    dv₂dp₁(t, q, p) * g²²(t, q) * dv₂dq₃(t, q, p) +
    dv₃dp₁(t, q, p) * g³³(t, q) * dv₃dq₃(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₃(t, q) * dv₁dp₁(t, q, p) +
    v₂(t, q, p) * dg²²dx₃(t, q) * dv₂dp₁(t, q, p) +
    v₃(t, q, p) * dg³³dx₃(t, q) * dv₃dp₁(t, q, p)
)

d²Hdq₃dp₂(t, q, p, params) = (
    dv₁dp₂(t, q, p) * g¹¹(t, q) * dv₁dq₃(t, q, p) +
    dv₂dp₂(t, q, p) * g²²(t, q) * dv₂dq₃(t, q, p) +
    dv₃dp₂(t, q, p) * g³³(t, q) * dv₃dq₃(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₃(t, q) * dv₁dp₂(t, q, p) +
    v₂(t, q, p) * dg²²dx₃(t, q) * dv₂dp₂(t, q, p) +
    v₃(t, q, p) * dg³³dx₃(t, q) * dv₃dp₂(t, q, p)
)

d²Hdq₃dp₃(t, q, p, params) = (
    dv₁dp₃(t, q, p) * g¹¹(t, q) * dv₁dq₃(t, q, p) +
    dv₂dp₃(t, q, p) * g²²(t, q) * dv₂dq₃(t, q, p) +
    dv₃dp₃(t, q, p) * g³³(t, q) * dv₃dq₃(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₃(t, q) * dv₁dp₃(t, q, p) +
    v₂(t, q, p) * dg²²dx₃(t, q) * dv₂dp₃(t, q, p) +
    v₃(t, q, p) * dg³³dx₃(t, q) * dv₃dp₃(t, q, p)
)


d²Hdp₁dp₁(t, q, p, params) = (
    dv₁dp₁(t, q, p) * g¹¹(t, q) * dv₁dp₁(t, q, p) +
    dv₂dp₁(t, q, p) * g²²(t, q) * dv₂dp₁(t, q, p) +
    dv₃dp₁(t, q, p) * g³³(t, q) * dv₃dp₁(t, q, p)
)


d²Hdp₁dp₂(t, q, p, params) = (
    dv₁dp₂(t, q, p) * g¹¹(t, q) * dv₁dp₁(t, q, p) +
    dv₂dp₂(t, q, p) * g²²(t, q) * dv₂dp₁(t, q, p) +
    dv₃dp₂(t, q, p) * g³³(t, q) * dv₃dp₁(t, q, p)
)

d²Hdp₂dp₂(t, q, p, params) = (
    dv₁dp₂(t, q, p) * g¹¹(t, q) * dv₁dp₂(t, q, p) +
    dv₂dp₂(t, q, p) * g²²(t, q) * dv₂dp₂(t, q, p) +
    dv₃dp₂(t, q, p) * g³³(t, q) * dv₃dp₂(t, q, p)
)


d²Hdp₁dp₃(t, q, p, params) = (
    dv₁dp₃(t, q, p) * g¹¹(t, q) * dv₁dp₁(t, q, p) +
    dv₂dp₃(t, q, p) * g²²(t, q) * dv₂dp₁(t, q, p) +
    dv₃dp₃(t, q, p) * g³³(t, q) * dv₃dp₁(t, q, p)
)

d²Hdp₂dp₃(t, q, p, params) = (
    dv₁dp₃(t, q, p) * g¹¹(t, q) * dv₁dp₂(t, q, p) +
    dv₂dp₃(t, q, p) * g²²(t, q) * dv₂dp₂(t, q, p) +
    dv₃dp₃(t, q, p) * g³³(t, q) * dv₃dp₂(t, q, p)
)

d²Hdp₃dp₃(t, q, p, params) = (
    dv₁dp₃(t, q, p) * g¹¹(t, q) * dv₁dp₃(t, q, p) +
    dv₂dp₃(t, q, p) * g²²(t, q) * dv₂dp₃(t, q, p) +
    dv₃dp₃(t, q, p) * g³³(t, q) * dv₃dp₃(t, q, p)
)


# Indexed views of the Hessian blocks above. Only the upper triangle of the `qq` and `pp` blocks is
# written out, so those two symmetrise here; the `qp` block is not symmetric and has all nine.
for i in 1:3, j in 1:3
    @eval @inline d²Hdqᵢdqⱼ(::Val{$i}, ::Val{$j}, t, q, p, params) =
        $(Symbol("d²Hdq", INDEX_SUBSCRIPTS[min(i, j)], "dq", INDEX_SUBSCRIPTS[max(i, j)]))(t, q, p, params)
    @eval @inline d²Hdpᵢdpⱼ(::Val{$i}, ::Val{$j}, t, q, p, params) =
        $(Symbol("d²Hdp", INDEX_SUBSCRIPTS[min(i, j)], "dp", INDEX_SUBSCRIPTS[max(i, j)]))(t, q, p, params)
    @eval @inline d²Hdqᵢdpⱼ(::Val{$i}, ::Val{$j}, t, q, p, params) =
        $(Symbol("d²Hdq", INDEX_SUBSCRIPTS[i], "dp", INDEX_SUBSCRIPTS[j]))(t, q, p, params)
end


# Derivatives of the two brackets of `guiding_center_3d_equations.jl`. The `∂²gᵏ/∂pₗ∂pₘ` terms that
# the product rule would contribute vanish identically and are left out.
@inline dbracket_gg_dqⱼ(k₁::Val, k₂::Val, j::Val, t, q, p) = contract(l ->
    d²gᵏdqₗdqₘ(k₁, j, l, t, q, p) * dgᵏdpₗ(k₂, l, t, q, p) +
    dgᵏdqₗ(k₁, l, t, q, p) * d²gᵏdqₗdpₘ(k₂, j, l, t, q, p) -
    d²gᵏdqₗdpₘ(k₁, j, l, t, q, p) * dgᵏdqₗ(k₂, l, t, q, p) -
    dgᵏdpₗ(k₁, l, t, q, p) * d²gᵏdqₗdqₘ(k₂, j, l, t, q, p))

@inline dbracket_gg_dpⱼ(k₁::Val, k₂::Val, j::Val, t, q, p) = contract(l ->
    d²gᵏdqₗdpₘ(k₁, l, j, t, q, p) * dgᵏdpₗ(k₂, l, t, q, p) -
    dgᵏdpₗ(k₁, l, t, q, p) * d²gᵏdqₗdpₘ(k₂, l, j, t, q, p))

@inline dbracket_gH_dqⱼ(k::Val, j::Val, t, q, p, params) = contract(l ->
    d²gᵏdqₗdqₘ(k, j, l, t, q, p) * dHdpᵢ(l, t, q, p, params) +
    dgᵏdqₗ(k, l, t, q, p) * d²Hdqᵢdpⱼ(j, l, t, q, p, params) -
    d²gᵏdqₗdpₘ(k, j, l, t, q, p) * dHdqᵢ(l, t, q, p, params) -
    dgᵏdpₗ(k, l, t, q, p) * d²Hdqᵢdqⱼ(j, l, t, q, p, params))

@inline dbracket_gH_dpⱼ(k::Val, j::Val, t, q, p, params) = contract(l ->
    d²gᵏdqₗdpₘ(k, l, j, t, q, p) * dHdpᵢ(l, t, q, p, params) +
    dgᵏdqₗ(k, l, t, q, p) * d²Hdpᵢdpⱼ(j, l, t, q, p, params) -
    dgᵏdpₗ(k, l, t, q, p) * d²Hdqᵢdpⱼ(l, j, t, q, p, params))


dλₒdqⱼ(j::Val, t, q, p, c) = dbracket_gg_dqⱼ(c[1], c[2], j, t, q, p)
dλₒdpⱼ(j::Val, t, q, p, c) = dbracket_gg_dpⱼ(c[1], c[2], j, t, q, p)

dλ₁dqⱼ(j::Val, t, q, p, params, c) =
    -dλₒdqⱼ(j, t, q, p, c) * λ₁(t, q, p, params, c) / λₒ(t, q, p, c) +
    dbracket_gH_dqⱼ(c[2], j, t, q, p, params) / λₒ(t, q, p, c)

dλ₁dpⱼ(j::Val, t, q, p, params, c) =
    -dλₒdpⱼ(j, t, q, p, c) * λ₁(t, q, p, params, c) / λₒ(t, q, p, c) +
    dbracket_gH_dpⱼ(c[2], j, t, q, p, params) / λₒ(t, q, p, c)

# Both multipliers are a bracket over `λₒ`, so both carry the same `-λ ∂λₒ/λₒ` term regardless of the
# sign of the numerator: for `λ = N/λₒ`, `∂λ = ∂N/λₒ - λ ∂λₒ/λₒ`. `dλ₂` carried that term with a `+`,
# which is the sign the quotient rule gives before `-{g₁,H}` is folded back into `λ₂`. The error is
# proportional to `∂λₒ`, and the whole term is multiplied by a constraint in the right-hand side, so
# it vanished on the constraint manifold and showed up only as constraint drift — see the
# finite-difference test in `test/structure_tests.jl`.
dλ₂dqⱼ(j::Val, t, q, p, params, c) =
    -dλₒdqⱼ(j, t, q, p, c) * λ₂(t, q, p, params, c) / λₒ(t, q, p, c) -
    dbracket_gH_dqⱼ(c[1], j, t, q, p, params) / λₒ(t, q, p, c)

dλ₂dpⱼ(j::Val, t, q, p, params, c) =
    -dλₒdpⱼ(j, t, q, p, c) * λ₂(t, q, p, params, c) / λₒ(t, q, p, c) -
    dbracket_gH_dpⱼ(c[1], j, t, q, p, params) / λₒ(t, q, p, c)


function guiding_center_3d_canonical_v(v, t, q, p, params, c = default_constraint_pair())
    guiding_center_3d_v(v, t, q, p, params, c)

    g1 = gᵏ(c[1], t, q, p)
    g2 = gᵏ(c[2], t, q, p)

    addcomponents!(v, i -> dλ₁dpⱼ(i, t, q, p, params, c) * g1 +
                           dλ₂dpⱼ(i, t, q, p, params, c) * g2)
    nothing
end

function guiding_center_3d_canonical_f(f, t, q, p, params, c = default_constraint_pair())
    guiding_center_3d_f(f, t, q, p, params, c)

    g1 = gᵏ(c[1], t, q, p)
    g2 = gᵏ(c[2], t, q, p)

    addcomponents!(f, i -> -dλ₁dqⱼ(i, t, q, p, params, c) * g1 -
                            dλ₂dqⱼ(i, t, q, p, params, c) * g2)
    nothing
end


"""
    hodeproblem_canonical(q₀, p₀; kwargs...)
    hodeproblem_canonical(x₀ = qᵢ; kwargs...)
    hodeproblem_canonical(ics::NamedTuple; kwargs...)

The canonicalised guiding centre system, `H̃ = H + λ₁ g₁ + λ₂ g₂`, which agrees with
[`hodeproblem`](@ref) on the constraint manifold but is genuinely canonical. Same three argument
forms and the same `constraints` keyword as `hodeproblem`.
"""
function hodeproblem_canonical(q₀::AbstractVector, p₀::AbstractVector; timespan = DEFAULT_TIMESPAN,
                               timestep = DEFAULT_TIMESTEP, parameters = default_parameters(), periodic = true,
                               constraints = default_constraints())
    c = constraint_pair(constraints)

    _v(v, t, q, p, params) = guiding_center_3d_canonical_v(v, t, q, p, params, c)
    _f(f, t, q, p, params) = guiding_center_3d_canonical_f(f, t, q, p, params, c)
    _h(t, q, p, params) = hamiltonian_canonical(t, q, p, params, c)

    HODEProblem(
        _v,
        _f,
        _h,
        timespan, timestep, q₀, p₀;
        parameters=parameters,
        periodicity=guiding_center_3d_periodicity(q₀, periodic))
end

function hodeproblem_canonical(x₀::AbstractVector = qᵢ; timespan = DEFAULT_TIMESPAN, kwargs...)
    ics = initial_conditions(timespan[begin], x₀)
    hodeproblem_canonical(ics.q, ics.p; timespan = timespan, kwargs...)
end

hodeproblem_canonical(ics::NamedTuple; kwargs...) =
    hodeproblem_canonical(ics.q, ics.p; parameters = ics.params, kwargs...)
