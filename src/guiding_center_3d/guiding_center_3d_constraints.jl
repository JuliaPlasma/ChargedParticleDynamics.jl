#
# The constraints of the 3D guiding centre model, written once against a running index.
#
# The degenerate Legendre transform of the guiding centre Lagrangian gives `pᵢ = bᵢ (b·Ẋ) + Aᵢ`, so
# the three momenta are tied together by `v × b = 0` with `v = p - A`. Two of the three components
# are independent, and Li, Zhang & Liu label them
#
#     g¹ = b₁ v₃ - b₃ v₁ ,    g² = b₂ v₃ - b₃ v₂ ,    g³ = b₁ v₂ - b₂ v₁ .
#
# Any two of them may serve as the constraint pair. The Poisson bracket `{gⁱ, gʲ}` that divides both
# Lagrange multipliers carries a factor of the `b`-component belonging to the *omitted* constraint —
# `b₃` for the pair `(g¹, g²)`, `b₁` for `(g³, g¹)`, `b₂` for `(g², g³)` — so which pair is usable is
# a property of the equilibrium and of where the orbit sits in it. See `constraint_pair` below.
#
# The three are one indexed family, `gᵏ = b_i v_j - b_j v_i` with `(i, j)` from
# `constraint_indices`, and so are their derivatives. Writing them out per index is what made
# `guiding_center_3d_canonical.jl` seven hundred lines long and is why only one pair was ever
# implemented.
#
# `ElectromagneticFields.@code()` injects the field functions under names carrying literal
# subscripts — `b₁`, `db₂dx₃`, `d²A₁dx₂dx₃` — so the generic expressions need index accessors for
# them. Those take `Val`s rather than `Int`s: the index is then part of the signature, every call
# inlines to the injected function it stands for, and no dynamic dispatch survives into the
# right-hand side. Their names must avoid the injected ones, of which `A`, `a`, `b`, `c`, `B`, `E`,
# `g`, `φ`, `DF`, `J` and the `aₚ`/`a⃗` triads are all taken; hence the `ᵢ`/`ⱼ`/`ₖ` suffixes.
#

const INDEX_SUBSCRIPTS = ('₁', '₂', '₃')

for i in 1:3
    @eval @inline Aᵢ(::Val{$i}, t, q) = $(Symbol("A", INDEX_SUBSCRIPTS[i]))(t, q)
    @eval @inline bᵢ(::Val{$i}, t, q) = $(Symbol("b", INDEX_SUBSCRIPTS[i]))(t, q)

    for j in 1:3
        @eval @inline dAᵢdxⱼ(::Val{$i}, ::Val{$j}, t, q) = $(Symbol("dA", INDEX_SUBSCRIPTS[i], "dx", INDEX_SUBSCRIPTS[j]))(t, q)
        @eval @inline dbᵢdxⱼ(::Val{$i}, ::Val{$j}, t, q) = $(Symbol("db", INDEX_SUBSCRIPTS[i], "dx", INDEX_SUBSCRIPTS[j]))(t, q)

        for k in 1:3
            @eval @inline d²Aᵢdxⱼdxₖ(::Val{$i}, ::Val{$j}, ::Val{$k}, t, q) = $(Symbol("d²A", INDEX_SUBSCRIPTS[i], "dx", INDEX_SUBSCRIPTS[j], "dx", INDEX_SUBSCRIPTS[k]))(t, q)
            @eval @inline d²bᵢdxⱼdxₖ(::Val{$i}, ::Val{$j}, ::Val{$k}, t, q) = $(Symbol("d²b", INDEX_SUBSCRIPTS[i], "dx", INDEX_SUBSCRIPTS[j], "dx", INDEX_SUBSCRIPTS[k]))(t, q)
        end
    end
end


# Contraction over the three coordinate directions. Written as a sum of three calls rather than a
# loop or a generator so that the `Val` index stays a compile-time constant in each summand.
@inline contract(f) = f(Val(1)) + f(Val(2)) + f(Val(3))


vᵢ(::Val{i}, t, q, p) where {i} = p[i] - Aᵢ(Val(i), t, q)

dvᵢdqⱼ(i::Val, j::Val, t, q, p) = -dAᵢdxⱼ(i, j, t, q)
dvᵢdpⱼ(::Val{i}, ::Val{j}, t, q, p) where {i,j} = i == j ? one(eltype(p)) : zero(eltype(p))

d²vᵢdqⱼdqₖ(i::Val, j::Val, k::Val, t, q, p) = -d²Aᵢdxⱼdxₖ(i, j, k, t, q)


"""
    constraint_indices(::Val{k})

The pair `(i, j)` of field indices the `k`-th constraint `gᵏ = bᵢ vⱼ - bⱼ vᵢ` is built from, as
`Val`s so that the caller keeps compile-time indices. The labelling is that of Li, Zhang & Liu:
`g¹ = b₁v₃ - b₃v₁`, `g² = b₂v₃ - b₃v₂`, `g³ = b₁v₂ - b₂v₁`.
"""
constraint_indices(::Val{1}) = (Val(1), Val(3))
constraint_indices(::Val{2}) = (Val(2), Val(3))
constraint_indices(::Val{3}) = (Val(1), Val(2))


"""
    constraint_pair(constraints::Symbol)

The two constraints selected by the `constraints` keyword of the problem constructors, in the order
they enter the multipliers:

| symbol | pair | multiplier denominator |
|---|---|---|
| `:g31` | `(g³, g¹)` | `b₁ [B + (p-A)·(∇×b)]` |
| `:g12` | `(g¹, g²)` | `b₃ [B + (p-A)·(∇×b)]` |
| `:g23` | `(g², g³)` | `b₂ [B + (p-A)·(∇×b)]` |

`:g31` is the pair the package implemented before the others existed. Which of the three is usable
is a property of the equilibrium, so each one declares its own `default_constraints` beside its
`default_parameters`.
"""
constraint_pair(::Val{:g31}) = (Val(3), Val(1))
constraint_pair(::Val{:g12}) = (Val(1), Val(2))
constraint_pair(::Val{:g23}) = (Val(2), Val(3))

constraint_pair(constraints::Symbol) = constraint_pair(Val(constraints))

default_constraint_pair() = constraint_pair(default_constraints())


"""
    gᵏ(::Val{k}, t, q, p)

The `k`-th constraint, `gᵏ = bᵢ vⱼ - bⱼ vᵢ` with `(i, j) = constraint_indices(Val(k))`. All three
vanish on the constraint manifold, whichever pair a problem was built with.
"""
@inline function gᵏ(k::Val, t, q, p)
    i, j = constraint_indices(k)
    bᵢ(i, t, q) * vᵢ(j, t, q, p) - bᵢ(j, t, q) * vᵢ(i, t, q, p)
end

@inline function dgᵏdqₗ(k::Val, l::Val, t, q, p)
    i, j = constraint_indices(k)
    (dbᵢdxⱼ(i, l, t, q) * vᵢ(j, t, q, p) + bᵢ(i, t, q) * dvᵢdqⱼ(j, l, t, q, p)) -
    (dbᵢdxⱼ(j, l, t, q) * vᵢ(i, t, q, p) + bᵢ(j, t, q) * dvᵢdqⱼ(i, l, t, q, p))
end

@inline function dgᵏdpₗ(k::Val, l::Val, t, q, p)
    i, j = constraint_indices(k)
    bᵢ(i, t, q) * dvᵢdpⱼ(j, l, t, q, p) - bᵢ(j, t, q) * dvᵢdpⱼ(i, l, t, q, p)
end

@inline function d²gᵏdqₗdqₘ(k::Val, l::Val, m::Val, t, q, p)
    i, j = constraint_indices(k)
    (d²bᵢdxⱼdxₖ(i, l, m, t, q) * vᵢ(j, t, q, p) +
     dbᵢdxⱼ(i, l, t, q) * dvᵢdqⱼ(j, m, t, q, p) +
     dbᵢdxⱼ(i, m, t, q) * dvᵢdqⱼ(j, l, t, q, p) +
     bᵢ(i, t, q) * d²vᵢdqⱼdqₖ(j, l, m, t, q, p)) -
    (d²bᵢdxⱼdxₖ(j, l, m, t, q) * vᵢ(i, t, q, p) +
     dbᵢdxⱼ(j, l, t, q) * dvᵢdqⱼ(i, m, t, q, p) +
     dbᵢdxⱼ(j, m, t, q) * dvᵢdqⱼ(i, l, t, q, p) +
     bᵢ(j, t, q) * d²vᵢdqⱼdqₖ(i, l, m, t, q, p))
end

@inline function d²gᵏdqₗdpₘ(k::Val, l::Val, m::Val, t, q, p)
    i, j = constraint_indices(k)
    dbᵢdxⱼ(i, l, t, q) * dvᵢdpⱼ(j, m, t, q, p) - dbᵢdxⱼ(j, l, t, q) * dvᵢdpⱼ(i, m, t, q, p)
end

# ∂²gᵏ/∂pₗ∂pₘ vanishes identically — `gᵏ` is linear in `p`. Every contraction below therefore drops
# the terms it would contribute rather than multiplying by an exact zero, which also keeps a `NaN`
# elsewhere in the expression from spreading through a `0 * Inf`.


"""
    g₁(t, q, p)
    g₂(t, q, p)
    g₃(t, q, p)

The three constraints under the labelling of Li, Zhang & Liu, `gᵏ` evaluated at a fixed `k`. Note
that `g₁` and `g₂` used to name the members of the one implemented pair, which were `g³` and `g¹` in
this labelling; all three vanish on the constraint manifold either way.
"""
g₁(t, q, p) = gᵏ(Val(1), t, q, p)
g₂(t, q, p) = gᵏ(Val(2), t, q, p)
g₃(t, q, p) = gᵏ(Val(3), t, q, p)

# The constraints do not depend on the parameters, but `GeometricSolutions.compute_invariant` calls
# what it is given as `f(t, q, p, params)`, so these are the forms that can be passed to it — as
# `docs/src/examples/iter_cylindrical.md` does.
g₁(t, q, p, params) = g₁(t, q, p)
g₂(t, q, p, params) = g₂(t, q, p)
g₃(t, q, p, params) = g₃(t, q, p)
