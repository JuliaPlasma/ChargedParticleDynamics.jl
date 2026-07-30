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
# `+b₃` for the pair `(g¹, g²)`, `+b₁` for `(g³, g¹)`, `-b₂` for `(g², g³)` — so which pair is usable
# is a property of the equilibrium and of where the orbit sits in it. The sign varies with the pair's
# ordering and does not matter; see `constraint_pair` below.
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
const INDEX_SUPERSCRIPTS = ('¹', '²', '³')

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


#
# The injected field functions evaluated once per point.
#
# Everything below — the constraints, their derivatives, the Hamiltonian gradients, the Poisson
# brackets — is written against the accessors above, and each of them calls straight through to an
# `ElectromagneticFields.@code()`-injected function. Those are the expensive part, by a wide margin:
# `db₁dx₁` for the ITER Solov'ev X-point costs 57 ns against a `g¹¹` that folds to a constant, and a
# single evaluation of the Hamilton-Dirac right-hand side used to call `dbᵢdxⱼ` more than seventy
# times, because every one of the twelve Poisson bracket terms re-entered it from the top.
#
# It used to be worse. On `ElectromagneticFields` 0.6.2 that same `db₁dx₁` was 1905 statements
# containing 108 separate evaluations of `log(x₁)` and cost 549 ns, because `convert(Expr, ::Basic)`
# wrote out SymEngine's expanded tree verbatim and SymEngine shares nothing. 0.6.3 eliminates common
# subexpressions in what it emits, which is why the figure above is 57 and not 549; `Project.toml`
# requires it. The values are unchanged — every field function this package injects is bit-identical
# between the two versions.
#
# `FieldValues` holds one evaluation of each. The functions that consume it need no changes at all:
# they are already generic in their `q` argument, so passing a `FieldValues` where the coordinate
# vector would go picks up the methods defined below by dispatch, and every field read becomes a
# tuple index. That took `guiding_center_3d_v` from 24.7 µs to 3.6 µs — 6.9x, and 32x once 0.6.3 is
# counted too, at 770 ns. The values are identical, so the trajectories are as well, bit for bit.
#
# It is parametric because the integrator evaluates the right-hand side on `ForwardDiff.Dual` as well
# as `Float64` — four of the roughly twenty-two calls per step, for the Newton Jacobian.
#

"""
    FieldValues{T}

One evaluation of every injected field function the first-derivative right-hand sides read: the
vector potential `A`, the unit vector `b`, the diagonal metric `g`, the gradient of `B`, the electric
field `E`, and the first derivatives of `A`, `b` and `g`.

Built by [`fieldvalues`](@ref) and passed where those right-hand sides expect the coordinate vector
`q` — they are generic in that argument, so the accessors pick this up by dispatch and every field
read becomes a tuple index. [`SecondFieldValues`](@ref) adds the second derivatives that only the
canonicalised formulation needs.
"""
struct FieldValues{T}
    A ::NTuple{3,T}
    b ::NTuple{3,T}
    g ::NTuple{3,T}
    dB::NTuple{3,T}
    E ::NTuple{3,T}
    dA::NTuple{3,NTuple{3,T}}
    db::NTuple{3,NTuple{3,T}}
    dg::NTuple{3,NTuple{3,T}}
end

# An expression that names an injected field function directly, rather than going through one of the
# accessors below, reaches the two-argument `f(t, ξ)` form and indexes its second argument. That is a
# missing cached method, not a bug in the caller, and the `MethodError` it raises names `getindex`
# rather than the function that is missing one — so say which.
Base.getindex(::FieldValues, i::Integer) = error(
    "a `FieldValues` was passed to an `ElectromagneticFields` field function as a coordinate " *
    "vector. Some expression on the right-hand side names an injected function that has no " *
    "`FieldValues` method; add one beside the others in `guiding_center_3d_constraints.jl`.")

"""
    fieldvalues(t, q)

Every injected field function the first-derivative right-hand sides need, evaluated at `(t, q)`.
Pass the result where those right-hand sides expect `q`; see [`FieldValues`](@ref).
"""
@inline function fieldvalues(t, q)
    # The generated functions return an exact `0` for a vanishing component — an `Int` literal in a
    # cartesian chart, where `dg¹¹dx₁` and `E₁` fold away entirely — so the tuples need converting
    # rather than merely collecting, or `FieldValues` would not be concretely typed.
    T = eltype(q)
    c(x) = convert(T, x)

    FieldValues{T}(
        (c(A₁(t, q)),     c(A₂(t, q)),     c(A₃(t, q))),
        (c(b₁(t, q)),     c(b₂(t, q)),     c(b₃(t, q))),
        (c(g¹¹(t, q)),    c(g²²(t, q)),    c(g³³(t, q))),
        (c(dBdx₁(t, q)),  c(dBdx₂(t, q)),  c(dBdx₃(t, q))),
        (c(E₁(t, q)),     c(E₂(t, q)),     c(E₃(t, q))),
        ((c(dA₁dx₁(t, q)),  c(dA₁dx₂(t, q)),  c(dA₁dx₃(t, q))),
         (c(dA₂dx₁(t, q)),  c(dA₂dx₂(t, q)),  c(dA₂dx₃(t, q))),
         (c(dA₃dx₁(t, q)),  c(dA₃dx₂(t, q)),  c(dA₃dx₃(t, q)))),
        ((c(db₁dx₁(t, q)),  c(db₁dx₂(t, q)),  c(db₁dx₃(t, q))),
         (c(db₂dx₁(t, q)),  c(db₂dx₂(t, q)),  c(db₂dx₃(t, q))),
         (c(db₃dx₁(t, q)),  c(db₃dx₂(t, q)),  c(db₃dx₃(t, q)))),
        ((c(dg¹¹dx₁(t, q)), c(dg¹¹dx₂(t, q)), c(dg¹¹dx₃(t, q))),
         (c(dg²²dx₁(t, q)), c(dg²²dx₂(t, q)), c(dg²²dx₃(t, q))),
         (c(dg³³dx₁(t, q)), c(dg³³dx₂(t, q)), c(dg³³dx₃(t, q)))))
end

# The cached counterparts of the accessors above, and of the injected names the Hamiltonian gradients
# call directly. Same signatures with a `FieldValues` in place of the coordinate vector, so that
# every expression written against them works unchanged on either.
#
# Written out per index rather than as `Aᵢ(::Val{i}, t, F::FieldValues) where {i}`, which would be
# ambiguous against the per-index methods above: those fix the `Val`s and leave `q` untyped, so
# neither signature is more specific than the other and dispatch has no way to choose.
for i in 1:3
    @eval @inline Aᵢ(::Val{$i}, t, F::FieldValues) = F.A[$i]
    @eval @inline bᵢ(::Val{$i}, t, F::FieldValues) = F.b[$i]

    # `Aᵢ`/`bᵢ` do not cover every call: `u` and `ϑ` name `A₁`, `b₁` and their siblings directly, and
    # without these they reach the injected functions and index the `FieldValues` as if it were the
    # coordinate vector.
    @eval @inline $(Symbol("A", INDEX_SUBSCRIPTS[i]))(t, F::FieldValues) = F.A[$i]
    @eval @inline $(Symbol("b", INDEX_SUBSCRIPTS[i]))(t, F::FieldValues) = F.b[$i]

    @eval @inline $(Symbol("g", INDEX_SUPERSCRIPTS[i], INDEX_SUPERSCRIPTS[i]))(t, F::FieldValues) = F.g[$i]
    @eval @inline $(Symbol("dBdx", INDEX_SUBSCRIPTS[i]))(t, F::FieldValues) = F.dB[$i]
    @eval @inline $(Symbol("E", INDEX_SUBSCRIPTS[i]))(t, F::FieldValues) = F.E[$i]

    for j in 1:3
        @eval @inline dAᵢdxⱼ(::Val{$i}, ::Val{$j}, t, F::FieldValues) = F.dA[$i][$j]
        @eval @inline dbᵢdxⱼ(::Val{$i}, ::Val{$j}, t, F::FieldValues) = F.db[$i][$j]

        @eval @inline $(Symbol("dg", INDEX_SUPERSCRIPTS[i], INDEX_SUPERSCRIPTS[i], "dx", INDEX_SUBSCRIPTS[j]))(t, F::FieldValues) = F.dg[$i][$j]
    end
end


#
# The same for the second derivatives, which only the canonicalised formulation needs.
#
# `guiding_center_3d_canonical.jl` differentiates the Lagrange multipliers, so it reaches
# `d²Aᵢdxⱼdxₖ`, `d²bᵢdxⱼdxₖ`, `d²gⁱⁱdxⱼdxₖ`, `d²Bdxᵢdxⱼ` and `dEᵢdxⱼ` — and reaches them far more
# often than it has distinct values to read, because `dλ₁dqⱼ` and `dλ₂dqⱼ` recompute `λₒ` and its
# derivative for each of the three components. `d²b₁dx₁dx₁` costs 80 ns for the ITER Solov'ev
# X-point, and 1733 on `ElectromagneticFields` 0.6.2 where it was 6231 statements, so the repetition
# is the whole cost of that formulation.
#
# Held separately from `FieldValues` rather than merged into it: `hodeproblem` and
# `hodeproblem_compact` never touch a second derivative, and filling these ninety-nine slots would
# cost them more than the first-derivative caching saves.
#
"""
    SecondFieldValues{T,F}

A [`FieldValues`](@ref) together with every second derivative the canonicalised right-hand side
reads: the second derivatives of `A`, `b` and `g`, the Hessian of `B`, and the first derivatives of
`E`.

Built by [`secondfieldvalues`](@ref) and used the same way — it answers everything a `FieldValues`
does, and is held separately rather than merged into it because `hodeproblem` and
`hodeproblem_compact` never touch a second derivative.
"""
struct SecondFieldValues{T,F<:FieldValues{T}}
    first::F
    d²A::NTuple{3,NTuple{3,NTuple{3,T}}}
    d²b::NTuple{3,NTuple{3,NTuple{3,T}}}
    d²g::NTuple{3,NTuple{3,NTuple{3,T}}}
    d²B::NTuple{3,NTuple{3,T}}
    dE ::NTuple{3,NTuple{3,T}}
end

Base.getindex(::SecondFieldValues, i::Integer) = error(
    "a `SecondFieldValues` was passed to an `ElectromagneticFields` field function as a coordinate " *
    "vector. Some expression on the right-hand side names an injected function that has no " *
    "`SecondFieldValues` method; add one beside the others in `guiding_center_3d_constraints.jl`.")

"""
    secondfieldvalues(t, q)

[`fieldvalues`](@ref) together with every second derivative the canonicalised right-hand side needs.
Pass it where that right-hand side expects `q`; it answers everything a `FieldValues` does.

The mixed derivatives are stored for all twenty-seven index triples rather than for the eighteen the
symmetry of `∂²/∂xⱼ∂xₖ` leaves independent. `ElectromagneticFields` expands each one separately, so
`d²b₁dx₁dx₂` and `d²b₁dx₂dx₁` are different floating-point expressions of the same quantity and need
not agree in the last bit; reading one for the other would change the trajectory.
"""
@inline function secondfieldvalues(t, q)
    F = fieldvalues(t, q)
    T = eltype(q)
    c(x) = convert(T, x)

    SecondFieldValues{T,typeof(F)}(F,
        (((c(d²A₁dx₁dx₁(t, q)), c(d²A₁dx₁dx₂(t, q)), c(d²A₁dx₁dx₃(t, q))),
          (c(d²A₁dx₂dx₁(t, q)), c(d²A₁dx₂dx₂(t, q)), c(d²A₁dx₂dx₃(t, q))),
          (c(d²A₁dx₃dx₁(t, q)), c(d²A₁dx₃dx₂(t, q)), c(d²A₁dx₃dx₃(t, q)))),
         ((c(d²A₂dx₁dx₁(t, q)), c(d²A₂dx₁dx₂(t, q)), c(d²A₂dx₁dx₃(t, q))),
          (c(d²A₂dx₂dx₁(t, q)), c(d²A₂dx₂dx₂(t, q)), c(d²A₂dx₂dx₃(t, q))),
          (c(d²A₂dx₃dx₁(t, q)), c(d²A₂dx₃dx₂(t, q)), c(d²A₂dx₃dx₃(t, q)))),
         ((c(d²A₃dx₁dx₁(t, q)), c(d²A₃dx₁dx₂(t, q)), c(d²A₃dx₁dx₃(t, q))),
          (c(d²A₃dx₂dx₁(t, q)), c(d²A₃dx₂dx₂(t, q)), c(d²A₃dx₂dx₃(t, q))),
          (c(d²A₃dx₃dx₁(t, q)), c(d²A₃dx₃dx₂(t, q)), c(d²A₃dx₃dx₃(t, q))))),
        (((c(d²b₁dx₁dx₁(t, q)), c(d²b₁dx₁dx₂(t, q)), c(d²b₁dx₁dx₃(t, q))),
          (c(d²b₁dx₂dx₁(t, q)), c(d²b₁dx₂dx₂(t, q)), c(d²b₁dx₂dx₃(t, q))),
          (c(d²b₁dx₃dx₁(t, q)), c(d²b₁dx₃dx₂(t, q)), c(d²b₁dx₃dx₃(t, q)))),
         ((c(d²b₂dx₁dx₁(t, q)), c(d²b₂dx₁dx₂(t, q)), c(d²b₂dx₁dx₃(t, q))),
          (c(d²b₂dx₂dx₁(t, q)), c(d²b₂dx₂dx₂(t, q)), c(d²b₂dx₂dx₃(t, q))),
          (c(d²b₂dx₃dx₁(t, q)), c(d²b₂dx₃dx₂(t, q)), c(d²b₂dx₃dx₃(t, q)))),
         ((c(d²b₃dx₁dx₁(t, q)), c(d²b₃dx₁dx₂(t, q)), c(d²b₃dx₁dx₃(t, q))),
          (c(d²b₃dx₂dx₁(t, q)), c(d²b₃dx₂dx₂(t, q)), c(d²b₃dx₂dx₃(t, q))),
          (c(d²b₃dx₃dx₁(t, q)), c(d²b₃dx₃dx₂(t, q)), c(d²b₃dx₃dx₃(t, q))))),
        (((c(d²g¹¹dx₁dx₁(t, q)), c(d²g¹¹dx₁dx₂(t, q)), c(d²g¹¹dx₁dx₃(t, q))),
          (c(d²g¹¹dx₂dx₁(t, q)), c(d²g¹¹dx₂dx₂(t, q)), c(d²g¹¹dx₂dx₃(t, q))),
          (c(d²g¹¹dx₃dx₁(t, q)), c(d²g¹¹dx₃dx₂(t, q)), c(d²g¹¹dx₃dx₃(t, q)))),
         ((c(d²g²²dx₁dx₁(t, q)), c(d²g²²dx₁dx₂(t, q)), c(d²g²²dx₁dx₃(t, q))),
          (c(d²g²²dx₂dx₁(t, q)), c(d²g²²dx₂dx₂(t, q)), c(d²g²²dx₂dx₃(t, q))),
          (c(d²g²²dx₃dx₁(t, q)), c(d²g²²dx₃dx₂(t, q)), c(d²g²²dx₃dx₃(t, q)))),
         ((c(d²g³³dx₁dx₁(t, q)), c(d²g³³dx₁dx₂(t, q)), c(d²g³³dx₁dx₃(t, q))),
          (c(d²g³³dx₂dx₁(t, q)), c(d²g³³dx₂dx₂(t, q)), c(d²g³³dx₂dx₃(t, q))),
          (c(d²g³³dx₃dx₁(t, q)), c(d²g³³dx₃dx₂(t, q)), c(d²g³³dx₃dx₃(t, q))))),
        ((c(d²Bdx₁dx₁(t, q)), c(d²Bdx₁dx₂(t, q)), c(d²Bdx₁dx₃(t, q))),
         (c(d²Bdx₂dx₁(t, q)), c(d²Bdx₂dx₂(t, q)), c(d²Bdx₂dx₃(t, q))),
         (c(d²Bdx₃dx₁(t, q)), c(d²Bdx₃dx₂(t, q)), c(d²Bdx₃dx₃(t, q)))),
        ((c(dE₁dx₁(t, q)), c(dE₁dx₂(t, q)), c(dE₁dx₃(t, q))),
         (c(dE₂dx₁(t, q)), c(dE₂dx₂(t, q)), c(dE₂dx₃(t, q))),
         (c(dE₃dx₁(t, q)), c(dE₃dx₂(t, q)), c(dE₃dx₃(t, q)))))
end

# Building a cache from a cache is the identity, which is what lets the canonicalised right-hand side
# hand its `SecondFieldValues` straight to `guiding_center_3d_v` instead of having it fill a second
# one from the coordinate vector it no longer has.
@inline fieldvalues(t, F::FieldValues) = F
@inline fieldvalues(t, S::SecondFieldValues) = S

for i in 1:3
    @eval @inline Aᵢ(::Val{$i}, t, S::SecondFieldValues) = S.first.A[$i]
    @eval @inline bᵢ(::Val{$i}, t, S::SecondFieldValues) = S.first.b[$i]

    @eval @inline $(Symbol("A", INDEX_SUBSCRIPTS[i]))(t, S::SecondFieldValues) = S.first.A[$i]
    @eval @inline $(Symbol("b", INDEX_SUBSCRIPTS[i]))(t, S::SecondFieldValues) = S.first.b[$i]

    @eval @inline $(Symbol("g", INDEX_SUPERSCRIPTS[i], INDEX_SUPERSCRIPTS[i]))(t, S::SecondFieldValues) = S.first.g[$i]
    @eval @inline $(Symbol("dBdx", INDEX_SUBSCRIPTS[i]))(t, S::SecondFieldValues) = S.first.dB[$i]
    @eval @inline $(Symbol("E", INDEX_SUBSCRIPTS[i]))(t, S::SecondFieldValues) = S.first.E[$i]

    for j in 1:3
        @eval @inline dAᵢdxⱼ(::Val{$i}, ::Val{$j}, t, S::SecondFieldValues) = S.first.dA[$i][$j]
        @eval @inline dbᵢdxⱼ(::Val{$i}, ::Val{$j}, t, S::SecondFieldValues) = S.first.db[$i][$j]

        @eval @inline $(Symbol("dg", INDEX_SUPERSCRIPTS[i], INDEX_SUPERSCRIPTS[i], "dx", INDEX_SUBSCRIPTS[j]))(t, S::SecondFieldValues) = S.first.dg[$i][$j]

        @eval @inline $(Symbol("d²Bdx", INDEX_SUBSCRIPTS[i], "dx", INDEX_SUBSCRIPTS[j]))(t, S::SecondFieldValues) = S.d²B[$i][$j]
        @eval @inline $(Symbol("dE", INDEX_SUBSCRIPTS[i], "dx", INDEX_SUBSCRIPTS[j]))(t, S::SecondFieldValues) = S.dE[$i][$j]

        for k in 1:3
            @eval @inline d²Aᵢdxⱼdxₖ(::Val{$i}, ::Val{$j}, ::Val{$k}, t, S::SecondFieldValues) = S.d²A[$i][$j][$k]
            @eval @inline d²bᵢdxⱼdxₖ(::Val{$i}, ::Val{$j}, ::Val{$k}, t, S::SecondFieldValues) = S.d²b[$i][$j][$k]

            @eval @inline $(Symbol("d²g", INDEX_SUPERSCRIPTS[i], INDEX_SUPERSCRIPTS[i], "dx", INDEX_SUBSCRIPTS[j], "dx", INDEX_SUBSCRIPTS[k]))(t, S::SecondFieldValues) = S.d²g[$i][$j][$k]
        end
    end
end


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
    unval(::Val{i})

The integer a `Val` index carries. `constraint_pair` and `compact_index` hand out `Val`s so that the
right-hand side stays specialised on them, and nothing inside the right-hand side ever needs the
integer back; this is for the callers outside it — the test suite and the scripts — that have to index
a tuple of data series by a constraint index.

Deliberately not defined for the `nothing` that `compact_index(:parallel)` returns: `:parallel` retains
no pair, so there is no index to unwrap and asking for one is a mistake rather than a special case.
"""
unval(::Val{i}) where {i} = i


"""
    constraint_pair(constraints::Symbol)

The two constraints selected by the `constraints` keyword of the problem constructors, in the order
they enter the multipliers:

| symbol | pair | multiplier denominator `{g₁, g₂}` |
|---|---|---|
| `:g31` | `(g³, g¹)` | `+b₁ [B + (p-A)·(∇×b)]` |
| `:g12` | `(g¹, g²)` | `+b₃ [B + (p-A)·(∇×b)]` |
| `:g23` | `(g², g³)` | `-b₂ [B + (p-A)·(∇×b)]` |

The sign is not uniform, and the minus on `:g23` is not a typo: in the antisymmetric labelling
`c₁ = g²`, `c₂ = -g¹`, `c₃ = g³` the bracket over a cyclic pair is `+bₘ [ … ]`, and `(g², g³)` is the
one of the three whose ordering here reverses that cycle — `{g², g³} = {c₁, c₃} = -{c₃, c₁}`. It makes
no difference to the dynamics, since reversing a pair's order flips `λₒ` and both multipliers together
and the products `λ₁ ∂g₁ + λ₂ ∂g₂` are invariant, but it does mean `λₒ` and `bₘ` need not share a sign.
The measured values in `docs/src/findings.md` show it.

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
