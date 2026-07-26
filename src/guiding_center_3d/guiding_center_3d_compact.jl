#
# The compact form of the Hamilton-Dirac guiding centre equations, Eq. (29) of Li, Zhang & Liu.
#
# Eq. (11) — what `hodeproblem` builds — carries `bₘ` in the denominator of both Lagrange
# multipliers, with `m` the index of the constraint the pair omits. Rearranged, the right-hand side
# splits into a part that is regular there and a part proportional to the constraints `gʲ`, which
# vanish along the solution. Dropping the latter is Eq. (29), and removes all of the singularity bar
# one factor of `(pₘ-Aₘ)/bₘ`; the last paragraph below is about removing that one too.
#
# The paper writes that rearrangement with Cartesian vector identities (`∇×b`, `b×∇B`) and
# `H = ½Σ(pᵢ-Aᵢ)²`, which does not carry over: eight of the thirteen equilibria here are curvilinear
# and the Hamiltonian of this package carries `g¹¹, g²², g³³`. It is derived below instead from the
# objects the model already has, which makes it coordinate-general and, as it turns out, shorter.
#
# Write the constraints in the antisymmetric labelling `c = b × v`, so `cₗ = ε_{lab} bₐ v_b`, related
# to the paper's by `c₁ = g²`, `c₂ = -g¹`, `c₃ = g³`. Then
#
#     ∂cₖ/∂pₗ = ε_{kal} bₐ                                                                    (i)
#
# with no metric in it. For the pair retaining `(cᵢ, cⱼ)`, the cyclic complement of `m`, the
# multiplier contribution to `Ẋ` is
#
#     Cₗ = Σₖ λᵏ ∂cₖ/∂pₗ = -bₘ {cₗ, H} / {cᵢ, cⱼ} ,                                          (ii)
#
# uniformly in `l`: for `l ∈ {i, j}` this is exact algebra from (i), and for `l = m` it needs
# `Σₗ bₗ cₗ = b·(b×v) = 0`, which holds identically and therefore gives `Σₗ bₗ {cₗ, H} = 0` on the
# constraint manifold. That identity is the one thing Eq. (29) uses the constraints for.
#
# `{cᵢ, cⱼ} = bₘ D` with `D` the same for all three `m`, so the singular factor `bₘ / {cᵢ, cⱼ}`
# may be replaced by `1/D` with
#
#     D = Σₘ bₘ {cᵢ, cⱼ}(m) / Σₘ bₘ bₘ ,                                                    (iii)
#
# whose denominator cannot vanish, since `|b| = 1`. What is left is free of every `bₘ` denominator
# and no longer depends on which pair was chosen.
#
# The momentum equation follows the same way. Splitting `∂cₖ/∂qₗ` into its `∂b/∂x` and `∂A/∂x` parts,
# the second contracts straight back into `C` through (i), and the first does too once `vₐ = u bₐ` is
# used on the manifold:
#
#     ṗₗ = -∂H/∂qₗ + Σₐ (∂Aₐ/∂xₗ + u ∂bₐ/∂xₗ) Cₐ .                                          (iv)
#
# `u` is the parallel velocity. Eq. (29) as printed uses `(pₘ-Aₘ)/bₘ` for it, which is `0/0` exactly
# where the formulation was singular to begin with; `u(t, q, p)` is the same quantity on the manifold
# — `vₐ = bₐ u` and `|b| = 1` give `Σₐ vₐ gᵃᵃ bₐ = u` — and is regular everywhere. Passing
# `constraints = :parallel`, the default, uses it together with (iii). Passing one of the pair
# symbols reproduces the paper literally, singularity included, which is what
# `scripts/study_guiding_center_3d_conditioning.jl` measures.
#

export hodeproblem_compact


# `cₗ = σ gᵏ`: the antisymmetric labelling in terms of the paper's.
cross_constraint(::Val{1}) = (Val(2), +1)
cross_constraint(::Val{2}) = (Val(1), -1)
cross_constraint(::Val{3}) = (Val(3), +1)

# The two indices completing `m` to a positively oriented triple.
cyclic_complement(::Val{1}) = (Val(2), Val(3))
cyclic_complement(::Val{2}) = (Val(3), Val(1))
cyclic_complement(::Val{3}) = (Val(1), Val(2))


"""
    cₗ(::Val{l}, t, q, p)

The `l`-th component of `b × v`, which is the `l`-th constraint in the antisymmetric labelling.
Related to the `gᵏ` of Li, Zhang & Liu by `c₁ = g²`, `c₂ = -g¹`, `c₃ = g³`, so the two sets vanish
together.
"""
@inline function cₗ(l::Val, t, q, p)
    k, σ = cross_constraint(l)
    σ * gᵏ(k, t, q, p)
end

# {cₗ, H}
@inline function bracket_cH(l::Val, t, q, p, params)
    k, σ = cross_constraint(l)
    σ * bracket_gH(k, t, q, p, params)
end

# {cᵢ, cⱼ} over the cyclic complement of `m`, which equals `bₘ [B + (p-A)·(∇×b)]`.
@inline function bracket_cc(m::Val, t, q, p)
    i, j = cyclic_complement(m)
    kᵢ, σᵢ = cross_constraint(i)
    kⱼ, σⱼ = cross_constraint(j)
    σᵢ * σⱼ * bracket_gg(kᵢ, kⱼ, t, q, p)
end


"""
    compact_denominator(t, q, p)

The scalar `D = B + (p-A)·(∇×b)` of Eq. (29), obtained as `Σₘ bₘ {cᵢ, cⱼ}(m) / Σₘ bₘ bₘ` so that no
single component of `b` ever divides. The three ratios `{cᵢ, cⱼ}(m) / bₘ` agree on the constraint
manifold; this weighted mean of them is the combination that stays finite off it as well.
"""
compact_denominator(t, q, p) = contract(m -> bᵢ(m, t, q) * bracket_cc(m, t, q, p)) /
                               contract(m -> bᵢ(m, t, q)^2)


# The scalar the three brackets `{cₗ, H}` are divided by, Eq. (ii) above: regularised for
# `:parallel`, literal for a constraint pair whose omitted index is `m`.
@inline compact_scale(t, q, p, ::Nothing) = -inv(compact_denominator(t, q, p))
@inline compact_scale(t, q, p, m::Val) = -bᵢ(m, t, q) / bracket_cc(m, t, q, p)

"""
    compact_multipliers(t, q, p, params, m)

The three components of `Cₗ`, the multiplier contribution to `Ẋ`, as a tuple. All three share the
same denominator, and it is the expensive part of the expression — nine Poisson-bracket terms — so
they are computed together rather than one at a time.
"""
@inline function compact_multipliers(t, q, p, params, m)
    s = compact_scale(t, q, p, m)
    (s * bracket_cH(Val(1), t, q, p, params),
     s * bracket_cH(Val(2), t, q, p, params),
     s * bracket_cH(Val(3), t, q, p, params))
end

# Index a tuple with a compile-time index, so that it stays one.
@inline component(x::Tuple, ::Val{i}) where {i} = x[i]

# The parallel velocity entering Eq. (iv), likewise dispatching on the choice of form.
@inline compact_parallel_velocity(t, q, p, ::Nothing) = u(t, q, p)
@inline compact_parallel_velocity(t, q, p, m::Val{i}) where {i} = (p[i] - Aᵢ(m, t, q)) / bᵢ(m, t, q)


"""
    compact_index(constraints::Symbol)

The index of the constraint a pair omits, which is the component of `b` the multipliers of that pair
divide by: `:g31 → 1`, `:g23 → 2`, `:g12 → 3`. `:parallel` maps to `nothing`, selecting the
regularised, pair-independent form of the compact equations.
"""
compact_index(::Val{:parallel}) = nothing
compact_index(::Val{:g31}) = Val(1)
compact_index(::Val{:g23}) = Val(2)
compact_index(::Val{:g12}) = Val(3)

compact_index(constraints::Symbol) = compact_index(Val(constraints))


function guiding_center_3d_compact_v(v, t, q, p, params, m = nothing)
    C = compact_multipliers(t, q, p, params, m)

    components!(v, i -> dHdpᵢ(i, t, q, p, params) + component(C, i))
    nothing
end

function guiding_center_3d_compact_f(f, t, q, p, params, m = nothing)
    C = compact_multipliers(t, q, p, params, m)
    uₚ = compact_parallel_velocity(t, q, p, m)

    components!(f, l -> -dHdqᵢ(l, t, q, p, params) +
        contract(a -> (dAᵢdxⱼ(a, l, t, q) + uₚ * dbᵢdxⱼ(a, l, t, q)) * component(C, a)))
    nothing
end


"""
    hodeproblem_compact(q₀, p₀; kwargs...)
    hodeproblem_compact(x₀ = qᵢ; kwargs...)
    hodeproblem_compact(ics::NamedTuple; kwargs...)

The compact form of the guiding centre system, Eq. (29) of Li, Zhang & Liu, which drops the terms
proportional to the constraints from the right-hand side of [`hodeproblem`](@ref). The two agree on
the constraint manifold. Same three argument forms as `hodeproblem`.

`constraints` defaults to `:parallel`, the derivation at the top of
`src/guiding_center_3d/guiding_center_3d_compact.jl`: coordinate-general, independent of the
constraint pair, and free of the `bₘ = 0` singularity, which is the point of this formulation.
Passing `:g31`, `:g12` or `:g23` instead reproduces the paper's expressions literally for that pair,
`bₘ` in the denominator included.
"""
function hodeproblem_compact(q₀::AbstractVector, p₀::AbstractVector; timespan = DEFAULT_TIMESPAN,
                             timestep = DEFAULT_TIMESTEP, parameters = default_parameters(), periodic = true,
                             constraints = :parallel)
    m = compact_index(constraints)

    _v(v, t, q, p, params) = guiding_center_3d_compact_v(v, t, q, p, params, m)
    _f(f, t, q, p, params) = guiding_center_3d_compact_f(f, t, q, p, params, m)

    HODEProblem(
        _v,
        _f,
        hamiltonian,
        timespan, timestep, q₀, p₀;
        parameters=parameters,
        periodicity=guiding_center_3d_periodicity(q₀, periodic))
end

function hodeproblem_compact(x₀::AbstractVector = qᵢ; timespan = DEFAULT_TIMESPAN, kwargs...)
    ics = initial_conditions(timespan[begin], x₀)
    hodeproblem_compact(ics.q, ics.p; timespan = timespan, kwargs...)
end

hodeproblem_compact(ics::NamedTuple; kwargs...) =
    hodeproblem_compact(ics.q, ics.p; parameters = ics.params, kwargs...)
