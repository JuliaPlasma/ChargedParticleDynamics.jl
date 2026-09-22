@doc raw"""
The field quantities a model reads, evaluated once at a point of its phasespace.

`ElectromagneticFields` answers each quantity as a tensor-valued generic taking the field first —
`DA♭(field, t, ξ)` is the whole matrix ``\partial_j A_i`` — while the equations of this package are
written component by component, `dA₁dx₂(t, q)`, in the notation of the papers they come from.
A [`FieldPoint`](@ref) joins the two. It evaluates each tensor a model needs once, and the scalar
accessors defined here read one entry each.

A `FieldPoint` goes where the phasespace point `q` would. The accessors dispatch on it, and
indexing it indexes the point it was built from, so an expression that reads both the field and
`q[4]` works unchanged. Every right-hand side builds one on entry, from `params.field`, and hands it
on; the functions it calls never see the field itself. See [`fieldpoint`](@ref).

The accessors ignore their `t`: a `FieldPoint` holds the field at the time it was built for, and
every caller builds it at the time it then passes on.
"""
module FieldPoints

import ElectromagneticFields as EMF
using StaticArrays: SVector

export FieldPoint, fieldpoint

"""
    FieldPoint{T}

A phasespace point `q` together with the field tensors evaluated at its first three components.
Built by [`fieldpoint`](@ref). It is an `AbstractVector` whose entries are those of `q`, so it can
be indexed, multiplied and broadcast wherever `q` can.
"""
struct FieldPoint{T, Q <: AbstractVector{T}, V <: NamedTuple} <: AbstractVector{T}
    q::Q
    values::V
end

Base.size(P::FieldPoint) = size(P.q)
Base.getindex(P::FieldPoint, i::Int) = P.q[i]
Base.IndexStyle(::Type{<:FieldPoint}) = IndexLinear()

# One call per requested tensor, spelled out so that the result is a concrete `NamedTuple`.
@generated function fieldvalues(field, t, x, ::Val{names}) where {names}
    calls = (:(EMF.$(name)(field, t, x)) for name in names)
    :(NamedTuple{$names}(($(calls...),)))
end

"""
    fieldpoint(field, t, q, Val(names))

Evaluate the `ElectromagneticFields` generics `names` — a tuple such as `(:A♭, :DA♭, :g♭)` — once
at `(t, q[1:3])`, and return them with `q` as a [`FieldPoint`](@ref). Each model family declares the
tuple it needs beside its equations.

A `FieldPoint` passed in is returned as it is, so a right-hand side that calls another hands its
point on rather than evaluating the field twice. It must then hold every tensor the callee reads:
an accessor for a tensor that was not requested fails on the missing `NamedTuple` field.
"""
@inline function fieldpoint(field, t, q::AbstractVector, names::Val)
    FieldPoint(q, fieldvalues(field, t, SVector(q[1], q[2], q[3]), names))
end

@inline fieldpoint(field, t, P::FieldPoint, ::Val) = P

#
# The scalar accessors, under the names `ElectromagneticFields` 0.8 generated into each module:
# `A₁` is the first covariant component of `A♭`, `b¹` the first contravariant one of `b♯`,
# `dA₁dx₂` is `DA♭[1, 2]`, `dg¹¹dx₂` is `Dg♯[1, 1, 2]`, and so on. The metric appears only through
# its diagonal, since every chart this package uses is orthogonal.
#

const SUB = ('₁', '₂', '₃')
const SUP = ('¹', '²', '³')

accessor(name, tensor, index...) = quote
    export $name
    @inline $name(t, P::FieldPoint) = P.values.$tensor[$(index...)]
end

for i in 1:3
    eval(accessor(Symbol("A", SUB[i]), :A♭, i))
    eval(accessor(Symbol("B", SUB[i]), :B♭, i))
    eval(accessor(Symbol("b", SUB[i]), :b♭, i))
    eval(accessor(Symbol("b", SUP[i]), :b♯, i))
    eval(accessor(Symbol("E", SUB[i]), :E♭, i))
    eval(accessor(Symbol("dBdx", SUB[i]), :DB, i))
    eval(accessor(Symbol("g", SUB[i], SUB[i]), :g♭, i, i))
    eval(accessor(Symbol("g", SUP[i], SUP[i]), :g♯, i, i))

    for j in 1:3
        eval(accessor(Symbol("dA", SUB[i], "dx", SUB[j]), :DA♭, i, j))
        eval(accessor(Symbol("db", SUB[i], "dx", SUB[j]), :Db♭, i, j))
        eval(accessor(Symbol("dE", SUB[i], "dx", SUB[j]), :DE♭, i, j))
        eval(accessor(Symbol("d²Bdx", SUB[i], "dx", SUB[j]), :DDB, i, j))
        eval(accessor(Symbol("dg", SUB[i], SUB[i], "dx", SUB[j]), :Dg♭, i, i, j))
        eval(accessor(Symbol("dg", SUP[i], SUP[i], "dx", SUB[j]), :Dg♯, i, i, j))

        for k in 1:3
            eval(accessor(
                Symbol("d²A", SUB[i], "dx", SUB[j], "dx", SUB[k]), :DDA♭, i, j, k))
            eval(accessor(
                Symbol("d²b", SUB[i], "dx", SUB[j], "dx", SUB[k]), :DDb♭, i, j, k))
            eval(accessor(
                Symbol("d²g", SUP[i], SUP[i], "dx", SUB[j], "dx", SUB[k]), :DDg♯, i, i, j, k))
        end
    end
end

eval(accessor(:B, :B))
eval(accessor(:φ, :φ))

end
