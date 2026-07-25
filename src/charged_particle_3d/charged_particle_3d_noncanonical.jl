import GeometricEquations
using LinearAlgebra: I
using GeometricEquations: ODEProblem, IODEProblem, LODEProblem, SODEProblem
using GeometricSolutions: GeometricSolution, DataSeries, ScalarDataSeries, TimeSeries
using GeometricSolutions: compute_invariant, compute_invariant_error


const Δt = 0.01
const tspan = (0.0, 10.0)

ϑ₁(t, q) = g₁₁(t,q) * q[4] + A₁(t,q)
ϑ₂(t, q) = g₂₂(t,q) * q[5] + A₂(t,q)
ϑ₃(t, q) = g₃₃(t,q) * q[6] + A₃(t,q)

# Derivatives of the one-form with respect to the *position* coordinates. These were used by `ω`
# and `dϑ` but defined nowhere, which stayed hidden only because neither function was reachable:
# `ω` had its output argument in the wrong slot, so the `LODEProblem` built from it could never be
# evaluated. The metric depends on position, hence the gᵢᵢ derivative terms.
dϑ₁dx₁(t, q) = dg₁₁dx₁(t,q) * q[4] + dA₁dx₁(t,q)
dϑ₁dx₂(t, q) = dg₁₁dx₂(t,q) * q[4] + dA₁dx₂(t,q)
dϑ₁dx₃(t, q) = dg₁₁dx₃(t,q) * q[4] + dA₁dx₃(t,q)

dϑ₂dx₁(t, q) = dg₂₂dx₁(t,q) * q[5] + dA₂dx₁(t,q)
dϑ₂dx₂(t, q) = dg₂₂dx₂(t,q) * q[5] + dA₂dx₂(t,q)
dϑ₂dx₃(t, q) = dg₂₂dx₃(t,q) * q[5] + dA₂dx₃(t,q)

dϑ₃dx₁(t, q) = dg₃₃dx₁(t,q) * q[6] + dA₃dx₁(t,q)
dϑ₃dx₂(t, q) = dg₃₃dx₂(t,q) * q[6] + dA₃dx₂(t,q)
dϑ₃dx₃(t, q) = dg₃₃dx₃(t,q) * q[6] + dA₃dx₃(t,q)


function ϑ(θ, t, q)
    θ[1] = ϑ₁(t,q)
    θ[2] = ϑ₂(t,q)
    θ[3] = ϑ₃(t,q)
    θ[4] = zero(eltype(q))
    θ[5] = zero(eltype(q))
    θ[6] = zero(eltype(q))
    nothing
end


@doc raw"""
The symplectic two-form of the noncanonical formulation,
``\Omega_{ij} = \partial \vartheta_{i} / \partial z^{j} - \partial \vartheta_{j} / \partial z^{i}``
on the phase space ``z = (x, v)``, which is the convention shared with the 4D guiding centre and
with `GeometricProblems`.

The velocity block is ``\partial \vartheta_{i} / \partial v^{j} = g_{ij}``, the metric — not the
identity, which is what it reduces to in cartesian coordinates only.
"""
function ω(Β, t, q, params)
    Β .= 0

    Β[1,2] = dϑ₁dx₂(t,q) - dϑ₂dx₁(t,q)
    Β[1,3] = dϑ₁dx₃(t,q) - dϑ₃dx₁(t,q)
    Β[1,4] = g₁₁(t,q)

    Β[2,1] = dϑ₂dx₁(t,q) - dϑ₁dx₂(t,q)
    Β[2,3] = dϑ₂dx₃(t,q) - dϑ₃dx₂(t,q)
    Β[2,5] = g₂₂(t,q)

    Β[3,1] = dϑ₃dx₁(t,q) - dϑ₁dx₃(t,q)
    Β[3,2] = dϑ₃dx₂(t,q) - dϑ₂dx₃(t,q)
    Β[3,6] = g₃₃(t,q)

    Β[4,1] = -g₁₁(t,q)
    Β[5,2] = -g₂₂(t,q)
    Β[6,3] = -g₃₃(t,q)

    nothing
end

# LODE calls the two-form with an additional velocity slot; ω depends only on q.
ω(Β, t, q, v, params) = ω(Β, t, q, params)


@doc raw"""
The Jacobian of the one-form on the phase space ``z = (x, v)``,
``(\mathrm{d}\vartheta)_{ij} = \partial \vartheta_{i} / \partial z^{j}``.

The velocity block is ``\partial \vartheta_{i} / \partial v^{j} = g_{ij}``; it was previously
set to zero, which is wrong in every coordinate system, cartesian included.
"""
function dϑ(dϑ, t, q, params)
    dϑ .= 0

    dϑ[1,1] = dϑ₁dx₁(t,q)
    dϑ[1,2] = dϑ₁dx₂(t,q)
    dϑ[1,3] = dϑ₁dx₃(t,q)
    dϑ[1,4] = g₁₁(t,q)

    dϑ[2,1] = dϑ₂dx₁(t,q)
    dϑ[2,2] = dϑ₂dx₂(t,q)
    dϑ[2,3] = dϑ₂dx₃(t,q)
    dϑ[2,5] = g₂₂(t,q)

    dϑ[3,1] = dϑ₃dx₁(t,q)
    dϑ[3,2] = dϑ₃dx₂(t,q)
    dϑ[3,3] = dϑ₃dx₃(t,q)
    dϑ[3,6] = g₃₃(t,q)

    nothing
end


β₁(t,q) = dϑ₃dx₂(t,q) - dϑ₂dx₃(t,q)
β₂(t,q) = dϑ₁dx₃(t,q) - dϑ₃dx₁(t,q)
β₃(t,q) = dϑ₂dx₁(t,q) - dϑ₁dx₂(t,q)

β(t,q) = sqrt(β₁(t,q)^2 + β₂(t,q)^2 + β₃(t,q)^2)


hamiltonian(t,q) = 0.5 * ( g₁₁(t,q) * q[4]^2 + g₂₂(t,q) * q[5]^2 + g₃₃(t,q) * q[6]^2 ) + φ(t,q)
hamiltonian(t,q,params) = hamiltonian(t,q)
hamiltonian(t,q,p,params) = hamiltonian(t,q)
lagrangian(t,q,v,params) = ϑ₁(t, q) * v[1] + ϑ₂(t, q) * v[2] + ϑ₃(t, q) * v[3] - hamiltonian(t,q)


dHdx₁(t, q) = dφdx₁(t,q)
dHdx₂(t, q) = dφdx₂(t,q)
dHdx₃(t, q) = dφdx₃(t,q)
dHdx₄(t, q) = q[4]
dHdx₅(t, q) = q[5]
dHdx₆(t, q) = q[6]

function dH(dH, t, q)
    dH[1] = dHdx₁(t, q)
    dH[2] = dHdx₂(t, q)
    dH[3] = dHdx₃(t, q)
    dH[4] = dHdx₄(t, q)
    dH[5] = dHdx₅(t, q)
    dH[6] = dHdx₆(t, q)
    nothing
end


v₁(t, q, v) = q[4]
v₂(t, q, v) = q[5]
v₃(t, q, v) = q[6]
v₄(t, q, v) = E₁(t,q) + q[5] * B₃(t,q) - q[6] * B₂(t,q)
v₅(t, q, v) = E₂(t,q) + q[6] * B₁(t,q) - q[4] * B₃(t,q)
v₆(t, q, v) = E₃(t,q) + q[4] * B₂(t,q) - q[5] * B₁(t,q)


# `GeometricEquations` recognises periodicity only as a `(xmin, xmax)` tuple of arrays
# (`hasperiodicity(::GEperType{<:Tuple{AT,AT}})`); a single vector of periods is silently ignored,
# which is what this returned before.
#
# The periodic range of the position coordinates comes from the coordinate system, via the
# `rangemin`/`rangemax` injected with the field code — the same mechanism the guiding centre models
# use. Hard-coding the third coordinate to [0, 2π) would have wrapped `z` in the cartesian
# equilibria once the periodicity started taking effect. The velocities are never periodic.
function charged_particle_3d_periodicity(qᵢ, periodic=true)
    T    = eltype(qᵢ)
    xmin = -T(Inf) * ones(T, size(qᵢ, 1))
    xmax = +T(Inf) * ones(T, size(qᵢ, 1))

    if periodic
        xmin[1:3] .= rangemin(xmin[1:3])
        xmax[1:3] .= rangemax(xmax[1:3])
    end

    return (xmin, xmax)
end


function charged_particle_3d_pᵢ(tᵢ, qᵢ)
    pᵢ = zero(qᵢ)
    ϑ(pᵢ, tᵢ, qᵢ)
    return pᵢ
end


function charged_particle_3d_v(v, t, q, params)
    v[1] = v₁(t, q, v)
    v[2] = v₂(t, q, v)
    v[3] = v₃(t, q, v)
    v[4] = v₄(t, q, v)
    v[5] = v₅(t, q, v)
    v[6] = v₆(t, q, v)

    nothing
end

charged_particle_3d_v(v, t, q, p, params) = charged_particle_3d_v(v, t, q, params)


charged_particle_3d_iode_ϑ(θ, t, q, v, params) = ϑ(θ, t, q)


function charged_particle_3d_iode_f(f, t, q, v, params)
    f[1] = dA₁dx₁(t,q) * v[1] + dA₂dx₁(t,q) * v[2] + dA₃dx₁(t,q) * v[3] + E₁(t,q)
    f[2] = dA₁dx₂(t,q) * v[1] + dA₂dx₂(t,q) * v[2] + dA₃dx₂(t,q) * v[3] + E₂(t,q)
    f[3] = dA₁dx₃(t,q) * v[1] + dA₂dx₃(t,q) * v[2] + dA₃dx₃(t,q) * v[3] + E₃(t,q)
    f[4] = v[1] - q[4]
    f[5] = v[2] - q[5]
    f[6] = v[3] - q[6]
    nothing
end

# The projection force is gᵢ = (∂ϑⱼ/∂zⁱ) λʲ on the phase space z = (x, v), i.e. the transpose of
# `dϑ` above contracted with λ, evaluated at the same point — the velocities carried in `q[4:6]`.
#
# Both blocks used to be written for a one-form without the metric: the position block omitted the
# ∂gⱼⱼ/∂xⁱ vʲ terms that `dϑⱼdxᵢ` carries, and the velocity block was the identity where
# ∂ϑⱼ/∂vⁱ = gᵢᵢ δᵢⱼ. That is the same convention `dϑ` and `ω` were corrected away from, so `g` had
# been left describing a different one-form than the two functions above it.
function charged_particle_3d_iode_g(g, t, q, v, λ, params)
    g[1] = dϑ₁dx₁(t,q) * λ[1] + dϑ₂dx₁(t,q) * λ[2] + dϑ₃dx₁(t,q) * λ[3]
    g[2] = dϑ₁dx₂(t,q) * λ[1] + dϑ₂dx₂(t,q) * λ[2] + dϑ₃dx₂(t,q) * λ[3]
    g[3] = dϑ₁dx₃(t,q) * λ[1] + dϑ₂dx₃(t,q) * λ[2] + dϑ₃dx₃(t,q) * λ[3]
    g[4] = g₁₁(t,q) * λ[1]
    g[5] = g₂₂(t,q) * λ[2]
    g[6] = g₃₃(t,q) * λ[3]
    nothing
end


@doc raw"""
The drift half of the splitting, at frozen velocity: ``\dot{x} = v``, ``\dot{v} = 0``, whose flow
is exact.

`GeometricEquations` calls the solution map of an `SODE` substep as
`q(q₁, t₁, q₀, t₀, params)`; the trailing `params` used to be missing from both maps here, so
`check_methods` rejected the problem.
"""
function charged_particle_3d_sode_fx(q₁::AbstractArray{DT}, t₁, q₀::AbstractArray{DT}, t₀, params) where {DT}
    @assert axes(q₁) == axes(q₀)

    x = q₀[1:3]
    v = q₀[4:6]

    q₁[1:3] .= x .+ (t₁ - t₀) .* v
    q₁[4:6] .= v

    nothing
end

# The vector field of the substep above, for integrators that want it rather than the exact flow.
function charged_particle_3d_sode_vx(v, t, q, params)
    v[1] = q[4]
    v[2] = q[5]
    v[3] = q[6]
    v[4] = zero(eltype(q))
    v[5] = zero(eltype(q))
    v[6] = zero(eltype(q))

    nothing
end

@doc raw"""
The velocity half of the splitting, at frozen position: ``\dot{x} = 0``,
``\dot{v} = E(x) + v \times B(x)``.

Writing ``\hat{B}_{ij} = - \varepsilon_{ijk} B_{k}``, so that ``\hat{B} v = v \times B``, the
implicit midpoint rule applied to the frozen-position equation gives the update

```math
\bigg( \mathbb{1} - \frac{h}{2} \hat{B} \bigg) v_{1}
    = \bigg( \mathbb{1} + \frac{h}{2} \hat{B} \bigg) v_{0} + h \, E(x) ,
```

which is the Boris push. For ``E = 0`` it reduces to the Cayley transform of ``\hat{B}``, an exact
rotation, and then preserves ``\vert v \vert`` exactly; with ``E \neq 0`` the electric field does
work on the particle and ``\vert v \vert`` is not conserved.

Together with `charged_particle_3d_sode_fx`, which advances ``\dot{x} = v`` at frozen
velocity, this covers the whole of `charged_particle_3d_v`. The electric field used to be missing
from both halves, so the splitting integrated the ``E = 0`` system regardless of the equilibrium.

As with `charged_particle_3d_v`, this is the cartesian Lorentz force; see the model audit for the
note on curvilinear coordinates.
"""
function charged_particle_3d_sode_fv(q₁::AbstractArray{DT}, t₁, q₀::AbstractArray{DT}, t₀, params) where {DT}
    @assert axes(q₁) == axes(q₀)

    x = @view q₀[1:3]
    v = @view q₀[4:6]

    local lB₁ = B₁(t₀, x)
    local lB₂ = B₂(t₀, x)
    local lB₃ = B₃(t₀, x)

    # B̂ with B̂ v = v × B
    local B̂ = zeros(DT, 3, 3)

    B̂[1,2] = + lB₃
    B̂[1,3] = - lB₂
    B̂[2,1] = - lB₃
    B̂[2,3] = + lB₁
    B̂[3,1] = + lB₂
    B̂[3,2] = - lB₁

    local lE = DT[E₁(t₀, x), E₂(t₀, x), E₃(t₀, x)]

    local h = t₁ - t₀
    v₁ = (I - h / 2 .* B̂) \ ((I + h / 2 .* B̂) * v .+ h .* lE)

    q₁[1:3] .= x
    q₁[4:6] .= v₁

    nothing
end

# The vector field of the substep above, for integrators that want it rather than the exact flow.
function charged_particle_3d_sode_vv(v, t, q, params)
    v[1] = zero(eltype(q))
    v[2] = zero(eltype(q))
    v[3] = zero(eltype(q))
    v[4] = v₄(t, q, v)
    v[5] = v₅(t, q, v)
    v[6] = v₆(t, q, v)

    nothing
end


function charged_particle_3d_ode(qᵢ=qᵢ; tspan=tspan, tstep=Δt, periodic=true)
    if periodic
        ODEProblem(charged_particle_3d_v, tspan, tstep, qᵢ; invariants=(h=hamiltonian,), periodicity=charged_particle_3d_periodicity(qᵢ))
    else
        ODEProblem(charged_particle_3d_v, tspan, tstep, qᵢ; invariants=(h=hamiltonian,))
    end
end

# Both the vector fields and their exact solution maps are supplied, index by index, which is the
# form `SODE` is happiest with: a composition method then uses the exact flows, while anything that
# wants to sub-integrate a substep has the vector field available.
#
# Passing `nothing` for the vector fields, as this used to, does not work: the typed constructor
# requires `v::Tuple`, and the fallback `SODEProblem(v, args...) = SODEProblem(v, nothing, args...)`
# then recurses on itself forever, so the constructor overflowed the stack rather than building.
function charged_particle_3d_sode(qᵢ=qᵢ; tspan=tspan, tstep=Δt, periodic=true)
    if periodic
        SODEProblem((charged_particle_3d_sode_vv, charged_particle_3d_sode_vx),
                    (charged_particle_3d_sode_fv, charged_particle_3d_sode_fx),
                    tspan, tstep, qᵢ, periodicity=charged_particle_3d_periodicity(qᵢ))
    else
        SODEProblem((charged_particle_3d_sode_vv, charged_particle_3d_sode_vx),
                    (charged_particle_3d_sode_fv, charged_particle_3d_sode_fx),
                    tspan, tstep, qᵢ)
    end
end

function charged_particle_3d_iode(qᵢ=qᵢ; tspan=tspan, tstep=Δt)
    IODEProblem(
        charged_particle_3d_iode_ϑ,
        charged_particle_3d_iode_f,
        charged_particle_3d_iode_g,
        tspan, tstep, qᵢ, charged_particle_3d_pᵢ(tspan[begin], qᵢ);
        invariants = (h = hamiltonian,),
        v̄ = charged_particle_3d_v)
end

function charged_particle_3d_lode(qᵢ=qᵢ; tspan=tspan, tstep=Δt)
    LODEProblem(
        charged_particle_3d_iode_ϑ,
        charged_particle_3d_iode_f,
        charged_particle_3d_iode_g,
        ω, lagrangian,
        tspan, tstep, qᵢ, charged_particle_3d_pᵢ(tspan[begin], qᵢ);
        invariants = (h = hamiltonian,),
        v̄ = charged_particle_3d_v)
end


# `Union{TimeSeries, ScalarDataSeries}`: since GeometricSolutions 0.6 a solution's `sol.t` is a
# `ScalarDataSeries`, not a `TimeSeries`, so annotating only the latter made the `GeometricSolution`
# methods below fail to dispatch.
compute_energy(t::Union{TimeSeries, ScalarDataSeries}, q::DataSeries, params) = compute_invariant(t, q, params, hamiltonian)
compute_energy(sol::GeometricSolution) = compute_energy(sol.t, sol.q, GeometricEquations.parameters(sol.problem))

compute_energy_error(t::Union{TimeSeries, ScalarDataSeries}, q::DataSeries, params) = compute_invariant_error(t, q, params, hamiltonian)
compute_energy_error(sol::GeometricSolution) = compute_energy_error(sol.t, sol.q, GeometricEquations.parameters(sol.problem))
