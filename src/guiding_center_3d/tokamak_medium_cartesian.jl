"""
Analytic axisymmetric medium-size tokamak equilibrium in cartesian coordinates.
"""
module TokamakMediumCartesian

import ElectromagneticFields.AxisymmetricTokamakCartesian

export initial_conditions_barely_passing, initial_conditions_barely_trapped,
    initial_conditions_deeply_passing, initial_conditions_deeply_trapped

export hamiltonian

AxisymmetricTokamakCartesian.@code(2.0, 5.0, 2.0) # inject magnetic field code

# `hodeproblem` and `hodeproblem_compact` hold this thousand-step example on the module's `:g23` pair,
# at 2.6E-6 and 2.2E-9 in `max|gᵏ|` for `barely_passing`. `hodeproblem_canonical` does not: it loses
# `barely_trapped` and `deeply_passing` outright over it, at 4.1E+07 and 8.0E+10 in relative energy.
#
# It holds a thousand steps at a tenth of the step, and therefore a tenth of the span — `Δt = 0.01` over
# `t ∈ [0, 10]` gives 2.8E-13, 2.6E-13, 9.6E-08 and 6.7E-15 — so the exception it carries in
# `test/guiding_center_3d_tests.jl` is a step and a span together, not a step alone. `(g¹, g²)` is not
# an alternative: it meets a `NaN` or a `SingularException` at every step size from 0.1 down to 0.01.
# The canonicalised form is the fragile one here, as it is on `SolovevSymmetricField`.
const DEFAULT_TIMESTEP = 0.1
const DEFAULT_TIMESPAN = (0.0, 1E2)

# `y = 0`, matching this equilibrium's `GuidingCenter4d` counterpart and every other module in the
# package. This read `y = 0.1` while the module defaulted to `:g31`, and had to: in a cartesian chart
# `B_x` vanishes on `y = z = 0`, so `b₁ = 0` exactly here and the `(g³, g¹)` pair — which divides by it
# — is singular. The offset bought that one pair, at the cost of being the only initial condition in the
# package that did not agree with its 4D and Pauli siblings.
#
# Moving to `y = 0` and defaulting to `(g², g³)` instead is the better trade: `b₂ = 0.9923` is the
# largest component of `b` here, so it is also the best conditioned pair — `λₒ = -3.85` against
# `:g31`'s `-0.154` at `y = 0.1` — and `hodeproblem` conserves the constraints about three times better
# on it. What it costs is that this equilibrium is no longer one where all three pairs are regular at
# once; `SolovevIterXpoint` is now the sharpest case for that comparison. See `default_constraints`
# below and the conditioning table in `docs/src/findings.md`.
initial_conditions_barely_passing() = merge(initial_conditions(0, [2.5, 0.0, 0.0, 3.425E-1]), (params=(μ=1E-2,),))
initial_conditions_barely_trapped() = merge(initial_conditions(0, [2.5, 0.0, 0.0, 3.375E-1]), (params=(μ=1E-2,),))
initial_conditions_deeply_passing() = merge(initial_conditions(0, [2.5, 0.0, 0.0, 5E-1]), (params=(μ=1E-2,),))
initial_conditions_deeply_trapped() = merge(initial_conditions(0, [2.5, 0.0, 0.0, 1E-1]), (params=(μ=1E-2,),))

μ_loop() = 1E-3
μ_surface() = 1E-3

function f_loop(s)
    X0 = 1.75
    Y0 = 0.0
    Z0 = 0.0
    u0 = 0.5
    rx = 0.1
    ry = 0.1

    Xs = X0 + rx * cos(2π * s)
    Zs = Z0 + ry * sin(2π * s)

    qs = [Xs, Y0, Zs, u0]

    return qs
end

function f_surface(s, t)
    X0 = 1.75
    Y0 = 0.0
    Z0 = 0.0
    u0 = 0.5
    r0 = 0.1

    Xt = X0 + r0 * (s - 0.5)
    Zt = Z0 + r0 * (t - 0.5)

    qt = [Xt, Y0, Zt, u0]

    return qt
end

export default_parameters, default_constraints

"The magnetic moment μ this equilibrium is set up for."
default_parameters(::Type{T}=Float64) where {T} = (μ = T(1E-2),)

"""
The constraint pair [`hodeproblem`](@ref) and its siblings use here by default.

`(g², g³)`, which divides by `b₂ = 0.9923` — the largest component of `b` at this initial condition, so
the best conditioned of the three pairs. `(g³, g¹)` divides by `b₁`, which vanishes identically on
`y = z = 0` in a cartesian chart, so it is not merely badly conditioned here but singular.

This is one of two equilibria whose default is chosen on conditioning rather than on bare regularity;
`Dipole3d` is the other, and its reason is different — there all three pairs are regular at the initial
condition and it is the *orbit* that takes `b₁` and `b₂` through zero. See `docs/src/findings.md`.

!!! note "`hodeproblem_canonical` needs a smaller step here"
    On this pair the canonicalised form loses `barely_trapped` and `deeply_passing` at the module's
    `Δt = 0.1`, and holds all four at `Δt = 0.01`. `(g¹, g²)` is not an alternative — it meets a `NaN`
    or a `SingularException` at every step size tried, including 0.01. The step, not the pair, is what
    that formulation needs here.
"""
default_constraints() = :g23


include("guiding_center_3d_equations.jl")
include("guiding_center_3d_canonical.jl")
include("guiding_center_3d_diagnostics.jl")

end
