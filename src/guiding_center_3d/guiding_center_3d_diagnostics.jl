
import GeometricEquations
using GeometricSolutions
using GeometricSolutions: compute_error_drift, compute_momentum_error, compute_one_form

export compute_energy, compute_energy_error,
       compute_toroidal_momentum, compute_toroidal_momentum_error,
       compute_constraints,
       compute_momentum_error, compute_one_form, compute_error_drift


# `TimeSeries` accepts the time axis of a `TimeSeries` object, `ScalarDataSeries` the one a
# `GeometricSolution` stores in `sol.t`; see the note in the 4D diagnostics for why both are needed.
const SolutionTimes = Union{TimeSeries, ScalarDataSeries}


# The 3D model is a `HODEProblem`, so unlike the 4D one its invariants are functions of `(q, p)`
# and the diagnostics take both series.

compute_energy(t::SolutionTimes, q::DataSeries, p::DataSeries, params) =
    compute_invariant(t, q, p, params, hamiltonian)
compute_energy(sol::GeometricSolution, params=GeometricEquations.parameters(sol.problem)) =
    compute_energy(sol.t, sol.q, sol.p, params)

"""
    compute_energy_error(sol)

!!! note "Returns a `(value, error)` pair"
    Like every `compute_*_error` here, this forwards to `GeometricSolutions.compute_invariant_error`
    and so returns **two** series, not one: the invariant itself and its relative error. Destructure
    it — `h, e = compute_energy_error(sol)` — rather than passing the result somewhere a single
    series is expected. [`compute_constraints`](@ref) is the exception: the constraints start at
    zero, so it returns the values alone.
"""
compute_energy_error(t::SolutionTimes, q::DataSeries, p::DataSeries, params) =
    compute_invariant_error(t, q, p, params, hamiltonian)
compute_energy_error(sol::GeometricSolution, params=GeometricEquations.parameters(sol.problem)) =
    compute_energy_error(sol.t, sol.q, sol.p, params)


# Defined only by the equilibria that have a symmetry direction; calling these on one that does not
# raises an `UndefVarError` for `toroidal_momentum`, as it does for the 4D model.
compute_toroidal_momentum(t::SolutionTimes, q::DataSeries, p::DataSeries, params) =
    compute_invariant(t, q, p, params, (t, q, p, params) -> toroidal_momentum(t, q, p))
compute_toroidal_momentum(sol::GeometricSolution, params=GeometricEquations.parameters(sol.problem)) =
    compute_toroidal_momentum(sol.t, sol.q, sol.p, params)

"""
    compute_toroidal_momentum_error(sol)

Returns a `(value, error)` pair; see [`compute_energy_error`](@ref).
"""
compute_toroidal_momentum_error(t::SolutionTimes, q::DataSeries, p::DataSeries, params) =
    compute_invariant_error(t, q, p, params, (t, q, p, params) -> toroidal_momentum(t, q, p))
compute_toroidal_momentum_error(sol::GeometricSolution, params=GeometricEquations.parameters(sol.problem)) =
    compute_toroidal_momentum_error(sol.t, sol.q, sol.p, params)


@doc raw"""
    compute_constraints(sol)

All three constraints ``g^{1}``, ``g^{2}``, ``g^{3}`` that confine the solution to the manifold on
which the momentum equals the guiding centre one-form, returned as a named tuple of three
`DataSeries`.

Any two of the three make up the constraint pair a problem was built with — see
[`constraint_pair`](@ref) — but all three vanish along the flow whichever pair that was, so
reporting all three is independent of the choice and shows a drift that is invisible in the retained
pair alone.

They vanish identically along the continuous flow, so their magnitude measures how far a numerical
solution has drifted off that manifold — the quantity the approximately symplectic methods of Li,
Zhang and Liu are designed to keep bounded, and the one that decides whether `hodeproblem_canonical`
and `hodeproblem_compact` agree with `hodeproblem`. Since they start at zero, the absolute value is
the meaningful one; there is no relative-error variant.
"""
function compute_constraints(t::SolutionTimes, q::DataSeries, p::DataSeries)
    (g₁=DataSeries([g₁(t[i], q[i], p[i]) for i in eachindex(t)]),
     g₂=DataSeries([g₂(t[i], q[i], p[i]) for i in eachindex(t)]),
     g₃=DataSeries([g₃(t[i], q[i], p[i]) for i in eachindex(t)]))
end

compute_constraints(sol::GeometricSolution) = compute_constraints(sol.t, sol.q, sol.p)
