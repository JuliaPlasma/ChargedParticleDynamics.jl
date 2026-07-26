
import GeometricEquations
using GeometricSolutions
using GeometricSolutions: compute_error_drift, compute_momentum_error, compute_one_form

export compute_energy, compute_toroidal_momentum,
       compute_energy_error, compute_toroidal_momentum_error,
       compute_momentum_error, compute_one_form, compute_error_drift


# `TimeSeries` accepts the time axis of a `TimeSeries` object, `ScalarDataSeries` the one a
# `GeometricSolution` stores in `sol.t` — since GeometricSolutions 0.6 restructured the solution to
# hold a vector of states, those are two different types, and annotating only the former made the
# `GeometricSolution` methods below fail to dispatch. This is the same union
# `GeometricSolutions.compute_invariant` and the GeometricProblems diagnostics accept.
const SolutionTimes = Union{TimeSeries, ScalarDataSeries}

"""
    compute_energy(sol)

The Hamiltonian along the solution, as a `DataSeries`.
"""
compute_energy(t::SolutionTimes, q::DataSeries, params) = compute_invariant(t, q, params, hamiltonian)
compute_energy(sol::GeometricSolution, params = GeometricEquations.parameters(sol.problem)) = compute_energy(sol.t, sol.q, params)

"""
    compute_energy_error(sol)

!!! note "Returns a `(value, error)` pair"
    Like every `compute_*_error` here, this forwards to `GeometricSolutions.compute_invariant_error`
    and so returns **two** series, not one: the invariant itself and its relative error. Destructure
    it — `h, e = compute_energy_error(sol)` — rather than passing the result somewhere a single
    series is expected.
"""
compute_energy_error(t::SolutionTimes, q::DataSeries, params) = compute_invariant_error(t, q, params, hamiltonian)
compute_energy_error(sol::GeometricSolution, params = GeometricEquations.parameters(sol.problem)) = compute_energy_error(sol.t, sol.q, params)

"""
    compute_toroidal_momentum(sol)

The canonical toroidal momentum `ϑ₃` along the solution, as a `DataSeries`. It is conserved in an
axisymmetric equilibrium.
"""
compute_toroidal_momentum(t::SolutionTimes, q::DataSeries, params) = compute_invariant(t, q, params, (t, q, params) -> toroidal_momentum(t, q))
compute_toroidal_momentum(sol::GeometricSolution, params = GeometricEquations.parameters(sol.problem)) = compute_toroidal_momentum(sol.t, sol.q, params)

"""
    compute_toroidal_momentum_error(sol)

Returns a `(value, error)` pair; see [`compute_energy_error`](@ref).
"""
compute_toroidal_momentum_error(t::SolutionTimes, q::DataSeries, params) = compute_invariant_error(t, q, params, (t, q, params) -> toroidal_momentum(t, q))
compute_toroidal_momentum_error(sol::GeometricSolution, params = GeometricEquations.parameters(sol.problem)) = compute_toroidal_momentum_error(sol.t, sol.q, params)
