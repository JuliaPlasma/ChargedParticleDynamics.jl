
import GeometricEquations
using GeometricSolutions
using GeometricSolutions: compute_error_drift, compute_momentum_error, compute_one_form

export compute_energy, compute_toroidal_momentum,
       compute_energy_error, compute_toroidal_momentum_error,
       compute_momentum_error, compute_one_form, compute_error_drift


compute_energy(t::TimeSeries, q::DataSeries, params) = compute_invariant(t, q, params, hamiltonian)
compute_energy(sol::GeometricSolution, params = GeometricEquations.parameters(sol.problem)) = compute_energy(sol.t, sol.q, params)

compute_energy_error(t::TimeSeries, q::DataSeries, params) = compute_invariant_error(t, q, params, hamiltonian)
compute_energy_error(sol::GeometricSolution, params = GeometricEquations.parameters(sol.problem)) = compute_energy_error(sol.t, sol.q, params)

compute_toroidal_momentum(t::TimeSeries, q::DataSeries, params) = compute_invariant(t, q, params, (t, q, params) -> toroidal_momentum(t, q))
compute_toroidal_momentum(sol::GeometricSolution, params = GeometricEquations.parameters(sol.problem)) = compute_toroidal_momentum(sol.t, sol.q, params)

compute_toroidal_momentum_error(t::TimeSeries, q::DataSeries, params) = compute_invariant_error(t, q, params, (t, q, params) -> toroidal_momentum(t, q))
compute_toroidal_momentum_error(sol::GeometricSolution, params = GeometricEquations.parameters(sol.problem)) = compute_toroidal_momentum_error(sol.t, sol.q, params)
