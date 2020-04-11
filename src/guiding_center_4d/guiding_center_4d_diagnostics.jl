
using GeometricIntegrators.Solutions

import GeometricProblems.Diagnostics: compute_invariant, compute_one_form,
       compute_invariant_error, compute_momentum_error, compute_error_drift

export compute_energy, compute_toroidal_momentum,
       compute_energy_error, compute_toroidal_momentum_error,
       compute_momentum_error, compute_one_form, compute_error_drift


function convert_coordinates_RZphi_to_xyz(q::DataSeries{DT,2}) where {DT}
    x = SDataSeries(DT, q.nd, q.nt, q.ni)
    for k in 1:q.ni
        for j in 0:q.nt
            for i in 1:q.nd
                x[1,j,k] = q[1,j,k] * cos(q[3,j,k])
                x[2,j,k] = q[1,j,k] * sin(q[3,j,k])
                x[3,j,k] = q[2,j,k]
                x[4,j,k] = q[4,j,k]
            end
        end
    end
    return x
end


compute_energy(t::TimeSeries, q::DataSeries, params::NamedTuple) = compute_invariant(t, q, (t,q) -> hamiltonian(t, q, params))
compute_energy(sol::Solution, params::NamedTuple) = compute_energy(sol.t, sol.q, params)

compute_energy_error(t::TimeSeries, q::DataSeries, params::NamedTuple) = compute_invariant_error(t, q, (t,q) -> hamiltonian(t, q, params))
compute_energy_error(sol::Solution, params::NamedTuple) = compute_energy_error(sol.t, sol.q, params)

compute_toroidal_momentum(t::TimeSeries, q::DataSeries) = compute_invariant(t, q, toroidal_momentum)
compute_toroidal_momentum(sol::Solution) = compute_toroidal_momentum(sol.t, sol.q)

compute_toroidal_momentum_error(t::TimeSeries, q::DataSeries) = compute_invariant_error(t, q, toroidal_momentum)
compute_toroidal_momentum_error(sol::Solution) = compute_toroidal_momentum_error(sol.t, sol.q)

compute_momentum_error(t::TimeSeries, q::DataSeries, p::DataSeries) = compute_momentum_error(t, q, p, ϑ)
compute_momentum_error(sol::Solution) = compute_momentum_error(sol.t, sol.q, sol.p)
