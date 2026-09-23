# The problem constructors of one equilibrium module. They forward to the family-level
# constructors in `GyroKinetics4d`, with the module's own `qᵢ`, `default_parameters()` — which
# carries its field — `DEFAULT_TIMESPAN` and `DEFAULT_TIMESTEP` as the defaults. Included into each
# equilibrium module.

import ..GyroKinetics4d
using ..GyroKinetics4d: hamiltonian, ϑ, ω, ωabs, β, γ, v
using ..GyroKinetics4d: transform_q_to_q̃!, transform_q_to_q̃, transform_q̃_to_q!,
                        transform_q̃_to_q

export hamiltonian, ϑ, ω, ωabs, β, γ, v
export transform_q_to_q̃!, transform_q_to_q̃, transform_q̃_to_q!, transform_q̃_to_q
export odeproblem, sodeproblem

for problem in (:odeproblem, :sodeproblem)
    @eval begin
        function $problem(
                q₀ = qᵢ; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP,
                parameters = default_parameters(), kwargs...)
            GyroKinetics4d.$problem(q₀; timespan = timespan, timestep = timestep,
                parameters = parameters, kwargs...)
        end

        # `initial_conditions_*` returns `(q = …, params = …)`; take it directly.
        $problem(ics::NamedTuple; kwargs...) = $problem(ics.q; parameters = ics.params, kwargs...)
    end
end
