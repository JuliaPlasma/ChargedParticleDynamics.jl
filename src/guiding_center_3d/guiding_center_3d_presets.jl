# The initial conditions and problem constructors of one equilibrium module. They forward to the
# family-level functions in `GuidingCenter3d`, with the module's own `qᵢ`, `default_parameters()` —
# which carries its field — `default_constraints()`, `DEFAULT_TIMESPAN` and `DEFAULT_TIMESTEP` as
# the defaults. Included into each equilibrium module, after its `qᵢ`.

import ..GuidingCenter3d
using ..GuidingCenter3d: hamiltonian, hamiltonian_canonical, g₁, g₂, g₃

export hamiltonian, hamiltonian_canonical
export hodeproblem, hodeproblem_canonical, hodeproblem_compact

# The momentum is the one-form at `Qᵢ = (x, u)`, which reads the field and not `μ`.
function initial_conditions(tᵢ, Qᵢ)
    GuidingCenter3d.initial_conditions(tᵢ, Qᵢ, default_parameters())
end

for problem in (:hodeproblem, :hodeproblem_canonical, :hodeproblem_compact)
    # The compact form is independent of the pair by default; the other two retain this
    # equilibrium's own. See `default_constraints`.
    constraints = problem == :hodeproblem_compact ? QuoteNode(:parallel) :
                  :(default_constraints())

    @eval begin
        function $problem(
                q₀::AbstractVector, p₀::AbstractVector; timespan = DEFAULT_TIMESPAN,
                timestep = DEFAULT_TIMESTEP, parameters = default_parameters(),
                constraints = $constraints, kwargs...)
            GuidingCenter3d.$problem(q₀, p₀; timespan = timespan, timestep = timestep,
                parameters = parameters, constraints = constraints, kwargs...)
        end

        function $problem(x₀::AbstractVector; timespan = DEFAULT_TIMESPAN,
                timestep = DEFAULT_TIMESTEP, parameters = default_parameters(),
                constraints = $constraints, kwargs...)
            GuidingCenter3d.$problem(x₀; timespan = timespan, timestep = timestep,
                parameters = parameters, constraints = constraints, kwargs...)
        end

        # the named tuple every `initial_conditions_*` returns, which carries its own `μ`
        function $problem(ics::NamedTuple; kwargs...)
            $problem(ics.q, ics.p; parameters = ics.params, kwargs...)
        end
    end

    # No arguments: the module's own `qᵢ`. The Poincaré-invariant fixtures have a loop or a
    # surface instead of a point, so they have no `qᵢ` and no zero-argument method.
    if isdefined(@__MODULE__, :qᵢ)
        @eval $problem(; kwargs...) = $problem(qᵢ; kwargs...)
    end
end
