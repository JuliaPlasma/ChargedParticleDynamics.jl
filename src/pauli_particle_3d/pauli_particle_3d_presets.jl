# The initial conditions and problem constructors of one equilibrium module. They forward to the
# family-level functions in `PauliParticle3d`, with the module's own field, `qᵢ`, `vᵢ`,
# `default_parameters()`, `DEFAULT_TIMESPAN` and `DEFAULT_TIMESTEP` as the defaults. Included into
# each equilibrium module.

import ..PauliParticle3d
using ..PauliParticle3d: hamiltonian

function initial_conditions(x₀, v₀::AbstractVector)
    PauliParticle3d.initial_conditions(FIELD, x₀, v₀)
end
initial_conditions(x₀, u₀::Real, μ) = PauliParticle3d.initial_conditions(FIELD, x₀, u₀, μ)

for problem in (:podeproblem, :hodeproblem, :iodeproblem)
    @eval begin
        function $problem(
                q₀::AbstractVector, v₀::AbstractVector; timespan = DEFAULT_TIMESPAN,
                timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
            PauliParticle3d.$problem(q₀, v₀; timespan = timespan, timestep = timestep,
                parameters = parameters)
        end

        # `μ` on its own, for callers that have the moment rather than a parameter tuple
        function $problem(q₀::AbstractVector, v₀::AbstractVector, μ::Real; kwargs...)
            $problem(q₀, v₀; parameters = (field = FIELD, μ = μ), kwargs...)
        end

        # the named tuple `initial_conditions` returns, whose `v` is already the parallel velocity
        function $problem(ics::NamedTuple; kwargs...)
            $problem(ics.q, ics.v; parameters = ics.params, kwargs...)
        end

        # no arguments: split the module's own `(qᵢ, vᵢ)`. This must go through
        # `initial_conditions` rather than defaulting `v₀ = vᵢ` — `vᵢ` is the *full* velocity, and
        # handing it to the constructor as if it were the parallel one puts the particle on a
        # trajectory the solver cannot follow.
        $problem(; kwargs...) = $problem(initial_conditions(qᵢ, vᵢ); kwargs...)
    end
end
