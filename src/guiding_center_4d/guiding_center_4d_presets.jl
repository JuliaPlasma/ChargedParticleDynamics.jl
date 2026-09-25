# The problem constructors of one equilibrium module. They forward to the family-level
# constructors in `GuidingCenter4d`, with the module's own `qᵢ`, `default_parameters()` — which
# carries its field — `DEFAULT_TIMESPAN` and `DEFAULT_TIMESTEP` as the defaults. Included into each
# equilibrium module, after its `qᵢ` if it has one.

import ..GuidingCenter4d
using ..GuidingCenter4d: hamiltonian, u, ω, ϑ, ϑ₁, ϑ₂, ϑ₃, ϑ₄, dϑ, β₁, β₂, β₃, dH

export hamiltonian, u, ω, ϑ, ϑ₁, ϑ₂, ϑ₃, ϑ₄, dϑ, β₁, β₂, β₃, dH
export odeproblem, iodeproblem, iodeproblem_λ,
       lodeproblem,
       iodeproblem_dg, lodeproblem_formal_lagrangian

for problem in (:odeproblem, :iodeproblem, :iodeproblem_λ, :lodeproblem,
    :iodeproblem_dg, :lodeproblem_formal_lagrangian)
    @eval begin
        function $problem(
                q₀; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP,
                parameters = default_parameters(), kwargs...)
            GuidingCenter4d.$problem(q₀; timespan = timespan, timestep = timestep,
                parameters = parameters, kwargs...)
        end

        # Every `initial_conditions_*` returns `(q = …, params = …)`, so each constructor takes that
        # named tuple directly. The parameters travel with the initial condition because `μ`
        # differs between them.
        $problem(ics::NamedTuple; kwargs...) = $problem(ics.q; parameters = ics.params, kwargs...)
    end

    # No arguments: a copy of the module's own `qᵢ`, which the problem would otherwise share. The
    # Poincaré-invariant fixtures have a loop or a surface instead of a point, so they have no `qᵢ`
    # and no zero-argument method.
    if isdefined(@__MODULE__, :qᵢ)
        @eval $problem(; kwargs...) = $problem(copy(qᵢ); kwargs...)
    end
end
