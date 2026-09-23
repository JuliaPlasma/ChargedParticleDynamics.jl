# The problem constructors of one noncanonical equilibrium module. They forward to the
# family-level constructors in `Noncanonical`, with the module's own `qᵢ` and
# `default_parameters()` — which carries its field — as the defaults. Included into each
# noncanonical equilibrium module.

import ..Noncanonical
using ..Noncanonical: hamiltonian

for problem in (:odeproblem, :sodeproblem, :iodeproblem, :lodeproblem)
    @eval begin
        function $problem(q₀ = qᵢ; parameters = default_parameters(), kwargs...)
            Noncanonical.$problem(q₀; parameters = parameters, kwargs...)
        end
        $problem(ics::NamedTuple; kwargs...) = Noncanonical.$problem(ics; kwargs...)
    end
end
