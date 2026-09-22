# The problem constructors of one canonical equilibrium module. They forward to the family-level
# constructors in `Canonical`, with the module's own `qᵢ`, `pᵢ` and `default_parameters()` — which
# carries its field — as the defaults. Included into each canonical equilibrium module.

import ..Canonical
using ..Canonical: hamiltonian, lagrangian, toroidal_momentum
using ..Canonical: compute_energy, compute_energy_error

for problem in (:podeproblem, :iodeproblem, :lodeproblem)
    @eval begin
        function $problem(q₀ = qᵢ, p₀ = pᵢ; parameters = default_parameters(), kwargs...)
            Canonical.$problem(q₀, p₀; parameters = parameters, kwargs...)
        end
        $problem(ics::NamedTuple; kwargs...) = Canonical.$problem(ics; kwargs...)
    end
end
