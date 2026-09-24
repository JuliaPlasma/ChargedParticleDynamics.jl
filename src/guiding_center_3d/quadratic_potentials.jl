"""
Analytic toy problem with quadratic potentials.
"""
module QuadraticPotentials3d

using ElectromagneticFields: FieldFunctions, QuadraticPotentialsField

export initial_conditions_quadratic

export hamiltonian

const FIELD = FieldFunctions(QuadraticPotentialsField())

# A thousand steps, over which all three formulations sit between 3.1E-9 and 7.4E-9 in relative
# energy. A span of 2.5E4 at this step would be fifty thousand steps — a study rather than an
# example — and `hodeproblem_canonical` diverges to `Inf` over it.
const DEFAULT_TIMESTEP = 0.5
const DEFAULT_TIMESPAN = (0.0, 5E2)

# const DEFAULT_TIMESTEP = 0.35
# const DEFAULT_TIMESPAN = (0.0, 12250.)

const qᵢ = [0.3, 0.2, -1.4, 0.3]

function initial_conditions_quadratic()
    merge(initial_conditions(0.0, [0.3, 0.2, -1.4, 0.3]), (params = (
        field = FIELD, μ = 2.5E-3),))
end

export default_parameters, default_constraints

"The field, and the magnetic moment μ this equilibrium is set up for."
default_parameters(::Type{T} = Float64) where {T} = (field = FIELD, μ = T(2.5E-3))

"""
The constraint pair [`hodeproblem`](@ref ChargedParticleDynamics.GuidingCenter3d.hodeproblem) and its siblings use here by default.

All three pairs are regular at the initial condition, and this is the worst conditioned of them:
`b₁ = 2E-3` gives `λₒ = 0.2`, against 100 for `(g¹, g²)` and 0.3 for `(g², g³)`. It is kept because it
is the pair the package has always used here and because the numbers in `docs/src/findings.md` were
measured with it, not because it is the best available — see the conditioning table there, and the open
item in `TODO.md` about choosing the defaults on conditioning along the orbit.
"""
default_constraints() = :g31

include("guiding_center_3d_presets.jl")
include("guiding_center_3d_diagnostics.jl")

end
