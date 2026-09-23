"""
Analytic, quadratic Solov'ev equilibrium.
"""
module SolovevSymmetricField

using ElectromagneticFields: FieldFunctions, SolovevSymmetricEquilibrium
using ...ChargedParticleDynamics: check_chart

export initial_conditions_barely_passing, initial_conditions_barely_trapped,
       initial_conditions_deeply_passing, initial_conditions_deeply_trapped

export hamiltonian, toroidal_momentum

const FIELD = FieldFunctions(SolovevSymmetricEquilibrium(2.0, 5.0, 1.0, 1.0))

const DEFAULT_TIMESTEP = 1E0
const DEFAULT_TIMESPAN = (0.0, 1E3)

function initial_conditions_barely_passing()
    (q = [2.5, 0.0, 0.0, 3.425E-1], params = (field = FIELD, μ = 1E-2))
end # Δt=2.5, nt=50
function initial_conditions_barely_trapped()
    (q = [2.5, 0.0, 0.0, 3.375E-1], params = (field = FIELD, μ = 1E-2))
end # Δt=3.0, nt=100
function initial_conditions_deeply_passing()
    (q = [2.5, 0.0, 0.0, 5E-1], params = (field = FIELD, μ = 1E-2))
end     # Δt=2.5, nt=25
function initial_conditions_deeply_trapped()
    (q = [2.5, 0.0, 0.0, 1E-1], params = (field = FIELD, μ = 1E-2))
end     # Δt=5.0, nt=50

export default_parameters

"The field, and the magnetic moment μ this equilibrium is set up for."
default_parameters(::Type{T} = Float64) where {T} = (field = FIELD, μ = T(1E-2))

include("guiding_center_4d_presets.jl")

# The canonical toroidal momentum is the covariant φ-component of the one-form, ϑ₃, and not
# R ϑ₃, which destroys the conservation: on the small tokamak the
# relative variation over 10³ time units is 2e-13 for ϑ₃ and 3e-3 for R ϑ₃.
function toroidal_momentum(t, q, params)
    check_chart(params.field, FIELD)
    ϑ₃(t, q, params)
end

include("guiding_center_4d_diagnostics.jl")

end
