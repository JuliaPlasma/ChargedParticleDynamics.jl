"""
Analytic ITER-like Solov'ev equilibrium with X-point.
"""
module SolovevIterXpoint

using ElectromagneticFields: FieldFunctions, SolovevXpointEquilibriumITER, from_cartesian
using ...ChargedParticleDynamics: check_chart

export initial_conditions_barely_passing, initial_conditions_barely_trapped,
       initial_conditions_deeply_passing, initial_conditions_deeply_trapped,
       initial_conditions_trapped

export hamiltonian, toroidal_momentum

const FIELD = FieldFunctions(SolovevXpointEquilibriumITER())

const DEFAULT_TIMESTEP = 1.0
const DEFAULT_TIMESPAN = (0.0, 1000.0)

const xᵢ = [7.0-1.4, 0.0, 0.0]
const qᵢ = [from_cartesian(FIELD, 0, xᵢ)..., 2.8166280889939737]

export default_parameters

"The field, and the magnetic moment μ this equilibrium is set up for."
function default_parameters(::Type{T} = Float64) where {T}
    (field = FIELD, μ = T(4.607782183567846))
end

include("guiding_center_4d_presets.jl")

# The canonical toroidal momentum is the covariant φ-component of the one-form, ϑ₃. It was
# previously multiplied by R, which destroys the conservation: on the small tokamak the
# relative variation over 10³ time units is 2e-13 for ϑ₃ and 3e-3 for R ϑ₃.
function toroidal_momentum(t, q, params)
    check_chart(params.field, FIELD)
    ϑ₃(t, q, params)
end

include("guiding_center_4d_diagnostics.jl")

function initial_conditions_barely_passing()
    (q = [from_cartesian(FIELD, 0, [2.5, 0.0, 0.0])..., 3.425E-1],
        params = (field = FIELD, μ = 1E-2))
end
function initial_conditions_barely_trapped()
    (q = [from_cartesian(FIELD, 0, [2.5, 0.0, 0.0])..., 3.375E-1],
        params = (field = FIELD, μ = 1E-2))
end
function initial_conditions_deeply_passing()
    (q = [from_cartesian(FIELD, 0, [2.5, 0.0, 0.0])..., 5E-1],
        params = (field = FIELD, μ = 1E-2))
end
function initial_conditions_deeply_trapped()
    (q = [from_cartesian(FIELD, 0, [2.5, 0.0, 0.0])..., 1E-1],
        params = (field = FIELD, μ = 1E-2))
end
function initial_conditions_trapped()
    (q = [from_cartesian(FIELD, 0, [7.0, 0.0, 0.0])..., -2E-3],
        params = (field = FIELD, μ = 1.88E-7))
end

end
