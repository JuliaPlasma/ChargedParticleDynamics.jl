@doc raw"""
Analytic ITER-like Solov'ev equilibrium, in cylindrical coordinates ``(R, Z, \varphi)``.

The gyrokinetic counterpart of `GuidingCenter4d.SolovevIter`, sharing its initial
conditions. Its independent variable is the rescaled time ``s``, with
``dt = B^{\star}_{\parallel} \, ds``; here ``B^{\star}_{\parallel} \approx 2 \times 10^{2}`` at the
deeply passing initial condition, so the defaults below are the 4D guiding centre's divided by that.
"""
module GuidingCenter4dSolovevIter

using ElectromagneticFields: FieldFunctions, SolovevEquilibriumITER, from_cartesian
using ...ChargedParticleDynamics: check_chart
using ..GyroKinetics4d: fieldpoint, ϑ₃

export initial_conditions_barely_passing, initial_conditions_barely_trapped,
       initial_conditions_deeply_passing, initial_conditions_deeply_trapped,
       initial_conditions_trapped

export hamiltonian, toroidal_momentum

const FIELD = FieldFunctions(SolovevEquilibriumITER())

include("gc_presets.jl")

# Rescaled time: the 4D guiding centre uses Δt = 1.0 over (0, 10³) and B*∥ ≈ 204 here, which
# puts the equivalent step at 5E-3. The splitting does not hold together over the full span at
# that step, so this matches its X-point sibling instead — the same Solov'ev field with the
# same B*∥, whose 1E-3 is the value that equilibrium was already shipped with.
const DEFAULT_TIMESTEP = 1E-3
const DEFAULT_TIMESPAN = (0.0, 1.0)

const x₀ = from_cartesian(FIELD, 0, [2.5, 0.0, 0.0])
const μ₀ = 1E-2

export default_parameters

"""
The field, and the magnetic moment `μ` of the deeply passing initial condition.
"""
default_parameters(::Type{T} = Float64) where {T} = (field = FIELD, μ = T(μ₀))

const qᵢ = [x₀..., 5.0E-1]

function toroidal_momentum(t, q, params)
    check_chart(params.field, FIELD)
    q = fieldpoint(params.field, t, q)
    ϑ₃(t, q)
end

initial_conditions_barely_passing() = (
    q = [x₀..., 3.425E-1], params = (field = FIELD, μ = μ₀))
initial_conditions_barely_trapped() = (
    q = [x₀..., 3.375E-1], params = (field = FIELD, μ = μ₀))
initial_conditions_deeply_passing() = (
    q = [x₀..., 5.0E-1], params = (field = FIELD, μ = μ₀))
initial_conditions_deeply_trapped() = (
    q = [x₀..., 1.0E-1], params = (field = FIELD, μ = μ₀))

initial_conditions_trapped() = (
    q = [from_cartesian(FIELD, 0, [7.0, 0.0, 0.0])..., -2E-3], params = (
        field = FIELD, μ = 1.88E-7))

include("../guiding_center_4d/guiding_center_4d_diagnostics.jl")

end
