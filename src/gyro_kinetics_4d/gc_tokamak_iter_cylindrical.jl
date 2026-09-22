@doc raw"""
ITER-size axisymmetric tokamak equilibrium, in cylindrical coordinates ``(R, Z, \varphi)``.

The gyrokinetic counterpart of `GuidingCenter4d.TokamakIterCylindrical`, sharing its initial
conditions. The independent variable is the rescaled time ``s``, with
``dt = B^{\star}_{\parallel} \, ds`` and ``B^{\star}_{\parallel} \approx 36`` at the deeply passing
initial condition.
"""
module GuidingCenter4dTokamakIterCylindrical

using ElectromagneticFields: FieldFunctions, AxisymmetricTokamakCylindricalITER,
                             from_cartesian

export initial_conditions_barely_passing, initial_conditions_barely_trapped,
       initial_conditions_deeply_passing, initial_conditions_deeply_trapped,
       initial_conditions_trapped

export hamiltonian, toroidal_momentum

const FIELD = FieldFunctions(AxisymmetricTokamakCylindricalITER())

include("gc_presets.jl")

# Rescaled time: the 4D guiding centre uses Δt = 1.0 over (0, 10³) with B*∥ ≈ 36 here.
const DEFAULT_TIMESTEP = 2.5E-2
const DEFAULT_TIMESPAN = (0.0, 2.5E1)

const x₀ = from_cartesian(FIELD, 0, [2.5, 0.0, 0.0])

export default_parameters

"""
The field, and the magnetic moment `μ` of the shipped initial conditions.
"""
default_parameters(::Type{T} = Float64) where {T} = (field = FIELD, μ = T(1E-2))

const qᵢ = [x₀..., 5E-1]

function toroidal_momentum(t, q, params)
    q = fieldpoint(params.field, t, q)
    ϑ₃(t, q)
end

initial_conditions_barely_passing() = (
    q = [x₀..., 3.425E-1], params = (field = FIELD, μ = 1E-2))
initial_conditions_barely_trapped() = (
    q = [x₀..., 3.375E-1], params = (field = FIELD, μ = 1E-2))
initial_conditions_deeply_passing() = (
    q = [x₀..., 5E-1], params = (field = FIELD, μ = 1E-2))
initial_conditions_deeply_trapped() = (
    q = [x₀..., 1E-1], params = (field = FIELD, μ = 1E-2))

initial_conditions_trapped() = (
    q = [from_cartesian(FIELD, 0, [7.0, 0.0, 0.0])..., -2E-3], params = (
        field = FIELD, μ = 1.88E-7))

include("../guiding_center_4d/guiding_center_4d_diagnostics.jl")

end
