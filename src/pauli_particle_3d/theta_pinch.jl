"""
Charged Particle in an uniform magnetic field of the form
``B(x,y,z) = B_0 e_z``.
"""
module ThetaPinchField

using ElectromagneticFields: FieldFunctions, ThetaPinchEquilibrium

export podeproblem, hamiltonian, angular_momentum

const FIELD = FieldFunctions(ThetaPinchEquilibrium())

const qᵢ = [1.0, 0.0, 0.0]
const vᵢ = [0.0, 1.0, 1.0]

const DEFAULT_TIMESTEP = 10.0
const DEFAULT_TIMESPAN = (0.0, 1E4)

include("pauli_particle_3d_presets.jl")

export default_parameters

"""
The field, and the magnetic moment μ of the default initial condition `(qᵢ, vᵢ)`, obtained by
splitting `vᵢ` into its parallel and perpendicular parts at `qᵢ`.
"""
function default_parameters(::Type{T} = Float64) where {T}
    (field = FIELD, μ = T(initial_conditions(qᵢ, vᵢ).params.μ))
end

end
