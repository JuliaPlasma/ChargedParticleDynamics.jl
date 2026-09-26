using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "Plots" include("Plots.jl")
end
if "slow" in GROUPS
    @safetestset "Structure" include("integration/structure.jl")
    @safetestset "Charged particle 3D" include("ChargedParticle3d.jl")
    @safetestset "Guiding centre 3D" include("GuidingCenter3d.jl")
    @safetestset "Guiding centre 4D" include("GuidingCenter4d.jl")
    @safetestset "Gyrokinetics 4D" include("GyroKinetics4d.jl")
    @safetestset "Pauli particle 3D" include("PauliParticle3d.jl")
    @safetestset "Poincaré invariants" include("integration/poincare_invariants.jl")
    @safetestset "Model agreement" include("integration/model_agreement.jl")
end
