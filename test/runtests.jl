include("quiet_solver_warnings.jl")
quiet_solver_warnings!()

include("structure_tests.jl")

include("charged_particle_3d_tests.jl")
include("guiding_center_3d_tests.jl")
include("guiding_center_4d_tests.jl")
include("gyro_kinetics_4d_tests.jl")
include("pauli_particle_3d_tests.jl")
include("plots_tests.jl")

if suppressed_warning_count() > 0
    @info "Suppressed $(suppressed_warning_count()) nonlinear-solver warnings " *
          "(see test/quiet_solver_warnings.jl)"
end
