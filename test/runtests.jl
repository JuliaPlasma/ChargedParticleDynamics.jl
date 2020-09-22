
# solver settings
using GeometricIntegrators.Config

set_config(:nls_atol_break, 1E3)
set_config(:nls_rtol_break, 1E3)
set_config(:nls_stol_break, 1E3)

set_config(:nls_atol, 1E-14)
set_config(:nls_rtol, 1E-14)
set_config(:nls_stol, 1E-14)


include("charged_particle_3d_tests.jl")
include("guiding_center_4d_tests.jl")
include("gyro_kinetics_4d_tests.jl")
include("pauli_particle_3d_tests.jl")
