
# solver settings
using GeometricIntegrators.Config
set_config(:nls_stol_break, Inf)


include("charged_particle_3d_tests.jl")
include("guiding_center_4d_tests.jl")
include("gyro_kinetics_4d_tests.jl")
