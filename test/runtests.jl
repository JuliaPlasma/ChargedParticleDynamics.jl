
using SafeTestsets
using GeometricIntegrators.Config

# solver settings
set_config(:nls_stol_break, Inf)


include("charged_particle_3d_tests.jl")
include("guiding_center_4d_tests.jl")
