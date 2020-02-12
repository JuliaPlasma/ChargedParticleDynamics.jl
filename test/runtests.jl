
using GeometricIntegrators
using Test

# solver settings
# set_config(:nls_atol_break, 1E-3)
# set_config(:nls_rtol_break, 1E-3)
set_config(:nls_stol_break, 1E3)


include("charged_particle_3d_tests.jl")
include("guiding_center_4d_tests.jl")
