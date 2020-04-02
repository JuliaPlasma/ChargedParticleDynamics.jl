
using SafeTestsets
using GeometricIntegrators.Config

# solver settings
set_config(:nls_stol_break, Inf)


@safetestset "Charged Particle Dynamics in 3D                                                 " begin include("charged_particle_3d_tests.jl") end
@safetestset "Guiding Centre Dynamics in 4D                                                   " begin include("guiding_center_4d_tests.jl") end
