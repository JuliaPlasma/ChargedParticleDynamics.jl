using Test

include("quiet_solver_warnings.jl")
quiet_solver_warnings!()

const TEST_FILES = ("structure_tests.jl",
                    "charged_particle_3d_tests.jl",
                    "guiding_center_3d_tests.jl",
                    "guiding_center_4d_tests.jl",
                    "gyro_kinetics_4d_tests.jl",
                    "pauli_particle_3d_tests.jl",
                    "model_agreement_tests.jl",
                    "plots_tests.jl")

for f in TEST_FILES
    current_test_file!(f)
    include(f)
end

# Name any file that warned, and say how much, *before* asserting. A failing `@testset` throws at
# its end, so anything reported after it would be unreachable in exactly the case where it is needed.
for (f, n) in sort!(collect(suppressed_warning_counts()))
    n > 0 && @info "Suppressed $n nonlinear-solver warnings from $f " *
                   "(see test/quiet_solver_warnings.jl)"
end

# Every file is expected to be silent, and a nonlinear-solver warning from any of them is a
# regression — a tolerance that has slipped below a residual floor, or a time step at which the
# integrators stop converging. See the header of `quiet_solver_warnings.jl` for why this is
# assertable at all, and for what keeps the count at zero.
@testset "Nonlinear solver stays quiet                                                                        " begin
    for f in TEST_FILES
        @test get(suppressed_warning_counts(), f, 0) == 0
    end
end
