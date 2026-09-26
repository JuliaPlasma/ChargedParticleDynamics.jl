
using SafeTestsets
using Test

include("helpers/quiet_solver_warnings.jl")
quiet_solver_warnings!()
current_test_file!("ChargedParticle3d.jl")

module ChargedParticle3dTests

using GeometricIntegrators
using Test

export test_charged_particle_3d

# Asserts that the integration returns a solution rather than `@test_nowarn`; see the comment
# on the corresponding functions in `GuidingCenter3d.jl` for why.
function test_charged_particle_3d(equ::ODEProblem)
    @test integrate(equ, Gauss(1)) isa GeometricSolution
end

function test_charged_particle_3d(equ::Union{HODEProblem, PODEProblem})
    @test integrate(equ, PartitionedGauss(1)) isa GeometricSolution
end

function test_charged_particle_3d(equ::Union{IODEProblem, LODEProblem})
    @test integrate(equ, VPRKGauss(1)) isa GeometricSolution
end

end

@safetestset "Charged Particle Dynamics in 3D in Symmetric Magnetic Field                                         " begin
    using ChargedParticleDynamics.ChargedParticle3d.SymmetricField
    using ..ChargedParticle3dTests

    test_charged_particle_3d(SymmetricField.iodeproblem())
end

@safetestset "Charged Particle Dynamics in 3D in Singular Magnetic Field                                          " begin
    using ChargedParticleDynamics.ChargedParticle3d.SingularField
    using ..ChargedParticle3dTests

    test_charged_particle_3d(SingularField.iodeproblem())
end

@safetestset "Charged Particle Dynamics in 3D in Theta Pinch                                                      " begin
    using ChargedParticleDynamics.ChargedParticle3d.ThetaPinchCanonical
    using ..ChargedParticle3dTests

    test_charged_particle_3d(ThetaPinchCanonical.iodeproblem())
end

@safetestset "Charged Particle Dynamics in 3D in Theta Pinch (noncanonical formulation)                           " begin
    using ChargedParticleDynamics.ChargedParticle3d.ThetaPinchNoncanonical
    using ..ChargedParticle3dTests

    test_charged_particle_3d(ThetaPinchNoncanonical.iodeproblem())
end

@safetestset "Charged Particle Dynamics in 3D in Tokamak (cartesian coordinates)                                  " begin
    using ChargedParticleDynamics.ChargedParticle3d.TokamakSmallCartesian
    using ..ChargedParticle3dTests

    test_charged_particle_3d(TokamakSmallCartesian.podeproblem())
    test_charged_particle_3d(TokamakSmallCartesian.iodeproblem())
end

@safetestset "Charged Particle Dynamics in 3D in Tokamak (cylindrical coordinates)                                " begin
    using ChargedParticleDynamics.ChargedParticle3d.TokamakSmallCylindrical
    using ..ChargedParticle3dTests

    test_charged_particle_3d(TokamakSmallCylindrical.podeproblem())
    test_charged_particle_3d(TokamakSmallCylindrical.iodeproblem())
end

@safetestset "Charged Particle Dynamics in 3D in Tokamak (toroidal coordinates)                                   " begin
    using ChargedParticleDynamics.ChargedParticle3d.TokamakSmallToroidal
    using ..ChargedParticle3dTests

    test_charged_particle_3d(TokamakSmallToroidal.podeproblem())
    test_charged_particle_3d(TokamakSmallToroidal.iodeproblem())
end

@safetestset "Charged Particle Dynamics in 3D in Tokamak (noncanonical formulation)                               " begin
    using ChargedParticleDynamics.ChargedParticle3d.TokamakSmallNoncanonical
    using ..ChargedParticle3dTests

    test_charged_particle_3d(TokamakSmallNoncanonical.iodeproblem())
end

# This file is expected to be silent, and a nonlinear-solver warning from it is a regression — a
# tolerance that has slipped below a residual floor, or a time step at which the integrators stop
# converging. See the header of `helpers/quiet_solver_warnings.jl` for why this is assertable at
# all, and for what keeps the count at zero. The count is reported *before* asserting, because a
# failing `@testset` throws at its end.
suppressed_warning_count() > 0 &&
    @info "Suppressed $(suppressed_warning_count()) nonlinear-solver warnings " *
          "(see test/helpers/quiet_solver_warnings.jl)" suppressed_warning_counts()

# For non-negative counts the two assertions are one statement; see K1 in `KNOWN_ISSUES.md`.
@testset "Nonlinear solver stays quiet" begin
    @test all(iszero, values(suppressed_warning_counts()))
    @test suppressed_warning_count() == 0
end
