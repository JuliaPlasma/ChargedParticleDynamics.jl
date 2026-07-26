
using SafeTestsets

module PauliParticle3dTests

using GeometricIntegrators
using Test

export test_pauli_particle_3d

# See `guiding_center_3d_tests.jl` for why `f_abstol` has to stay above the residual's round-off
# floor and why `f_reltol` is left at its default.
const options = (f_abstol=1E-12, max_iterations=50, warn_iterations=50)

# Asserts that the integration returns a solution rather than `@test_nowarn`; see the comment on
# the corresponding functions in `guiding_center_3d_tests.jl` for why.
function test_pauli_particle_3d(equ::Union{HODEProblem,PODEProblem}; kwargs...)
    @test integrate(equ, PartitionedGauss(2); options..., kwargs...) isa GeometricSolution
end

function test_pauli_particle_3d(equ::Union{IODEProblem,LODEProblem}; kwargs...)
    @test integrate(equ, VPRKGauss(2); options..., kwargs...) isa GeometricSolution
end

end


@safetestset "Pauli Particle in 3D with ITER-like Solov'ev Equilibrium                                            " begin

    using ChargedParticleDynamics.PauliParticle3d.SolovevIter
    using ..PauliParticle3dTests

    test_pauli_particle_3d(SolovevIter.hodeproblem(initial_conditions_barely_passing(); timespan=(0.0, 1E2), timestep=0.1))
    test_pauli_particle_3d(SolovevIter.hodeproblem(initial_conditions_barely_trapped(); timespan=(0.0, 1E2), timestep=0.1))
    test_pauli_particle_3d(SolovevIter.hodeproblem(initial_conditions_deeply_passing(); timespan=(0.0, 1E2), timestep=0.1))
    test_pauli_particle_3d(SolovevIter.hodeproblem(initial_conditions_deeply_trapped(); timespan=(0.0, 1E2), timestep=0.1))

    test_pauli_particle_3d(SolovevIter.iodeproblem(initial_conditions_barely_passing(); timespan=(0.0, 1E2), timestep=0.1))
    test_pauli_particle_3d(SolovevIter.iodeproblem(initial_conditions_barely_trapped(); timespan=(0.0, 1E2), timestep=0.1))
    test_pauli_particle_3d(SolovevIter.iodeproblem(initial_conditions_deeply_passing(); timespan=(0.0, 1E2), timestep=0.1))
    test_pauli_particle_3d(SolovevIter.iodeproblem(initial_conditions_deeply_trapped(); timespan=(0.0, 1E2), timestep=0.1))

end


@safetestset "Pauli Particle in 3D with ITER-like Solov'ev Equilibrium with X-Point                               " begin

    using ChargedParticleDynamics.PauliParticle3d.SolovevIterXpoint
    using ..PauliParticle3dTests

    test_pauli_particle_3d(SolovevIterXpoint.hodeproblem(initial_conditions_barely_passing(); timespan=(0.0, 1E2), timestep=0.1))
    test_pauli_particle_3d(SolovevIterXpoint.hodeproblem(initial_conditions_barely_trapped(); timespan=(0.0, 1E2), timestep=0.1))
    test_pauli_particle_3d(SolovevIterXpoint.hodeproblem(initial_conditions_deeply_passing(); timespan=(0.0, 1E2), timestep=0.1))
    test_pauli_particle_3d(SolovevIterXpoint.hodeproblem(initial_conditions_deeply_trapped(); timespan=(0.0, 1E2), timestep=0.1))

    test_pauli_particle_3d(SolovevIterXpoint.iodeproblem(initial_conditions_barely_passing(); timespan=(0.0, 1E2), timestep=0.1))
    test_pauli_particle_3d(SolovevIterXpoint.iodeproblem(initial_conditions_barely_trapped(); timespan=(0.0, 1E2), timestep=0.1))
    test_pauli_particle_3d(SolovevIterXpoint.iodeproblem(initial_conditions_deeply_passing(); timespan=(0.0, 1E2), timestep=0.1))
    test_pauli_particle_3d(SolovevIterXpoint.iodeproblem(initial_conditions_deeply_trapped(); timespan=(0.0, 1E2), timestep=0.1))

end


@safetestset "Pauli Particle in 3D in Symmetric Magnetic Field                                                    " begin

    using ChargedParticleDynamics.PauliParticle3d.SymmetricField
    using ..PauliParticle3dTests

    test_pauli_particle_3d(SymmetricField.podeproblem())
    test_pauli_particle_3d(SymmetricField.hodeproblem())
    test_pauli_particle_3d(SymmetricField.iodeproblem())

end


@safetestset "Pauli Particle in 3D in Theta Pinch                                                                 " begin

    using GeometricIntegrators: MidpointExtrapolation
    using ChargedParticleDynamics.PauliParticle3d.ThetaPinchField
    using ..PauliParticle3dTests

    # the momentum is constant along the initial trajectory, so the default Hermite
    # extrapolation of the initial guess is degenerate
    test_pauli_particle_3d(ThetaPinchField.podeproblem(); initialguess=MidpointExtrapolation(5))
    test_pauli_particle_3d(ThetaPinchField.hodeproblem(); initialguess=MidpointExtrapolation(5))
    test_pauli_particle_3d(ThetaPinchField.iodeproblem(); initialguess=MidpointExtrapolation(5))

end


@safetestset "Pauli Particle in 3D in ITER Equilibrium in Cylindrical Coordinates                                 " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakIterCylindrical
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakIterCylindrical.podeproblem())
    test_pauli_particle_3d(TokamakIterCylindrical.hodeproblem())
    test_pauli_particle_3d(TokamakIterCylindrical.iodeproblem())

end

@safetestset "Pauli Particle in 3D in Tokamak Equilibrium in Cartesian Coordinates                                " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakSmallCartesian
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakSmallCartesian.podeproblem(timespan=(0.0, 1E3), timestep=10.0))
    test_pauli_particle_3d(TokamakSmallCartesian.hodeproblem(timespan=(0.0, 1E3), timestep=10.0))
    test_pauli_particle_3d(TokamakSmallCartesian.iodeproblem())

end

@safetestset "Pauli Particle in 3D in Tokamak Equilibrium in Cylindrical Coordinates                              " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakSmallCylindrical
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakSmallCylindrical.podeproblem(timespan=(0.0, 1E3), timestep=10.0))
    test_pauli_particle_3d(TokamakSmallCylindrical.hodeproblem(timespan=(0.0, 1E3), timestep=10.0))
    test_pauli_particle_3d(TokamakSmallCylindrical.iodeproblem(timespan=(0.0, 1E3), timestep=10.0))

end

@safetestset "Pauli Particle in 3D in Tokamak Equilibrium in Toroidal Coordinates                                 " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakSmallToroidal
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakSmallToroidal.podeproblem(timespan=(0.0, 1E3), timestep=10.0))
    test_pauli_particle_3d(TokamakSmallToroidal.hodeproblem(timespan=(0.0, 1E3), timestep=10.0))
    test_pauli_particle_3d(TokamakSmallToroidal.iodeproblem(timespan=(0.0, 1E3), timestep=10.0))

end
