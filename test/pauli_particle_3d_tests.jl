
using SafeTestsets

module PauliParticle3dTests

using GeometricIntegrators
using Test

export test_pauli_particle_3d

const options = (f_abstol=1E-15, f_reltol=1E-15)

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

    test_pauli_particle_3d(SolovevIter.pauli_particle_3d_hode(initial_conditions_barely_passing()...; tspan=(0.0, 1E2), tstep=0.1))
    test_pauli_particle_3d(SolovevIter.pauli_particle_3d_hode(initial_conditions_barely_trapped()...; tspan=(0.0, 1E2), tstep=0.1))
    test_pauli_particle_3d(SolovevIter.pauli_particle_3d_hode(initial_conditions_deeply_passing()...; tspan=(0.0, 1E2), tstep=0.1))
    test_pauli_particle_3d(SolovevIter.pauli_particle_3d_hode(initial_conditions_deeply_trapped()...; tspan=(0.0, 1E2), tstep=0.1))

    test_pauli_particle_3d(SolovevIter.pauli_particle_3d_iode(initial_conditions_barely_passing()...; tspan=(0.0, 1E2), tstep=0.1))
    test_pauli_particle_3d(SolovevIter.pauli_particle_3d_iode(initial_conditions_barely_trapped()...; tspan=(0.0, 1E2), tstep=0.1))
    test_pauli_particle_3d(SolovevIter.pauli_particle_3d_iode(initial_conditions_deeply_passing()...; tspan=(0.0, 1E2), tstep=0.1))
    test_pauli_particle_3d(SolovevIter.pauli_particle_3d_iode(initial_conditions_deeply_trapped()...; tspan=(0.0, 1E2), tstep=0.1))

end


@safetestset "Pauli Particle in 3D with ITER-like Solov'ev Equilibrium with X-Point                               " begin

    using ChargedParticleDynamics.PauliParticle3d.SolovevIterXpoint
    using ..PauliParticle3dTests

    test_pauli_particle_3d(SolovevIterXpoint.pauli_particle_3d_hode(initial_conditions_barely_passing()...; tspan=(0.0, 1E2), tstep=0.1))
    test_pauli_particle_3d(SolovevIterXpoint.pauli_particle_3d_hode(initial_conditions_barely_trapped()...; tspan=(0.0, 1E2), tstep=0.1))
    test_pauli_particle_3d(SolovevIterXpoint.pauli_particle_3d_hode(initial_conditions_deeply_passing()...; tspan=(0.0, 1E2), tstep=0.1))
    test_pauli_particle_3d(SolovevIterXpoint.pauli_particle_3d_hode(initial_conditions_deeply_trapped()...; tspan=(0.0, 1E2), tstep=0.1))

    test_pauli_particle_3d(SolovevIterXpoint.pauli_particle_3d_iode(initial_conditions_barely_passing()...; tspan=(0.0, 1E2), tstep=0.1))
    test_pauli_particle_3d(SolovevIterXpoint.pauli_particle_3d_iode(initial_conditions_barely_trapped()...; tspan=(0.0, 1E2), tstep=0.1))
    test_pauli_particle_3d(SolovevIterXpoint.pauli_particle_3d_iode(initial_conditions_deeply_passing()...; tspan=(0.0, 1E2), tstep=0.1))
    test_pauli_particle_3d(SolovevIterXpoint.pauli_particle_3d_iode(initial_conditions_deeply_trapped()...; tspan=(0.0, 1E2), tstep=0.1))

end


@safetestset "Pauli Particle in 3D in Symmetric Magnetic Field                                                    " begin

    using ChargedParticleDynamics.PauliParticle3d.SymmetricField
    using ..PauliParticle3dTests

    test_pauli_particle_3d(SymmetricField.pauli_particle_3d_pode())
    test_pauli_particle_3d(SymmetricField.pauli_particle_3d_hode())
    test_pauli_particle_3d(SymmetricField.pauli_particle_3d_iode())

end


@safetestset "Pauli Particle in 3D in Theta Pinch                                                                 " begin

    using GeometricIntegrators: MidpointExtrapolation
    using ChargedParticleDynamics.PauliParticle3d.ThetaPinchField
    using ..PauliParticle3dTests

    # the momentum is constant along the initial trajectory, so the default Hermite
    # extrapolation of the initial guess is degenerate
    test_pauli_particle_3d(ThetaPinchField.pauli_particle_3d_pode(); initialguess=MidpointExtrapolation(5))
    test_pauli_particle_3d(ThetaPinchField.pauli_particle_3d_hode(); initialguess=MidpointExtrapolation(5))
    test_pauli_particle_3d(ThetaPinchField.pauli_particle_3d_iode(); initialguess=MidpointExtrapolation(5))

end


@safetestset "Pauli Particle in 3D in ITER Equilibrium in Cylindrical Coordinates                                 " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakIterCylindrical
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakIterCylindrical.pauli_particle_3d_pode())
    test_pauli_particle_3d(TokamakIterCylindrical.pauli_particle_3d_hode())
    test_pauli_particle_3d(TokamakIterCylindrical.pauli_particle_3d_iode())

end

@safetestset "Pauli Particle in 3D in Tokamak Equilibrium in Cartesian Coordinates                                " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakSmallCartesian
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakSmallCartesian.pauli_particle_3d_pode(tspan=(0.0, 1E3), tstep=10.0))
    test_pauli_particle_3d(TokamakSmallCartesian.pauli_particle_3d_hode(tspan=(0.0, 1E3), tstep=10.0))
    test_pauli_particle_3d(TokamakSmallCartesian.pauli_particle_3d_iode())

end

@safetestset "Pauli Particle in 3D in Tokamak Equilibrium in Cylindrical Coordinates                              " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakSmallCylindrical
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakSmallCylindrical.pauli_particle_3d_pode(tspan=(0.0, 1E3), tstep=10.0))
    test_pauli_particle_3d(TokamakSmallCylindrical.pauli_particle_3d_hode(tspan=(0.0, 1E3), tstep=10.0))
    test_pauli_particle_3d(TokamakSmallCylindrical.pauli_particle_3d_iode(tspan=(0.0, 1E3), tstep=10.0))

end

@safetestset "Pauli Particle in 3D in Tokamak Equilibrium in Toroidal Coordinates                                 " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakSmallToroidal
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakSmallToroidal.pauli_particle_3d_pode(tspan=(0.0, 1E3), tstep=10.0))
    test_pauli_particle_3d(TokamakSmallToroidal.pauli_particle_3d_hode(tspan=(0.0, 1E3), tstep=10.0))
    test_pauli_particle_3d(TokamakSmallToroidal.pauli_particle_3d_iode(tspan=(0.0, 1E3), tstep=10.0))

end
