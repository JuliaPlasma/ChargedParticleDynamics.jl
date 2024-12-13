
using SafeTestsets

module PauliParticle3dTests

    using GeometricIntegrators
    using SimpleSolvers: Options
    using Test

    export test_pauli_particle_3d

    const options = Options(x_reltol = 1E-14, f_abstol = 1E-14, f_reltol = 1E-14)

    function test_pauli_particle_3d(ode::Union{HODEProblem,PODEProblem})
        @test_nowarn integrate(ode, PartitionedGauss(2); options = options)
    end

    function test_pauli_particle_3d(ode::Union{IODEProblem,LODEProblem})
        @test_nowarn integrate(ode, VPRKGauss(2); options = options)
    end

end



@safetestset "Pauli Particle in 3D in Symmetric Magnetic Field                                                    " begin

    using ChargedParticleDynamics.PauliParticle3d.SymmetricField
    using ..PauliParticle3dTests

    test_pauli_particle_3d(SymmetricField.pauli_particle_3d_pode(tspan = (0., 10.)))
    test_pauli_particle_3d(SymmetricField.pauli_particle_3d_hode(tspan = (0., 10.)))
    test_pauli_particle_3d(SymmetricField.pauli_particle_3d_iode())

end


@safetestset "Pauli Particle in 3D in Theta Pinch                                                                 " begin

    using ChargedParticleDynamics.PauliParticle3d.ThetaPinchField
    using ..PauliParticle3dTests

    test_pauli_particle_3d(ThetaPinchField.pauli_particle_3d_iode())

end


@safetestset "Pauli Particle in 3D in ITER Equilibrium in Cylindrical Coordinates                                 " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakIterCylindrical
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakIterCylindrical.pauli_particle_3d_pode(tspan = (0., 10.)))
    test_pauli_particle_3d(TokamakIterCylindrical.pauli_particle_3d_hode(tspan = (0., 10.)))
    test_pauli_particle_3d(TokamakIterCylindrical.pauli_particle_3d_iode())

end

@safetestset "Pauli Particle in 3D in Tokamak Equilibrium in Cartesian Coordinates                                " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakSmallCartesian
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakSmallCartesian.pauli_particle_3d_pode(tspan = (0., 200)))
    test_pauli_particle_3d(TokamakSmallCartesian.pauli_particle_3d_hode(tspan = (0., 200)))
    test_pauli_particle_3d(TokamakSmallCartesian.pauli_particle_3d_iode())

end

@safetestset "Pauli Particle in 3D in Tokamak Equilibrium in Cylindrical Coordinates                              " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakSmallCylindrical
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakSmallCylindrical.pauli_particle_3d_pode(tspan = (0., 200.)))
    test_pauli_particle_3d(TokamakSmallCylindrical.pauli_particle_3d_hode(tspan = (0., 200.)))
    test_pauli_particle_3d(TokamakSmallCylindrical.pauli_particle_3d_iode())

end

@safetestset "Pauli Particle in 3D in Tokamak Equilibrium in Toroidal Coordinates                                 " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakSmallToroidal
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakSmallToroidal.pauli_particle_3d_pode(tspan = (0., 200.)))
    test_pauli_particle_3d(TokamakSmallToroidal.pauli_particle_3d_hode(tspan = (0., 200.)))
    test_pauli_particle_3d(TokamakSmallToroidal.pauli_particle_3d_iode())

end
