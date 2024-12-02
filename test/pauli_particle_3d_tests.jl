
using SafeTestsets

module PauliParticle3dTests

    using Test
    using GeometricIntegrators

    export test_pauli_particle_3d

    function test_pauli_particle_3d(ode::PODEProblem)
        sol = integrate(ode, PartitionedGauss(2))
        @test true
    end

    function test_pauli_particle_3d(ode::IODEProblem)
        sol = integrate(ode, VPRKGauss(2))
        @test true
    end

end



@safetestset "Pauli Particle in 3D in Symmetric Magnetic Field                                                    " begin

    using ChargedParticleDynamics.PauliParticle3d.SymmetricField
    using ..PauliParticle3dTests

    test_pauli_particle_3d(SymmetricField.pauli_particle_3d_pode(tspan = (0., 10.)))
    test_pauli_particle_3d(SymmetricField.pauli_particle_3d_iode(tspan = (0., 10.)))

end


@safetestset "Pauli Particle in 3D in Theta Pinch                                                                 " begin

    using ChargedParticleDynamics.PauliParticle3d.ThetaPinchField
    using ..PauliParticle3dTests

    test_pauli_particle_3d(ThetaPinchField.pauli_particle_3d_pode(tspan = (0., 10.)))
    test_pauli_particle_3d(ThetaPinchField.pauli_particle_3d_iode(tspan = (0., 10.)))

end


@safetestset "Pauli Particle in 3D in ITER Equilibrium in Cylindrical Coordinates                                 " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakIterCylindrical
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakIterCylindrical.pauli_particle_3d_pode(tspan = (0., 10.)))
    test_pauli_particle_3d(TokamakIterCylindrical.pauli_particle_3d_iode(tspan = (0., 10.)))

end

@safetestset "Pauli Particle in 3D in Tokamak Equilibrium in Cartesian Coordinates                                " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakSmallCartesian
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakSmallCartesian.pauli_particle_3d_pode(tspan = (0., 4E3)))
    test_pauli_particle_3d(TokamakSmallCartesian.pauli_particle_3d_iode(tspan = (0., 4E3)))

end

@safetestset "Pauli Particle in 3D in Tokamak Equilibrium in Cylindrical Coordinates                              " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakSmallCylindrical
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakSmallCylindrical.pauli_particle_3d_pode(tspan = (0., 10.), tstep = 1.0))
    test_pauli_particle_3d(TokamakSmallCylindrical.pauli_particle_3d_iode(tspan = (0., 4E3)))

end

@safetestset "Pauli Particle in 3D in Tokamak Equilibrium in Toroidal Coordinates                                 " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakSmallToroidal
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakSmallToroidal.pauli_particle_3d_pode(tspan = (0., 10.), tstep = 1.0))
    test_pauli_particle_3d(TokamakSmallToroidal.pauli_particle_3d_iode(tspan = (0., 4E3)))

end
