
using SafeTestsets

module PauliParticle3dTests

    using Test
    using GeometricIntegrators

    const Δt = 0.1
    const nt = 1

    export test_pauli_particle_3d

    function test_pauli_particle_3d(ode::PODE)
        int = Integrator(ode, TableauIPRK(:pglrk, 2, getCoefficientsGLRK(1)), Δt)
        sol = integrate(ode, int, nt)
        @test true
    end

end



@safetestset "Pauli Particle in 3D in Uniform Magnetic Field                                                      " begin

    using ChargedParticleDynamics.PauliParticle3d.PauliParticle3dUniform
    using ..PauliParticle3dTests

    test_pauli_particle_3d(PauliParticle3dUniform.pauli_particle_3d_pode())

end


@safetestset "Pauli Particle in 3D in Symmetric Magnetic Field                                                    " begin

    using ChargedParticleDynamics.PauliParticle3d.PauliParticle3dSymmetric
    using ..PauliParticle3dTests

    test_pauli_particle_3d(PauliParticle3dSymmetric.pauli_particle_3d_pode())

end


@safetestset "Pauli Particle in 3D in ITER Equilibrium in Cylindrical Coordinates                                 " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakIterCylindrical
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakIterCylindrical.pauli_particle_3d_pode())

end

@safetestset "Pauli Particle in 3D in Tokamak Equilibrium in Cartesian Coordinates                                " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakSmallCartesian
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakSmallCartesian.pauli_particle_3d_pode())

end

@safetestset "Pauli Particle in 3D in Tokamak Equilibrium in Circular Coordinates                                 " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakSmallCircular
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakSmallCircular.pauli_particle_3d_pode())

end

@safetestset "Pauli Particle in 3D in Tokamak Equilibrium in Cylindrical Coordinates                              " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakSmallCylindrical
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakSmallCylindrical.pauli_particle_3d_pode())

end
