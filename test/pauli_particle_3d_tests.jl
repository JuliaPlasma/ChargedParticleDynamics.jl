
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


@safetestset "Pauli Particle in 3D in Tokamak                                                                     " begin

    using ChargedParticleDynamics.PauliParticle3d.PauliParticle3dTokamak
    using ..PauliParticle3dTests

    test_pauli_particle_3d(PauliParticle3dTokamak.pauli_particle_3d_pode())

end
