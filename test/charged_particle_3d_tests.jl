
using SafeTestsets

module ChargedParticle3dTests

    using Test
    using GeometricIntegrators

    const Δt = 0.01
    const nt = 1

    export test_charged_particle_3d

    function test_charged_particle_3d(ode::ODE)
        int = Integrator(ode, getTableauGLRK(1), Δt)
        sol = integrate(ode, int, nt)
        @test true
    end

    function test_charged_particle_3d(ode::PODE)
        int = Integrator(ode, TableauIPRK(:pglrk, 2, getCoefficientsGLRK(1)), Δt)
        sol = integrate(ode, int, nt)
        @test true
    end

    function test_charged_particle_3d(ode::IODE)
        int = Integrator(ode, getTableauVPGLRK(1), Δt)
        sol = integrate(ode, int, nt)
        @test true
    end

end



@safetestset "Charged Particle Dynamics in 3D in Uniform Magnetic Field                                           " begin

    using ChargedParticleDynamics.ChargedParticle3d.ChargedParticle3dUniform
    using ..ChargedParticle3dTests

    test_charged_particle_3d(ChargedParticle3dUniform.charged_particle_3d_iode())

end


@safetestset "Charged Particle Dynamics in 3D in Symmetric Magnetic Field                                         " begin

    using ChargedParticleDynamics.ChargedParticle3d.ChargedParticle3dSymmetric
    using ..ChargedParticle3dTests

    test_charged_particle_3d(ChargedParticle3dSymmetric.charged_particle_3d_iode())

end


@safetestset "Charged Particle Dynamics in 3D in Singular Magnetic Field                                          " begin

    using ChargedParticleDynamics.ChargedParticle3d.ChargedParticle3dSingular
    using ..ChargedParticle3dTests

    test_charged_particle_3d(ChargedParticle3dSingular.charged_particle_3d_iode())

end


@safetestset "Charged Particle Dynamics in 3D in Theta Pinch                                                      " begin

    using ChargedParticleDynamics.ChargedParticle3d.ThetaPinchNoncanonical
    using ..ChargedParticle3dTests

    test_charged_particle_3d(ThetaPinchNoncanonical.charged_particle_3d_iode())

end


@safetestset "Charged Particle Dynamics in 3D in Tokamak (cartesian coordinates)                                  " begin

    using ChargedParticleDynamics.ChargedParticle3d.TokamakSmallCartesian
    using ..ChargedParticle3dTests

    test_charged_particle_3d(TokamakSmallCartesian.charged_particle_3d_pode())

end


@safetestset "Charged Particle Dynamics in 3D in Tokamak (cylindrical coordinates)                                " begin

    using ChargedParticleDynamics.ChargedParticle3d.TokamakSmallCylindrical
    using ..ChargedParticle3dTests

    test_charged_particle_3d(TokamakSmallCylindrical.charged_particle_3d_pode())

end


@safetestset "Charged Particle Dynamics in 3D in Tokamak (toroidal coordinates)                                   " begin

    using ChargedParticleDynamics.ChargedParticle3d.TokamakSmallToroidal
    using ..ChargedParticle3dTests

    test_charged_particle_3d(TokamakSmallToroidal.charged_particle_3d_pode())

end


@safetestset "Charged Particle Dynamics in 3D in Tokamak (noncanonical formulation)                               " begin

    using ChargedParticleDynamics.ChargedParticle3d.TokamakSmallNoncanonical
    using ..ChargedParticle3dTests

    test_charged_particle_3d(TokamakSmallNoncanonical.charged_particle_3d_iode())

end
