
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



@safetestset "Charged Particle Dynamics in 3D in Symmetric Magnetic Field                                         " begin

    using ChargedParticleDynamics.ChargedParticle3d.SymmetricField
    using ..ChargedParticle3dTests

    test_charged_particle_3d(SymmetricField.charged_particle_3d_iode())

end


@safetestset "Charged Particle Dynamics in 3D in Singular Magnetic Field                                          " begin

    using ChargedParticleDynamics.ChargedParticle3d.SingularField
    using ..ChargedParticle3dTests

    test_charged_particle_3d(SingularField.charged_particle_3d_iode())

end


@safetestset "Charged Particle Dynamics in 3D in Theta Pinch                                                      " begin

    using ChargedParticleDynamics.ChargedParticle3d.ThetaPinchCanonical
    using ..ChargedParticle3dTests

    test_charged_particle_3d(ThetaPinchCanonical.charged_particle_3d_iode())

end


@safetestset "Charged Particle Dynamics in 3D in Theta Pinch (noncanonical formulation)                           " begin

    using ChargedParticleDynamics.ChargedParticle3d.ThetaPinchNoncanonical
    using ..ChargedParticle3dTests

    test_charged_particle_3d(ThetaPinchNoncanonical.charged_particle_3d_iode())

end


@safetestset "Charged Particle Dynamics in 3D in Tokamak (cartesian coordinates)                                  " begin

    using ChargedParticleDynamics.ChargedParticle3d.TokamakSmallCartesian
    using ..ChargedParticle3dTests

    test_charged_particle_3d(TokamakSmallCartesian.charged_particle_3d_pode())

end


@safetestset "Charged Particle Dynamics in 3D in Tokamak (circular coordinates)                                   " begin

    using ChargedParticleDynamics.ChargedParticle3d.TokamakSmallCircular
    using ..ChargedParticle3dTests

    test_charged_particle_3d(TokamakSmallCircular.charged_particle_3d_pode())

end


@safetestset "Charged Particle Dynamics in 3D in Tokamak (cylindrical coordinates)                                " begin

    using ChargedParticleDynamics.ChargedParticle3d.TokamakSmallCylindrical
    using ..ChargedParticle3dTests

    test_charged_particle_3d(TokamakSmallCylindrical.charged_particle_3d_pode())

end


@safetestset "Charged Particle Dynamics in 3D in Tokamak (noncanonical formulation)                               " begin

    using ChargedParticleDynamics.ChargedParticle3d.TokamakSmallNoncanonical
    using ..ChargedParticle3dTests

    test_charged_particle_3d(TokamakSmallNoncanonical.charged_particle_3d_iode())

end
