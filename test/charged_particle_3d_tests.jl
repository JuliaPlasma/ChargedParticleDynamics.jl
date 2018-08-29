

module ChargedParticle3dTests

    using GeometricIntegrators

    const Δt = 0.1
    const nt = 1

    export test_charged_particle_3d

    function test_charged_particle_3d(ode)
        int = Integrator(ode, getTableauVPGLRK(1), Δt)
        sol = integrate(int, nt)
    end

end


using .ChargedParticle3dTests

import ChargedParticleDynamics.ChargedParticle3d.ChargedParticle3dUniform
import ChargedParticleDynamics.ChargedParticle3d.ChargedParticle3dSymmetric
import ChargedParticleDynamics.ChargedParticle3d.ChargedParticle3dSingular

test_charged_particle_3d(ChargedParticle3dUniform.charged_particle_3d_iode())
test_charged_particle_3d(ChargedParticle3dSymmetric.charged_particle_3d_iode())
test_charged_particle_3d(ChargedParticle3dSingular.charged_particle_3d_iode())
