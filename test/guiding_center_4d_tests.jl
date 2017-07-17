

module GuidingCenter4dTests

    using GeometricIntegrators

    const Δt = 400.
    const nt = 1

    export test_guiding_center_4d_glrk, test_guiding_center_4d_vpglrk

    function test_guiding_center_4d_glrk(ode)
        int = Integrator(ode, getTableauGLRK(1), Δt)
        sol = integrate(int, nt)
    end

    function test_guiding_center_4d_vpglrk(ode)
        int = Integrator(ode, getTableauVPGLRK(1), Δt)
        sol = integrate(int, nt)
    end

end


using GuidingCenter4dTests

import ChargedParticleDynamics.GuidingCenter4d.TokamakPassing
import ChargedParticleDynamics.GuidingCenter4d.TokamakTrapped
import ChargedParticleDynamics.GuidingCenter4d.TokamakBarelyTrapped

test_guiding_center_4d_glrk(ChargedParticleDynamics.GuidingCenter4d.TokamakPassing.guiding_center_4d_ode())
test_guiding_center_4d_glrk(ChargedParticleDynamics.GuidingCenter4d.TokamakTrapped.guiding_center_4d_ode())
test_guiding_center_4d_glrk(ChargedParticleDynamics.GuidingCenter4d.TokamakBarelyTrapped.guiding_center_4d_ode())

test_guiding_center_4d_vpglrk(ChargedParticleDynamics.GuidingCenter4d.TokamakPassing.guiding_center_4d_iode())
test_guiding_center_4d_vpglrk(ChargedParticleDynamics.GuidingCenter4d.TokamakTrapped.guiding_center_4d_iode())
test_guiding_center_4d_vpglrk(ChargedParticleDynamics.GuidingCenter4d.TokamakBarelyTrapped.guiding_center_4d_iode())
