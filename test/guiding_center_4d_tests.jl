

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


using .GuidingCenter4dTests

import ChargedParticleDynamics.GuidingCenter4d.TokamakSlowDeeplyPassing
import ChargedParticleDynamics.GuidingCenter4d.TokamakSlowDeeplyTrapped
import ChargedParticleDynamics.GuidingCenter4d.TokamakSlowBarelyPassing
import ChargedParticleDynamics.GuidingCenter4d.TokamakSlowBarelyTrapped

test_guiding_center_4d_glrk(TokamakSlowDeeplyPassing.guiding_center_4d_ode())
test_guiding_center_4d_glrk(TokamakSlowDeeplyTrapped.guiding_center_4d_ode())
test_guiding_center_4d_glrk(TokamakSlowBarelyPassing.guiding_center_4d_ode())
test_guiding_center_4d_glrk(TokamakSlowBarelyTrapped.guiding_center_4d_ode())

test_guiding_center_4d_vpglrk(TokamakSlowDeeplyPassing.guiding_center_4d_iode())
test_guiding_center_4d_vpglrk(TokamakSlowDeeplyTrapped.guiding_center_4d_iode())
test_guiding_center_4d_vpglrk(TokamakSlowBarelyPassing.guiding_center_4d_iode())
test_guiding_center_4d_vpglrk(TokamakSlowBarelyTrapped.guiding_center_4d_iode())
