

module GuidingCenter4dTests

    using GeometricIntegrators

    const Δt = 1.
    const nt = 1

    export test_guiding_center_4d_glrk, test_guiding_center_4d_vpglrk

    function test_guiding_center_4d_glrk(ode; Δt = Δt)
        int = Integrator(ode, getTableauGLRK(1), Δt)
        sol = integrate(int, nt)
    end

    function test_guiding_center_4d_vpglrk(ode; Δt = Δt)
        int = Integrator(ode, getTableauVPGLRK(1), Δt)
        sol = integrate(int, nt)
    end

end


using .GuidingCenter4dTests

import ChargedParticleDynamics.GuidingCenter4d.TokamakSlowDeeplyPassing
import ChargedParticleDynamics.GuidingCenter4d.TokamakSlowDeeplyTrapped
import ChargedParticleDynamics.GuidingCenter4d.TokamakSlowBarelyPassing
import ChargedParticleDynamics.GuidingCenter4d.TokamakSlowBarelyTrapped

import ChargedParticleDynamics.GuidingCenter4d.TokamakFastDeeplyPassing
import ChargedParticleDynamics.GuidingCenter4d.TokamakFastDeeplyTrapped
import ChargedParticleDynamics.GuidingCenter4d.TokamakFastBarelyPassing
import ChargedParticleDynamics.GuidingCenter4d.TokamakFastBarelyTrapped

test_guiding_center_4d_glrk(TokamakSlowDeeplyPassing.guiding_center_4d_ode(), Δt=400.)
test_guiding_center_4d_glrk(TokamakSlowDeeplyTrapped.guiding_center_4d_ode(), Δt=400.)
test_guiding_center_4d_glrk(TokamakSlowBarelyPassing.guiding_center_4d_ode(), Δt=400.)
test_guiding_center_4d_glrk(TokamakSlowBarelyTrapped.guiding_center_4d_ode(), Δt=400.)

test_guiding_center_4d_glrk(TokamakFastDeeplyPassing.guiding_center_4d_ode(), Δt=1.)
test_guiding_center_4d_glrk(TokamakFastDeeplyTrapped.guiding_center_4d_ode(), Δt=1.)
test_guiding_center_4d_glrk(TokamakFastBarelyPassing.guiding_center_4d_ode(), Δt=1.)
test_guiding_center_4d_glrk(TokamakFastBarelyTrapped.guiding_center_4d_ode(), Δt=1.)

test_guiding_center_4d_vpglrk(TokamakSlowDeeplyPassing.guiding_center_4d_iode(), Δt=400.)
test_guiding_center_4d_vpglrk(TokamakSlowDeeplyTrapped.guiding_center_4d_iode(), Δt=400.)
test_guiding_center_4d_vpglrk(TokamakSlowBarelyPassing.guiding_center_4d_iode(), Δt=400.)
test_guiding_center_4d_vpglrk(TokamakSlowBarelyTrapped.guiding_center_4d_iode(), Δt=400.)

test_guiding_center_4d_vpglrk(TokamakFastDeeplyPassing.guiding_center_4d_iode(), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakFastDeeplyTrapped.guiding_center_4d_iode(), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakFastBarelyPassing.guiding_center_4d_iode(), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakFastBarelyTrapped.guiding_center_4d_iode(), Δt=1.)
