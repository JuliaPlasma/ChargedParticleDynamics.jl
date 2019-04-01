

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


nl = 100
nx = 10
ny = 10


using .GuidingCenter4dTests

import ChargedParticleDynamics.GuidingCenter4d.TokamakSlowDeeplyPassing
import ChargedParticleDynamics.GuidingCenter4d.TokamakSlowDeeplyTrapped
import ChargedParticleDynamics.GuidingCenter4d.TokamakSlowBarelyPassing
import ChargedParticleDynamics.GuidingCenter4d.TokamakSlowBarelyTrapped
import ChargedParticleDynamics.GuidingCenter4d.TokamakSlowLoop
import ChargedParticleDynamics.GuidingCenter4d.TokamakSlowSurface

import ChargedParticleDynamics.GuidingCenter4d.TokamakFastDeeplyPassing
import ChargedParticleDynamics.GuidingCenter4d.TokamakFastDeeplyTrapped
import ChargedParticleDynamics.GuidingCenter4d.TokamakFastBarelyPassing
import ChargedParticleDynamics.GuidingCenter4d.TokamakFastBarelyTrapped
import ChargedParticleDynamics.GuidingCenter4d.TokamakFastLoop
import ChargedParticleDynamics.GuidingCenter4d.TokamakFastSurface

import ChargedParticleDynamics.GuidingCenter4d.SymmetricLoop
import ChargedParticleDynamics.GuidingCenter4d.SymmetricSurface
import ChargedParticleDynamics.GuidingCenter4d.UniformLoop

test_guiding_center_4d_glrk(TokamakSlowDeeplyPassing.guiding_center_4d_ode(), Δt=400.)
test_guiding_center_4d_glrk(TokamakSlowDeeplyTrapped.guiding_center_4d_ode(), Δt=400.)
test_guiding_center_4d_glrk(TokamakSlowBarelyPassing.guiding_center_4d_ode(), Δt=400.)
test_guiding_center_4d_glrk(TokamakSlowBarelyTrapped.guiding_center_4d_ode(), Δt=400.)
test_guiding_center_4d_glrk(TokamakSlowLoop.guiding_center_4d_loop_ode(nl), Δt=400.)
test_guiding_center_4d_glrk(TokamakSlowSurface.guiding_center_4d_surface_ode(nx, ny), Δt=400.)

test_guiding_center_4d_glrk(TokamakFastDeeplyPassing.guiding_center_4d_ode(), Δt=1.)
test_guiding_center_4d_glrk(TokamakFastDeeplyTrapped.guiding_center_4d_ode(), Δt=1.)
test_guiding_center_4d_glrk(TokamakFastBarelyPassing.guiding_center_4d_ode(), Δt=1.)
test_guiding_center_4d_glrk(TokamakFastBarelyTrapped.guiding_center_4d_ode(), Δt=1.)
test_guiding_center_4d_glrk(TokamakFastLoop.guiding_center_4d_loop_ode(nl), Δt=1.)
test_guiding_center_4d_glrk(TokamakFastSurface.guiding_center_4d_surface_ode(nx, ny), Δt=1.)

test_guiding_center_4d_glrk(SymmetricLoop.guiding_center_4d_loop_ode(nl), Δt=1.)
test_guiding_center_4d_glrk(SymmetricSurface.guiding_center_4d_surface_ode(nx, ny), Δt=1.)
test_guiding_center_4d_glrk(UniformLoop.guiding_center_4d_loop_ode(nl), Δt=1.)

test_guiding_center_4d_vpglrk(TokamakSlowDeeplyPassing.guiding_center_4d_iode(), Δt=400.)
test_guiding_center_4d_vpglrk(TokamakSlowDeeplyTrapped.guiding_center_4d_iode(), Δt=400.)
test_guiding_center_4d_vpglrk(TokamakSlowBarelyPassing.guiding_center_4d_iode(), Δt=400.)
test_guiding_center_4d_vpglrk(TokamakSlowBarelyTrapped.guiding_center_4d_iode(), Δt=400.)
test_guiding_center_4d_vpglrk(TokamakSlowLoop.guiding_center_4d_loop_iode(nl), Δt=400.)
test_guiding_center_4d_vpglrk(TokamakSlowSurface.guiding_center_4d_surface_iode(nx, ny), Δt=400.)

test_guiding_center_4d_vpglrk(TokamakFastDeeplyPassing.guiding_center_4d_iode(), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakFastDeeplyTrapped.guiding_center_4d_iode(), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakFastBarelyPassing.guiding_center_4d_iode(), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakFastBarelyTrapped.guiding_center_4d_iode(), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakFastLoop.guiding_center_4d_loop_iode(nl), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakFastSurface.guiding_center_4d_surface_iode(nx, ny), Δt=1.)

test_guiding_center_4d_vpglrk(SymmetricLoop.guiding_center_4d_loop_iode(nl), Δt=1.)
test_guiding_center_4d_vpglrk(SymmetricSurface.guiding_center_4d_surface_iode(nx, ny), Δt=1.)
test_guiding_center_4d_vpglrk(UniformLoop.guiding_center_4d_loop_iode(nl), Δt=1.)
