

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

import ChargedParticleDynamics.GuidingCenter4d.TokamakCylindricalSlowDeeplyPassing
import ChargedParticleDynamics.GuidingCenter4d.TokamakCylindricalSlowDeeplyTrapped
import ChargedParticleDynamics.GuidingCenter4d.TokamakCylindricalSlowBarelyPassing
import ChargedParticleDynamics.GuidingCenter4d.TokamakCylindricalSlowBarelyTrapped
import ChargedParticleDynamics.GuidingCenter4d.TokamakCylindricalSlowLoop
import ChargedParticleDynamics.GuidingCenter4d.TokamakCylindricalSlowSurface

import ChargedParticleDynamics.GuidingCenter4d.TokamakCylindricalFastDeeplyPassing
import ChargedParticleDynamics.GuidingCenter4d.TokamakCylindricalFastDeeplyTrapped
import ChargedParticleDynamics.GuidingCenter4d.TokamakCylindricalFastBarelyPassing
import ChargedParticleDynamics.GuidingCenter4d.TokamakCylindricalFastBarelyTrapped
import ChargedParticleDynamics.GuidingCenter4d.TokamakCylindricalFastLoop
import ChargedParticleDynamics.GuidingCenter4d.TokamakCylindricalFastSurface

import ChargedParticleDynamics.GuidingCenter4d.TokamakCartesianFastDeeplyPassing
import ChargedParticleDynamics.GuidingCenter4d.TokamakCartesianFastDeeplyTrapped
import ChargedParticleDynamics.GuidingCenter4d.TokamakCartesianFastBarelyPassing
import ChargedParticleDynamics.GuidingCenter4d.TokamakCartesianFastBarelyTrapped
import ChargedParticleDynamics.GuidingCenter4d.TokamakCartesianFastLoop
import ChargedParticleDynamics.GuidingCenter4d.TokamakCartesianFastSurface

import ChargedParticleDynamics.GuidingCenter4d.SymmetricLoop
import ChargedParticleDynamics.GuidingCenter4d.SymmetricSurface
import ChargedParticleDynamics.GuidingCenter4d.UniformLoop

test_guiding_center_4d_glrk(TokamakCylindricalSlowDeeplyPassing.guiding_center_4d_ode(), Δt=400.)
test_guiding_center_4d_glrk(TokamakCylindricalSlowDeeplyTrapped.guiding_center_4d_ode(), Δt=400.)
test_guiding_center_4d_glrk(TokamakCylindricalSlowBarelyPassing.guiding_center_4d_ode(), Δt=400.)
test_guiding_center_4d_glrk(TokamakCylindricalSlowBarelyTrapped.guiding_center_4d_ode(), Δt=400.)
test_guiding_center_4d_glrk(TokamakCylindricalSlowLoop.guiding_center_4d_loop_ode(nl), Δt=400.)
test_guiding_center_4d_glrk(TokamakCylindricalSlowSurface.guiding_center_4d_surface_ode(nx, ny), Δt=400.)

test_guiding_center_4d_glrk(TokamakCylindricalFastDeeplyPassing.guiding_center_4d_ode(), Δt=1.)
test_guiding_center_4d_glrk(TokamakCylindricalFastDeeplyTrapped.guiding_center_4d_ode(), Δt=1.)
test_guiding_center_4d_glrk(TokamakCylindricalFastBarelyPassing.guiding_center_4d_ode(), Δt=1.)
test_guiding_center_4d_glrk(TokamakCylindricalFastBarelyTrapped.guiding_center_4d_ode(), Δt=1.)
test_guiding_center_4d_glrk(TokamakCylindricalFastLoop.guiding_center_4d_loop_ode(nl), Δt=1.)
test_guiding_center_4d_glrk(TokamakCylindricalFastSurface.guiding_center_4d_surface_ode(nx, ny), Δt=1.)

test_guiding_center_4d_glrk(TokamakCartesianFastDeeplyPassing.guiding_center_4d_ode(), Δt=1.)
test_guiding_center_4d_glrk(TokamakCartesianFastDeeplyTrapped.guiding_center_4d_ode(), Δt=1.)
test_guiding_center_4d_glrk(TokamakCartesianFastBarelyPassing.guiding_center_4d_ode(), Δt=1.)
test_guiding_center_4d_glrk(TokamakCartesianFastBarelyTrapped.guiding_center_4d_ode(), Δt=1.)
test_guiding_center_4d_glrk(TokamakCartesianFastLoop.guiding_center_4d_loop_ode(nl), Δt=1.)
test_guiding_center_4d_glrk(TokamakCartesianFastSurface.guiding_center_4d_surface_ode(nx, ny), Δt=1.)

test_guiding_center_4d_glrk(SymmetricLoop.guiding_center_4d_loop_ode(nl), Δt=1.)
test_guiding_center_4d_glrk(SymmetricSurface.guiding_center_4d_surface_ode(nx, ny), Δt=1.)
test_guiding_center_4d_glrk(UniformLoop.guiding_center_4d_loop_ode(nl), Δt=1.)

test_guiding_center_4d_vpglrk(TokamakCylindricalSlowDeeplyPassing.guiding_center_4d_iode(), Δt=400.)
test_guiding_center_4d_vpglrk(TokamakCylindricalSlowDeeplyTrapped.guiding_center_4d_iode(), Δt=400.)
test_guiding_center_4d_vpglrk(TokamakCylindricalSlowBarelyPassing.guiding_center_4d_iode(), Δt=400.)
test_guiding_center_4d_vpglrk(TokamakCylindricalSlowBarelyTrapped.guiding_center_4d_iode(), Δt=400.)
test_guiding_center_4d_vpglrk(TokamakCylindricalSlowLoop.guiding_center_4d_loop_iode(nl), Δt=400.)
test_guiding_center_4d_vpglrk(TokamakCylindricalSlowSurface.guiding_center_4d_surface_iode(nx, ny), Δt=400.)

test_guiding_center_4d_vpglrk(TokamakCylindricalFastDeeplyPassing.guiding_center_4d_iode(), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakCylindricalFastDeeplyTrapped.guiding_center_4d_iode(), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakCylindricalFastBarelyPassing.guiding_center_4d_iode(), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakCylindricalFastBarelyTrapped.guiding_center_4d_iode(), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakCylindricalFastLoop.guiding_center_4d_loop_iode(nl), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakCylindricalFastSurface.guiding_center_4d_surface_iode(nx, ny), Δt=1.)

test_guiding_center_4d_vpglrk(TokamakCartesianFastDeeplyPassing.guiding_center_4d_iode(), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakCartesianFastDeeplyTrapped.guiding_center_4d_iode(), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakCartesianFastBarelyPassing.guiding_center_4d_iode(), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakCartesianFastBarelyTrapped.guiding_center_4d_iode(), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakCartesianFastLoop.guiding_center_4d_loop_iode(nl), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakCartesianFastSurface.guiding_center_4d_surface_iode(nx, ny), Δt=1.)

test_guiding_center_4d_vpglrk(SymmetricLoop.guiding_center_4d_loop_iode(nl), Δt=1.)
test_guiding_center_4d_vpglrk(SymmetricSurface.guiding_center_4d_surface_iode(nx, ny), Δt=1.)
test_guiding_center_4d_vpglrk(UniformLoop.guiding_center_4d_loop_iode(nl), Δt=1.)
