
using Test

module GuidingCenter4dTests

    using GeometricIntegrators

    const Δt = 1.
    const nt = 1

    export test_guiding_center_4d_glrk, test_guiding_center_4d_vpglrk

    function test_guiding_center_4d_glrk(ode; Δt = Δt)
        int = Integrator(ode, getTableauGLRK(1), Δt)
        sol = integrate(ode, int, nt)
    end

    function test_guiding_center_4d_vpglrk(ode; Δt = Δt)
        int = Integrator(ode, getTableauVPGLRK(1), Δt)
        sol = integrate(ode, int, nt)
    end

end


nl = 100
nx = 10
ny = 10


using .GuidingCenter4dTests

import ChargedParticleDynamics.GuidingCenter4d.GuidingCenter4dSolovevQuadratic
import ChargedParticleDynamics.GuidingCenter4d.GuidingCenter4dSymmetricQuadratic
import ChargedParticleDynamics.GuidingCenter4d.GuidingCenter4dThetaPinch

import ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCartesian
import ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCylindrical
import ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCylindrical


test_guiding_center_4d_glrk(TokamakMediumCartesian.guiding_center_4d_ode(TokamakMediumCartesian.initial_conditions_barely_passing()...), Δt=1.)
test_guiding_center_4d_glrk(TokamakMediumCartesian.guiding_center_4d_ode(TokamakMediumCartesian.initial_conditions_barely_trapped()...), Δt=1.)
test_guiding_center_4d_glrk(TokamakMediumCartesian.guiding_center_4d_ode(TokamakMediumCartesian.initial_conditions_deeply_passing()...), Δt=1.)
test_guiding_center_4d_glrk(TokamakMediumCartesian.guiding_center_4d_ode(TokamakMediumCartesian.initial_conditions_deeply_trapped()...), Δt=1.)
test_guiding_center_4d_glrk(TokamakMediumCartesian.guiding_center_4d_loop_ode(nl), Δt=1.)
test_guiding_center_4d_glrk(TokamakMediumCartesian.guiding_center_4d_surface_ode(nx, ny), Δt=1.)

test_guiding_center_4d_glrk(TokamakMediumCylindrical.guiding_center_4d_ode(TokamakMediumCylindrical.initial_conditions_barely_passing()...), Δt=1.)
test_guiding_center_4d_glrk(TokamakMediumCylindrical.guiding_center_4d_ode(TokamakMediumCylindrical.initial_conditions_barely_trapped()...), Δt=1.)
test_guiding_center_4d_glrk(TokamakMediumCylindrical.guiding_center_4d_ode(TokamakMediumCylindrical.initial_conditions_deeply_passing()...), Δt=1.)
test_guiding_center_4d_glrk(TokamakMediumCylindrical.guiding_center_4d_ode(TokamakMediumCylindrical.initial_conditions_deeply_trapped()...), Δt=1.)
test_guiding_center_4d_glrk(TokamakMediumCylindrical.guiding_center_4d_loop_ode(nl), Δt=1.)
test_guiding_center_4d_glrk(TokamakMediumCylindrical.guiding_center_4d_surface_ode(nx, ny), Δt=1.)

test_guiding_center_4d_glrk(TokamakSmallCylindrical.guiding_center_4d_ode(TokamakSmallCylindrical.initial_conditions_barely_passing()...), Δt=400.)
test_guiding_center_4d_glrk(TokamakSmallCylindrical.guiding_center_4d_ode(TokamakSmallCylindrical.initial_conditions_barely_trapped()...), Δt=400.)
test_guiding_center_4d_glrk(TokamakSmallCylindrical.guiding_center_4d_ode(TokamakSmallCylindrical.initial_conditions_deeply_passing()...), Δt=400.)
test_guiding_center_4d_glrk(TokamakSmallCylindrical.guiding_center_4d_ode(TokamakSmallCylindrical.initial_conditions_deeply_trapped()...), Δt=400.)
test_guiding_center_4d_glrk(TokamakSmallCylindrical.guiding_center_4d_loop_ode(nl), Δt=400.)
test_guiding_center_4d_glrk(TokamakSmallCylindrical.guiding_center_4d_surface_ode(nx, ny), Δt=400.)

test_guiding_center_4d_glrk(GuidingCenter4dSolovevQuadratic.guiding_center_4d_ode(GuidingCenter4dSolovevQuadratic.initial_conditions_barely_passing()...), Δt=1.)
test_guiding_center_4d_glrk(GuidingCenter4dSolovevQuadratic.guiding_center_4d_ode(GuidingCenter4dSolovevQuadratic.initial_conditions_barely_trapped()...), Δt=1.)
test_guiding_center_4d_glrk(GuidingCenter4dSolovevQuadratic.guiding_center_4d_ode(GuidingCenter4dSolovevQuadratic.initial_conditions_deeply_passing()...), Δt=1.)
test_guiding_center_4d_glrk(GuidingCenter4dSolovevQuadratic.guiding_center_4d_ode(GuidingCenter4dSolovevQuadratic.initial_conditions_deeply_trapped()...), Δt=1.)

test_guiding_center_4d_glrk(GuidingCenter4dSymmetricQuadratic.guiding_center_4d_loop_ode(nl), Δt=1.)
test_guiding_center_4d_glrk(GuidingCenter4dSymmetricQuadratic.guiding_center_4d_surface_ode(nx, ny), Δt=1.)

test_guiding_center_4d_glrk(GuidingCenter4dThetaPinch.guiding_center_4d_loop_ode(nl), Δt=1.)

test_guiding_center_4d_vpglrk(TokamakMediumCartesian.guiding_center_4d_iode(TokamakMediumCartesian.initial_conditions_barely_passing()...), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakMediumCartesian.guiding_center_4d_iode(TokamakMediumCartesian.initial_conditions_barely_trapped()...), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakMediumCartesian.guiding_center_4d_iode(TokamakMediumCartesian.initial_conditions_deeply_passing()...), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakMediumCartesian.guiding_center_4d_iode(TokamakMediumCartesian.initial_conditions_deeply_trapped()...), Δt=1.)
# test_guiding_center_4d_vpglrk(TokamakMediumCartesian.guiding_center_4d_loop_iode(nl), Δt=1.)
# test_guiding_center_4d_vpglrk(TokamakMediumCartesian.guiding_center_4d_surface_iode(nx, ny), Δt=1.)

test_guiding_center_4d_vpglrk(TokamakMediumCylindrical.guiding_center_4d_iode(TokamakMediumCylindrical.initial_conditions_barely_passing()...), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakMediumCylindrical.guiding_center_4d_iode(TokamakMediumCylindrical.initial_conditions_barely_trapped()...), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakMediumCylindrical.guiding_center_4d_iode(TokamakMediumCylindrical.initial_conditions_deeply_passing()...), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakMediumCylindrical.guiding_center_4d_iode(TokamakMediumCylindrical.initial_conditions_deeply_trapped()...), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakMediumCylindrical.guiding_center_4d_loop_iode(nl), Δt=1.)
test_guiding_center_4d_vpglrk(TokamakMediumCylindrical.guiding_center_4d_surface_iode(nx, ny), Δt=1.)

test_guiding_center_4d_vpglrk(TokamakSmallCylindrical.guiding_center_4d_iode(TokamakSmallCylindrical.initial_conditions_barely_passing()...), Δt=400.)
test_guiding_center_4d_vpglrk(TokamakSmallCylindrical.guiding_center_4d_iode(TokamakSmallCylindrical.initial_conditions_barely_trapped()...), Δt=400.)
test_guiding_center_4d_vpglrk(TokamakSmallCylindrical.guiding_center_4d_iode(TokamakSmallCylindrical.initial_conditions_deeply_passing()...), Δt=400.)
test_guiding_center_4d_vpglrk(TokamakSmallCylindrical.guiding_center_4d_iode(TokamakSmallCylindrical.initial_conditions_deeply_trapped()...), Δt=400.)
test_guiding_center_4d_vpglrk(TokamakSmallCylindrical.guiding_center_4d_loop_iode(nl), Δt=400.)
test_guiding_center_4d_vpglrk(TokamakSmallCylindrical.guiding_center_4d_surface_iode(nx, ny), Δt=400.)

test_guiding_center_4d_vpglrk(GuidingCenter4dSolovevQuadratic.guiding_center_4d_iode(GuidingCenter4dSolovevQuadratic.initial_conditions_barely_passing()...), Δt=1.)
test_guiding_center_4d_vpglrk(GuidingCenter4dSolovevQuadratic.guiding_center_4d_iode(GuidingCenter4dSolovevQuadratic.initial_conditions_barely_trapped()...), Δt=1.)
test_guiding_center_4d_vpglrk(GuidingCenter4dSolovevQuadratic.guiding_center_4d_iode(GuidingCenter4dSolovevQuadratic.initial_conditions_deeply_passing()...), Δt=1.)
test_guiding_center_4d_vpglrk(GuidingCenter4dSolovevQuadratic.guiding_center_4d_iode(GuidingCenter4dSolovevQuadratic.initial_conditions_deeply_trapped()...), Δt=1.)

test_guiding_center_4d_vpglrk(GuidingCenter4dSymmetricQuadratic.guiding_center_4d_loop_iode(nl), Δt=1.)
test_guiding_center_4d_vpglrk(GuidingCenter4dSymmetricQuadratic.guiding_center_4d_surface_iode(nx, ny), Δt=1.)

test_guiding_center_4d_vpglrk(GuidingCenter4dThetaPinch.guiding_center_4d_loop_iode(nl), Δt=1.)
