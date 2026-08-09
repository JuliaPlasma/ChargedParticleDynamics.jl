
using SafeTestsets

module GuidingCenter4dTests

using GeometricIntegrators
using Test

const nl = 100
const nx = 10
const ny = 10

# See `guiding_center_3d_tests.jl` for why `f_abstol` has to stay above the residual's round-off
# floor and why `f_reltol` is left at its default.
const options = (f_abstol=1E-12, max_iterations=50, warn_iterations=50)

export test_guiding_center_4d
export nl, nx, ny

# ---------------------------------------------------------------------------------------------
# Time step and time span limits of the 4D guiding centre equilibria
# ---------------------------------------------------------------------------------------------
#
# Every call runs at its module's declared thousand-step example except the three noted here.
#
# **`SolovevSymmetricField` is the constrained one.** `‖ϑ‖ ≈ 256` against ≈ 15 for the ITER Solov'ev
# cases makes it the worst-conditioned equilibrium in the package. `Δt = 1E0` is the only stable rate
# measured for `Gauss(2)`: over the same physical span the maximum relative energy error is 1.9E-3
# there against 1.8E+1 at `Δt = 1E2` and 1.8E+11 at `Δt = 1E1`. Its variational form additionally
# needs a projected integrator; see the note on `test_guiding_center_4d` below.
#
# **Larger steps are not better here, and the fixed-span evidence that they are is an artifact.**
# Holding the span fixed and growing Δt takes fewer steps, so less instability accumulates — at
# Δt = 1E5 over t ∈ [0, 1E6] the run is ten steps, which cannot resolve the orbit at all. Holding the
# *step count* fixed at 1000 gives the ordinary monotone picture: 9.4E-8 at Δt = 1E-1, 1.9E-3 at 1E0,
# 1.8E+11 at 1E1, 1.3E+32 at 1E4.
#
# **Three exceptions, each owned by the formulation that blocks:**
#
#   SolovevIterXpoint.iodeproblem      Δt = 1E-1 over t ∈ [0, 1E1]. At the declared default the
#                                      variational and `ode` orbits agree — the projection made them
#                                      consistent — but their shared energy error is O(1), and
#                                      `deeply_passing` emits 39 line-search-at-round-off messages
#                                      that no f_abstol between 1E-12 and 1E-10 removes.
#   TokamakSmallCartesian.iodeproblem  Δt = 10 over t ∈ [0, 1E3]. The variational forms do not
#   TokamakSmallCylindrical.iodeproblem_dg   converge at the declared Δt = 500: `deeply_passing`
#                                      reaches 3.2E+41 with 495 capped steps, against 6.1E-4 and
#                                      silence at Δt = 10.
#   ThetaPinchField.loop_iodeproblem   *not* a time exception — it keeps the unprojected
#                                      `VPRKGauss(2)`, because `p` is an exact invariant here so
#                                      there is no drift to project, and the projection's inner
#                                      Hermite guess would re-introduce the degenerate-history
#                                      warning the `MidpointExtrapolation` there exists to avoid.
#
# The steps above were chosen against energy rather than against solver silence, which proves less
# than it looks; see the header of `quiet_solver_warnings.jl`.

# Asserts that the integration returns a solution rather than `@test_nowarn`; see the comment on
# the corresponding functions in `guiding_center_3d_tests.jl` for why.
function test_guiding_center_4d(equ::ODEProblem; kwargs...)
    @test integrate(equ, Gauss(2); options..., kwargs...) isa GeometricSolution
end

# The variational problems are integrated with a *projected* method throughout.
#
# The degenerate variational discretisation of the guiding centre carries a parasitic mode that
# plain `VPRKGauss(2)` does not control, and it grows as the time step shrinks — the reverse of the
# usual convergence expectation. `SolovevSymmetricField` is where that is unmissable (`‖ϑ‖ ≈ 256`,
# relative energy error 6.2E+25 at `Δt = 1E0`), but it is a property of the formulation rather than
# of that equilibrium, so every variational orbit here gets the same treatment.
# `SymmetricProjection` constrains the drift off `p = ϑ(q)` that the mode feeds on. See
# `docs/src/findings.md`.
function test_guiding_center_4d(equ::Union{IODEProblem,LODEProblem}; kwargs...)
    @test integrate(equ, SymmetricProjection(VPRKGauss(2)); options..., kwargs...) isa GeometricSolution
end

# Explicit-method form, for the few calls that need something other than the defaults above.
function test_guiding_center_4d(equ, method; kwargs...)
    @test integrate(equ, method; options..., kwargs...) isa GeometricSolution
end

end



@safetestset "Guiding Centre Dynamics in 4D with ITER-like Solov'ev Equilibrium with X-Point                      " begin

    using ChargedParticleDynamics.GuidingCenter4d.SolovevIterXpoint
    using ..GuidingCenter4dTests

    test_guiding_center_4d(odeproblem(initial_conditions_trapped()))
    test_guiding_center_4d(odeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_trapped()))
    # test_guiding_center_4d(loop_odeproblem(nl), Δt=1.)
    # test_guiding_center_4d(surface_odeproblem(nx, ny), Δt=1.)

    # `trapped` drops its `Δt = 1E4` override because the module default is simply better for it:
    # the relative energy error over the run is 8.8E-5 at `Δt = 1E0` against 2.2E-1 at `Δt = 1E4`,
    # for both formulations.
    #
    # The other four keep `Δt = 1E-1`. At the module default the variational orbits agree with the
    # `ode` ones to three digits in energy error (5.804 against 5.808 for `deeply_passing`) — that
    # is, the projection has made the two formulations consistent — but the shared error is O(1),
    # and `deeply_passing` emits 39 line-search-at-round-off messages there which no `f_abstol`
    # between 1E-12 and 1E-10 removes. `Δt = 1E-1` is where this equilibrium is actually resolved.
    test_guiding_center_4d(iodeproblem(initial_conditions_trapped()))
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_passing(); timestep=1E-1, timespan=(0, 1E1)))
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_trapped(); timestep=1E-1, timespan=(0, 1E1)))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_passing(); timestep=1E-1, timespan=(0, 1E1)))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_trapped(); timestep=1E-1, timespan=(0, 1E1)))
    # test_guiding_center_4d(loop_iodeproblem(nl), Δt=1.)
    # test_guiding_center_4d(surface_iodeproblem(nx, ny), Δt=1.)

end


@safetestset "Guiding Centre Dynamics in 4D with medium-size Tokamak Equilibrium in Cartesian Coordinates         " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCartesian
    using ..GuidingCenter4dTests

    test_guiding_center_4d(odeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_trapped()))
    # test_guiding_center_4d(loop_odeproblem(nl))
    # test_guiding_center_4d(surface_odeproblem(nx, ny))

    # At the module's own step the projected variational orbits are silent and track the `ode` ones
    # to within a factor of a few in energy error (4.7E-6 against 3.1E-6 for `barely_passing`).
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_trapped()))
    # test_guiding_center_4d(loop_iodeproblem(nl))
    # test_guiding_center_4d(surface_iodeproblem(nx, ny))

end


@safetestset "Guiding Centre Dynamics in 4D with medium-size Tokamak Equilibrium in Cylindrical Coordinates       " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCylindrical
    using ..GuidingCenter4dTests

    test_guiding_center_4d(odeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_trapped()))
    # test_guiding_center_4d(loop_odeproblem(nl))
    # test_guiding_center_4d(surface_odeproblem(nx, ny))

    test_guiding_center_4d(iodeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_trapped()))
    # test_guiding_center_4d(loop_iodeproblem(nl))
    # test_guiding_center_4d(surface_iodeproblem(nx, ny))

end


@safetestset "Guiding Centre Dynamics in 4D with small-size Tokamak Equilibrium in Cartesian Coordinates          " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCartesian
    using ..GuidingCenter4dTests

    test_guiding_center_4d(odeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_trapped()))
    # test_guiding_center_4d(loop_odeproblem(nl))
    # test_guiding_center_4d(surface_odeproblem(nx, ny))

    # The variational integrators do not converge at this module's default time step of Δt = 500,
    # and the projection does not rescue them: `deeply_passing` reaches a relative energy error of
    # 3.2E+41 with |R| = 2.3E+18 and 495 solver warnings there, against 6.1E-4 and silence at
    # Δt = 10. This is the one block in the file whose override is load-bearing.
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_passing(); timestep=10.0, timespan=(0, 1E3)))
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_trapped(); timestep=10.0, timespan=(0, 1E3)))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_passing(); timestep=10.0, timespan=(0, 1E3)))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_trapped(); timestep=10.0, timespan=(0, 1E3)))
    # test_guiding_center_4d(loop_iodeproblem(nl; timestep = 10., timespan = [0, 1E3]))
    # test_guiding_center_4d(surface_iodeproblem(nx, ny; timestep = 10., timespan = [0, 1E3]))

end


@safetestset "Guiding Centre Dynamics in 4D with small-size Tokamak Equilibrium in Cylindrical Coordinates        " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCylindrical
    using ..GuidingCenter4dTests

    test_guiding_center_4d(odeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_trapped()))
    # test_guiding_center_4d(loop_odeproblem(nl))
    # test_guiding_center_4d(surface_odeproblem(nx, ny))

    test_guiding_center_4d(iodeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_trapped()))
    # test_guiding_center_4d(loop_iodeproblem(nl; timestep = 10., timespan = [0, 1E3]))
    # test_guiding_center_4d(surface_iodeproblem(nx, ny; timestep = 10., timespan = [0, 1E3]))

    # The κ-dependent "dg" formulation. Its `ḡ` is checked against finite differences in
    # `structure_tests.jl`, but until now it had never been integrated — and `κ = 0`, which is the
    # default, reduces it to the plain `iodeproblem` and so would not exercise the κ terms at all.
    test_guiding_center_4d(iodeproblem_dg(initial_conditions_barely_passing(); κ=0.0, timespan=(0.0, 1E3), timestep=10.0))
    test_guiding_center_4d(iodeproblem_dg(initial_conditions_barely_passing(); κ=0.5, timespan=(0.0, 1E3), timestep=10.0))
    test_guiding_center_4d(iodeproblem_dg(initial_conditions_barely_passing(); κ=1.0, timespan=(0.0, 1E3), timestep=10.0))

end


@safetestset "Guiding Centre Dynamics in 4D with small-size Tokamak Equilibrium in Toroidal Coordinates           " begin

    using ChargedParticleDynamics.GuidingCenter4d.TokamakSmallToroidal
    using ..GuidingCenter4dTests

    test_guiding_center_4d(odeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_trapped()))
    # test_guiding_center_4d(loop_odeproblem(nl))
    # test_guiding_center_4d(surface_odeproblem(nx, ny))

    test_guiding_center_4d(iodeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_trapped()))
    # test_guiding_center_4d(loop_iodeproblem(nl))
    # test_guiding_center_4d(surface_iodeproblem(nx, ny))

end


@safetestset "Guiding Centre Dynamics in 4D with symmetric Solov'ev Equilibrium                                   " begin

    using ChargedParticleDynamics.GuidingCenter4d.SolovevSymmetricField
    using ..GuidingCenter4dTests

    # `Δt = 1E0` is not merely a convenient rate for `Gauss(2)`, it is the only stable one measured:
    # over the same physical span the maximum relative energy error is 1.9E-3 there against 1.8E+11
    # at `Δt = 1E1` and 1.8E+1 at `Δt = 1E2`. All four orbits converge in 2.1-2.3 Newton iterations
    # and reach the iteration cap on no step; at `Δt = 1E2`, `1E1`, `5E2` and `1E3` they cap on 45 %,
    # 4.6 %, 29 % and 14 % of steps respectively. A thousand small steps is close to free next to
    # that, since a step converging in two iterations is far cheaper than one grinding to fifty.
    test_guiding_center_4d(odeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(odeproblem(initial_conditions_deeply_trapped()))

    # The variational problem needs a *projected* integrator here, which is the one place in the
    # package where plain `VPRKGauss(2)` is not enough.
    #
    # `‖ϑ‖ ≈ 256` makes this equilibrium's degenerate variational discretisation carry a parasitic
    # mode that plain VPRK does not control, and it grows as the step shrinks: over `t ∈ [0, 1E4]`
    # the maximum relative energy error is 2.0E-1 at `Δt = 5E3` but 2.3E+3 at `Δt = 1E2`, 1.5E+18 at
    # `Δt = 1E1` and 6.2E+25 at `Δt = 1E0`, with `|R|` reaching 1.3E+5 against a machine radius of
    # 2.5. The solver then reaches the iteration cap on 97 % of steps — but that is the *symptom*.
    # It is being asked to solve for a state that has already left the device, and no solver setting
    # repairs it: `DogLeg` throws on the NaN direction for all four orbits, `StrongWolfe` and
    # `Static` throw a `SingularException` or cap just as hard, and `PartitionedGauss(2)` is singular.
    #
    # `SymmetricProjection` is what fixes it, because it constrains exactly the drift off `p = ϑ(q)`
    # that the parasitic mode feeds on. All four orbits then integrate silently and stay bounded at
    # `|R| ≤ 6.5`, the same excursion `Gauss(2)` gives. `MidpointProjection` is not sufficient: it
    # holds `barely_trapped` but throws on the NaN direction for the other three.
    #
    # Note what this does and does not buy. It restores *stability*, not accuracy: the projected
    # energy error is O(1) against the `odeproblem`'s 1.9E-3 over the same span. These calls assert
    # that the integration completes, and on this equilibrium the variational formulation should not
    # be read as a quantitative result. See `docs/src/findings.md`.
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_passing()))
    test_guiding_center_4d(iodeproblem(initial_conditions_barely_trapped()))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_passing()))
    test_guiding_center_4d(iodeproblem(initial_conditions_deeply_trapped()))

end


# These two equilibria carry only a Poincaré loop and surface parameterisation, no point initial
# condition, so what there is to test is the loop and surface problems. Both blocks used to be
# entirely commented out and ran zero tests: their calls passed the sampling counts positionally —
# `loop_odeproblem(nl)` — which the `PoincareInvariants` 0.5 rewrite replaced with keyword-only
# constructors that take the sample count at the invariant rather than at the problem.
@safetestset "Guiding Centre Dynamics in 4D with symmetric Equilibrium                                            " begin

    using ChargedParticleDynamics.GuidingCenter4d.SymmetricField
    using ..GuidingCenter4dTests

    test_guiding_center_4d(loop_odeproblem())
    test_guiding_center_4d(surface_odeproblem())

    test_guiding_center_4d(loop_iodeproblem())
    test_guiding_center_4d(surface_iodeproblem())

end


@safetestset "Guiding Centre Dynamics in 4D with Theta Pinch Equilibrium                                          " begin

    using GeometricIntegrators: MidpointExtrapolation, VPRKGauss
    using ChargedParticleDynamics.GuidingCenter4d.ThetaPinchField
    using ..GuidingCenter4dTests

    # The theta pinch defines a loop but no surface.
    test_guiding_center_4d(loop_odeproblem())

    # `B = B₀ e_z` is uniform, so ∇B = 0 and the force vanishes identically: the guiding centre
    # streams along z at constant u and never leaves its field line, which leaves x, y and u
    # bit-for-bit constant and makes p = A(x,y) + u b an exact invariant of every orbit in this
    # equilibrium — no choice of loop avoids it. The default Hermite extrapolation of the initial
    # guess then finds two identical `p` history entries on every step, and warns while falling back
    # to the constant that is in fact the exact answer. Same degeneracy and same remedy as the Pauli
    # theta pinch; see `pauli_particle_3d_tests.jl`. `MidpointExtrapolation` reads one history entry
    # rather than two, so the trajectory is bitwise unchanged. The `odeproblem` above needs none of
    # this: its `q` does move, so its Hermite guess is well posed.
    #
    # This is also the one variational call in the file that keeps the *unprojected* `VPRKGauss(2)`.
    # There is nothing here for `SymmetricProjection` to do — `p` is an exact invariant, so the
    # parasitic drift it exists to control is identically zero — and its inner projection solve
    # brings its own Hermite guess, which re-introduces the very degeneracy the
    # `MidpointExtrapolation` above is here to avoid: two identical `p` history entries, and two
    # warnings that the plain method does not produce.
    test_guiding_center_4d(loop_iodeproblem(), VPRKGauss(2); initialguess=MidpointExtrapolation(5))

end
