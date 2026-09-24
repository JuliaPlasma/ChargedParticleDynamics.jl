#
# Poincaré integral invariants of the guiding centre, on the loops and surfaces of the seven
# fixture equilibria.
#
# Every bound below is ten times the error `scripts/study_poincare_invariants.jl` measures at the
# same settings, rounded up, so a regression of a decade fails. The spans are about a hundred steps
# and lie inside the horizon that script measures for each invariant: past it the quadrature on the
# sheared loop or surface fails, not the method, and a larger sample count only postpones that.
#
# The small tokamaks step at Δt = 50 rather than their default 500, at which the 4D cartesian chart
# loses the invariant to a relative 4E-2. The 3D `TokamakMediumCartesian` stops at t = 2: from
# t ≈ 5 on, members of its loop and surface leave the region where its constraint pair is regular.
#
# The theta pinch loop has `I₁ = 0` identically — it lies in `y = 0` and the flow keeps x, y and u
# fixed — so its bound is absolute, and it is the one fixture that cannot tell a method that loses
# the invariant from one that keeps it.
#

using SafeTestsets

module PoincareInvariantsTests

using GeometricIntegrators
using PoincareInvariants
using Test

export invariant, invariant_error
export NLOOP, NSURFACE, SMALL, TEST_SETTINGS

const OPTIONS = (f_abstol = 1E-12, max_iterations = 50, warn_iterations = 50)

const NLOOP = 200
const NSURFACE = 351

const SMALL = (timestep = 50.0, timespan = (0.0, 5E3))

const TEST_SETTINGS = (
    SymmetricField = (G4 = (timespan = (0.0, 1E2),), G3 = (timespan = (0.0, 1E2),)),
    ThetaPinchField = (G4 = (timespan = (0.0, 1E2),), G3 = (timespan = (0.0, 1E2),)),
    TokamakMediumCartesian = (G4 = (timespan = (0.0, 1E2),), G3 = (timespan = (0.0, 2.0),)),
    TokamakMediumCylindrical = (
        G4 = (timespan = (0.0, 1E2),), G3 = (timespan = (0.0, 1E1),)),
    TokamakSmallCartesian = (G4 = SMALL, G3 = (timespan = (0.0, 1E1),)),
    TokamakSmallCylindrical = (G4 = SMALL, G3 = SMALL),
    TokamakSmallToroidal = (G4 = SMALL, G3 = SMALL))

"The invariant at every time step of the ensemble that `ensemble(prob, pinv)` builds."
function invariant(pinv, prob, ensemble, method; kwargs...)
    compute!(pinv, integrate(ensemble(prob, pinv), method; OPTIONS..., kwargs...),
        parameters(prob))
end

"""
The largest departure of the invariant from its initial value, relative to that value, or
absolute where the invariant vanishes.
"""
function invariant_error(args...; kwargs...)
    I = invariant(args...; kwargs...)
    maximum(abs, I .- I[begin]) / (iszero(I[begin]) ? one(eltype(I)) : abs(I[begin]))
end

end

@safetestset "Poincaré invariants: the 3D and 4D guiding centre agree on every loop and surface                   " begin
    using ChargedParticleDynamics: GuidingCenter3d, GuidingCenter4d
    using GeometricIntegrators: parameters
    using PoincareInvariants: compute!
    using Test
    using ..PoincareInvariantsTests

    # The two families reach the invariant by different routes — the 4D one integrates the
    # one-form ϑ(x, u) along the loop, the 3D one lifts every point to (q, p) and integrates the
    # canonical p·dq — so this is where a wrong lift, a wrong ϑ or a wrong sign shows.
    at0(pinv, prob, ensemble) = compute!(pinv,
        [[length(p.ics.q) < 4 ? [p.ics.q; p.ics.p] : p.ics.q] for p in ensemble(prob, pinv)],
        [0.0], parameters(prob))[1]

    for name in keys(TEST_SETTINGS)
        M₃ = getfield(GuidingCenter3d, name)
        M₄ = getfield(GuidingCenter4d, name)

        I₃ = at0(M₃.poincare_invariant_1st(NLOOP), M₃.loop_hodeproblem(), M₃.loop_ensemble)
        I₄ = at0(M₄.poincare_invariant_1st(NLOOP), M₄.loop_odeproblem(), M₄.loop_ensemble)
        @test I₃≈I₄ rtol=1E-12 atol=1E-15

        isdefined(M₄, :surface_odeproblem) || continue
        @test isdefined(M₃, :surface_hodeproblem)
        J₃ = at0(M₃.poincare_invariant_2nd(NSURFACE), M₃.surface_hodeproblem(),
            M₃.surface_ensemble)
        J₄ = at0(M₄.poincare_invariant_2nd(NSURFACE), M₄.surface_odeproblem(),
            M₄.surface_ensemble)
        @test J₃ ≈ J₄ rtol=1E-12
    end

    # The theta pinch loop vanishes identically; see the header of this file.
    M = GuidingCenter4d.ThetaPinchField
    @test at0(M.poincare_invariant_1st(NLOOP), M.loop_odeproblem(), M.loop_ensemble) == 0
end

@safetestset "Poincaré invariants of the 4D guiding centre are conserved                                          " begin
    using ChargedParticleDynamics: GuidingCenter4d
    using GeometricIntegrators: Gauss, MidpointExtrapolation, SymmetricProjection, VPRKGauss
    using Test
    using ..PoincareInvariantsTests

    # (odeproblem with Gauss(2), iodeproblem with the projected VPRKGauss(2)), each over loop and
    # surface. The variational formulation keeps the invariant far better in the cartesian charts,
    # where Gauss(2) on the noncanonical ODE loses it at the order of the method.
    bounds = (
        SymmetricField = (1E-12, 1E-12),
        ThetaPinchField = (1E-14, 1E-14),
        TokamakMediumCartesian = (5E-4, 2E-8),
        TokamakMediumCylindrical = (1E-7, 1E-11),
        TokamakSmallCartesian = (2E-5, 3E-9),
        TokamakSmallCylindrical = (2E-10, 2E-11),
        TokamakSmallToroidal = (3E-10, 1E-10))

    for name in keys(TEST_SETTINGS)
        M = getfield(GuidingCenter4d, name)
        kw = TEST_SETTINGS[name].G4
        ode, iode = bounds[name]

        # The theta pinch keeps the unprojected method and the one-entry initial guess of its
        # point test in `guiding_center_4d_tests.jl`: `p` is an exact invariant there, so the
        # projection has nothing to do, and a Hermite guess would read two identical history entries.
        vprk, vkw = name == :ThetaPinchField ?
                    (VPRKGauss(2), (initialguess = MidpointExtrapolation(5),)) :
                    (SymmetricProjection(VPRKGauss(2)), NamedTuple())

        @test invariant_error(M.poincare_invariant_1st(NLOOP), M.loop_odeproblem(; kw...),
            M.loop_ensemble, Gauss(2)) < ode
        @test invariant_error(M.poincare_invariant_1st(NLOOP), M.loop_iodeproblem(; kw...),
            M.loop_ensemble, vprk; vkw...) < iode

        isdefined(M, :surface_odeproblem) || continue
        @test invariant_error(
            M.poincare_invariant_2nd(NSURFACE), M.surface_odeproblem(; kw...),
            M.surface_ensemble, Gauss(2)) < ode
        @test invariant_error(
            M.poincare_invariant_2nd(NSURFACE), M.surface_iodeproblem(; kw...),
            M.surface_ensemble, vprk; vkw...) < iode
    end
end

@safetestset "Poincaré invariants of the 3D guiding centre are conserved                                          " begin
    using ChargedParticleDynamics: GuidingCenter3d
    using GeometricIntegrators: MidpointExtrapolation, PartitionedGauss
    using Test
    using ..PoincareInvariantsTests

    bounds = (
        SymmetricField = 1E-9,
        ThetaPinchField = 1E-14,
        TokamakMediumCartesian = 5E-10,
        TokamakMediumCylindrical = 1E-12,
        TokamakSmallCartesian = 1E-12,
        TokamakSmallCylindrical = 5E-12,
        TokamakSmallToroidal = 5E-12)

    for name in keys(TEST_SETTINGS)
        M = getfield(GuidingCenter3d, name)
        kw = TEST_SETTINGS[name].G3

        # `p` is an exact invariant of the theta pinch, as in the 4D model above
        ikw = name == :ThetaPinchField ? (initialguess = MidpointExtrapolation(5),) :
              NamedTuple()

        @test invariant_error(M.poincare_invariant_1st(NLOOP), M.loop_hodeproblem(; kw...),
            M.loop_ensemble, PartitionedGauss(2); ikw...) < bounds[name]

        isdefined(M, :surface_hodeproblem) || continue
        @test invariant_error(M.poincare_invariant_2nd(NSURFACE),
            M.surface_hodeproblem(; kw...), M.surface_ensemble, PartitionedGauss(2)) <
              bounds[name]
    end
end
