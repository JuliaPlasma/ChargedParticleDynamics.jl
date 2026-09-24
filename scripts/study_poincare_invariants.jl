#
# Over what time span can a Poincaré invariant of each equilibrium be computed at all?
#
# The flow shears the advected loop or surface, and the quadrature that evaluates the invariant on
# it — a Fourier plan on the loop, Chebyshev on Padua points on the surface — degrades as the
# neighbouring sample points drift apart. Past some time the computed invariant then departs from
# its initial value exponentially, whatever the integrator, and this reads exactly like the method
# losing the invariant. More points postpone that time; they do not remove it.
#
# So the span has to be measured, against a refined reference rather than a fixed tolerance: for
# every loop and surface this integrates the same ensemble at a coarse and a fine sample count and
# reports the first time at which the two invariants disagree by more than `AGREEMENT`, together
# with the conservation error of the fine one before that time. The test suite evaluates the
# invariants over a span inside it. The quadrature limit does not depend on the method, so each
# family is integrated with one: `Gauss(2)` on the 4D `ODEProblem`, `PartitionedGauss(2)` on the 3D
# `HODEProblem`.
#
# Run with:  julia --project=test scripts/study_poincare_invariants.jl [horizons] [tests] 2>/dev/null
#
# `horizons` is section 1 and `tests` section 2; with neither, both run. Name fixtures, for example
# `SymmetricField`, to run only those.
# The diverging members fill stderr with solver warnings; the table goes to stdout.
#

using ChargedParticleDynamics
using GeometricIntegrators
using GeometricIntegrators: MidpointExtrapolation
using PoincareInvariants

const G3 = ChargedParticleDynamics.GuidingCenter3d
const G4 = ChargedParticleDynamics.GuidingCenter4d

const OPTIONS = (f_abstol = 1E-12, max_iterations = 50, warn_iterations = 50)

# the relative disagreement between the coarse and the fine invariant that ends the usable span
const AGREEMENT = 1E-10

# coarse and fine sample counts; the surface counts are Padua numbers, 351 = 26·27/2, 861 = 41·42/2
const NLOOP = (200, 800)
const NSURFACE = (351, 861)

# Integrate each member of the ensemble on its own, as `compute!` takes one trajectory per sample
# point, and return the invariant at every time step. The 3D guiding centre is canonical, and its
# invariants are taken on the stacked `(q, p)`; the 4D one reads `q` alone.
function invariant(pinv, ensemble, method; kwargs...)
    sols = [integrate(prob, method; OPTIONS..., kwargs...) for prob in ensemble]
    nt = ntime(sols[begin])
    ts = [sols[begin].t[n] for n in 0:nt]
    point(s, n) = getdim(pinv) > length(s.q[n]) ? [s.q[n]; s.p[n]] : s.q[n]
    I = compute!(pinv, [[point(s, n) for n in 0:nt] for s in sols], ts,
        parameters(first(ensemble)))
    ts, I, count(diverged, sols)
end

# A member that leaves the device fails the quadrature at both sample counts at once, which the
# agreement test below cannot see, so these are counted separately.
function diverged(s)
    q₀ = maximum(abs, s.q[0])
    any(n -> !all(isfinite, s.q[n]) || maximum(abs, s.q[n]) > 1E3 * (1 + q₀), 0:ntime(s))
end

function horizon(name, kind, prob, pinvs, ensemble, method)
    ts, I₋, _ = invariant(pinvs[1], ensemble(prob, pinvs[1]), method)
    _, I₊, nd = invariant(pinvs[2], ensemble(prob, pinvs[2]), method)
    # The theta pinch loop has `I₁ = 0` identically, and its errors are absolute.
    scale = iszero(I₊[begin]) ? one(eltype(I₊)) : abs(I₊[begin])
    n = findfirst(i -> abs(I₊[i] - I₋[i]) > AGREEMENT * scale, eachindex(ts))
    last = n === nothing ? lastindex(ts) : n - 1
    err = maximum(abs, I₊[begin:last] .- I₊[begin]) / scale
    println(rpad(name, 26), rpad(kind, 9), "Δt = ", rpad(timestep(prob), 8),
        "span ", rpad(ts[last], 10), n === nothing ? "(whole run)" : "           ",
        iszero(I₊[begin]) ? "  abs." : "  rel.", " error before it ", round(err; sigdigits = 2),
        "   I(0) = ", I₊[begin],
        nd > 0 ? "   $nd of $(length(pinvs[2].points[:, 1])) members diverged" : "")
    flush(stdout)
end

const ALL_FIXTURES = (:SymmetricField, :ThetaPinchField, :TokamakMediumCartesian,
    :TokamakMediumCylindrical, :TokamakSmallCartesian, :TokamakSmallCylindrical,
    :TokamakSmallToroidal)

# the fixtures named on the command line, or all of them
const FIXTURES = let named = filter(n -> string(n) in ARGS, ALL_FIXTURES)
    isempty(named) ? ALL_FIXTURES : named
end

function guiding_centre_4d()
    println("4D guiding centre — ODEProblem, Gauss(2), over each module's default time span\n")
    for name in FIXTURES
        M = getfield(G4, name)
        horizon(name, "loop", M.loop_odeproblem(),
            [M.poincare_invariant_1st(n) for n in NLOOP], M.loop_ensemble, Gauss(2))
        isdefined(M, :surface_odeproblem) || continue
        horizon(name, "surface", M.surface_odeproblem(),
            [M.poincare_invariant_2nd(n) for n in NSURFACE], M.surface_ensemble, Gauss(2))
    end
end

function guiding_centre_3d()
    println("\n3D guiding centre — HODEProblem, PartitionedGauss(2), over each module's default " *
            "time span\n")
    for name in FIXTURES
        M = getfield(G3, name)
        horizon(name, "loop", M.loop_hodeproblem(),
            [M.poincare_invariant_1st(n) for n in NLOOP], M.loop_ensemble, PartitionedGauss(2))
        isdefined(M, :surface_hodeproblem) || continue
        horizon(name, "surface", M.surface_hodeproblem(),
            [M.poincare_invariant_2nd(n) for n in NSURFACE], M.surface_ensemble,
            PartitionedGauss(2))
    end
end

# Section 2: the settings the test suite runs, at the coarse sample counts it uses. About a hundred
# steps each, inside every horizon above. The small tokamaks step at Δt = 50: their default 500
# leaves the cartesian chart at a relative error of 4E-2, which is the method rather than the
# quadrature — the two sample counts still agree there. The 3D `TokamakMediumCartesian` stops at
# t = 2, before its first member leaves the region where its constraint pair is regular.
const SMALL = (timestep = 50.0, timespan = (0.0, 5E3))

const TEST_SETTINGS = (
    SymmetricField = (G4 = (timespan = (0.0, 1E2),), G3 = (timespan = (0.0, 1E2),)),
    ThetaPinchField = (G4 = (timespan = (0.0, 1E2),), G3 = (timespan = (0.0, 1E2),)),
    TokamakMediumCartesian = (G4 = (timespan = (0.0, 1E2),),
        G3 = (timespan = (0.0, 2.0),)),
    TokamakMediumCylindrical = (
        G4 = (timespan = (0.0, 1E2),), G3 = (timespan = (0.0, 1E1),)),
    TokamakSmallCartesian = (G4 = SMALL, G3 = (timespan = (0.0, 1E1),)),
    TokamakSmallCylindrical = (G4 = SMALL, G3 = SMALL),
    TokamakSmallToroidal = (G4 = SMALL, G3 = SMALL))

function conservation(name, kind, formulation, pinv, ensemble, method; kwargs...)
    t = @elapsed ts, I, nd = invariant(pinv, ensemble, method; kwargs...)
    scale = iszero(I[begin]) ? one(eltype(I)) : abs(I[begin])
    err = maximum(abs, I .- I[begin]) / scale
    println(rpad(name, 26), rpad(kind, 9), rpad(formulation, 24), "t ≤ ", rpad(ts[end], 8),
        iszero(I[begin]) ? "abs. " : "rel. ", rpad(round(err; sigdigits = 2), 10),
        rpad("$(round(t; digits = 1)) s", 9), nd > 0 ? "$nd diverged" : "")
    flush(stdout)
end

function test_settings()
    println("\nThe test suite's settings, at $(NLOOP[1]) loop and $(NSURFACE[1]) surface points\n")
    for name in FIXTURES
        M = getfield(G4, name)
        kw = TEST_SETTINGS[name].G4
        # The theta pinch keeps the unprojected variational method and a one-entry initial guess,
        # as in `test/poincare_invariants_tests.jl`: `p` is an exact invariant there.
        θ = name == :ThetaPinchField
        vprk = θ ? VPRKGauss(2) : SymmetricProjection(VPRKGauss(2))
        ikw = θ ? (initialguess = MidpointExtrapolation(5),) : NamedTuple()
        for (formulation, loop, surface, method, mkw) in (
            (
            "4D odeproblem", :loop_odeproblem, :surface_odeproblem, Gauss(2), NamedTuple()),
            ("4D iodeproblem", :loop_iodeproblem, :surface_iodeproblem, vprk, ikw))
            p1 = M.poincare_invariant_1st(NLOOP[1])
            conservation(name, "loop", formulation, p1,
                M.loop_ensemble(getfield(M, loop)(; kw...), p1), method; mkw...)
            isdefined(M, surface) || continue
            p2 = M.poincare_invariant_2nd(NSURFACE[1])
            conservation(name, "surface", formulation, p2,
                M.surface_ensemble(getfield(M, surface)(; kw...), p2), method; mkw...)
        end

        M = getfield(G3, name)
        kw = TEST_SETTINGS[name].G3
        p1 = M.poincare_invariant_1st(NLOOP[1])
        conservation(name, "loop", "3D hodeproblem", p1,
            M.loop_ensemble(M.loop_hodeproblem(; kw...), p1), PartitionedGauss(2); ikw...)
        isdefined(M, :surface_hodeproblem) || continue
        p2 = M.poincare_invariant_2nd(NSURFACE[1])
        conservation(name, "surface", "3D hodeproblem", p2,
            M.surface_ensemble(M.surface_hodeproblem(; kw...), p2), PartitionedGauss(2))
    end
end

const SECTIONS = filter(in(("horizons", "tests")), ARGS)

if isempty(SECTIONS) || "horizons" in SECTIONS
    guiding_centre_4d()
    guiding_centre_3d()
end
(isempty(SECTIONS) || "tests" in SECTIONS) && test_settings()
