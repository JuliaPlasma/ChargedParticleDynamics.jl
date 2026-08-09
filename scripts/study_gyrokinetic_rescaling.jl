#
# How do the gyrokinetic and the 4D guiding centre model relate?
#
# `GyroKinetics4d` and `GuidingCenter4d` describe the same physics — same equilibrium Hamiltonian,
# same generalised vector potential — but the gyrokinetic module integrates the characteristics
# with the phasespace Jacobian B*∥ cleared from the denominator, which is what makes its
# right-hand side divergence-free and its splitting volume preserving.
#
# The notes write that step as R̃ = B*∥ R, which reads like a change of variables. It is not one:
# it is the time reparametrisation dt = ωabs ds, and the two vector fields are simply proportional.
# This script establishes that, first pointwise and then by integrating both models and comparing
# where they end up.
#
# The factor is `ωabs = J B*∥`, the phasespace Jacobian, and it is positive here as in every chart.
# It is not `B*∥`: it contracts the coordinate curl `∂₂ϑ₃ - ∂₃ϑ₂`, which carries the *signed*
# Jacobian, against b, so the bare contraction inherits the handedness of the chart — and the
# Solov'ev chart is left-handed. Each module therefore applies the `orientation()` that
# `ElectromagneticFields` generates for its chart, which restores it and is what keeps the factor a
# measure density rather than a signed one.
#
# That sign was not cosmetic. Until it was applied, the vector field carried the negative factor and
# this model integrated the six left-handed-chart equilibria *backwards* — the same physical particle
# circulated one way in the cartesian chart of a tokamak and the other way in its cylindrical chart.
# ElectromagneticFields 0.7.0 is what exposed it: before that the reversed b cancelled the sign and
# the bare contraction came out positive everywhere, which is why it used to be written as B*∥ below.
#
# See docs/src/findings.md, "The gyrokinetic model is the guiding centre model in rescaled time" and
# "The orientation shows through anywhere a coordinate curl appears".
#
# Run with:  julia --project=scripts scripts/study_gyrokinetic_rescaling.jl
#

using ChargedParticleDynamics
using GeometricIntegrators
using Printf

const GK = ChargedParticleDynamics.GyroKinetics4d.GuidingCenter4dSolovevIterXpoint
const GC = ChargedParticleDynamics.GuidingCenter4d.SolovevIterXpoint

"v_gyrokinetic(q) = ωabs(q) · v_guiding_centre(q), pointwise. See the header on ωabs's orientation."
function pointwise()
    println("Pointwise:  max |v_gk - ωabs · v_gc| / |v_gk|\n")
    @printf("  %-38s %-14s %s\n", "q = (R₁, R₂, R₃, u)", "ωabs", "rel. difference")

    params = (μ = 1E-2,)
    for q in ([6.2, 0.3, 0.0, 3.4E-1], [5.5, -0.8, 1.1, -2.0E-1], [7.0, 0.5, 2.0, 5.0E-1])
        vgk = zeros(4); GK.v(vgk, 0.0, q, params)
        vgc = zeros(4); GC.guiding_center_4d_v(vgc, 0.0, q, params)
        ωfac = GK.ωabs(0.0, q)
        rel = maximum(abs.(vgk .- ωfac .* vgc) ./ max.(abs.(vgk), 1E-30))
        @printf("  %-38s %-14.6e %.2e\n", string(round.(q, digits = 2)), ωfac, rel)
    end
end

"""
The six subsystems of the splitting must reconstruct the full vector field exactly — this is a
property of the decomposition, not an approximation.
"""
function splitting()
    params = (μ = 1E-2,)
    q = [6.2, 0.3, 0.0, 3.4E-1]

    vfull = zeros(4); GK.v(vfull, 0.0, q, params)
    vsum  = zeros(4)
    for V in (GK.v₁, GK.v₂, GK.v₃, GK.v₄, GK.v₅, GK.v₆)
        vᵢ = zeros(4); V(vᵢ, 0.0, q, params); vsum .+= vᵢ
    end

    @printf("\nSplitting:  max |Σᵢ vᵢ - v| = %.2e\n", maximum(abs.(vsum .- vfull)))
end

"""
Integrate the gyrokinetic model over an interval of s and the guiding centre model over
t = ωabs(q₀) s, and compare the endpoints.

The agreement is limited by the variation of ωabs along the orbit, not by the integrator: the exact
relation is t = ∫₀ˢ ωabs(q(s')) ds', so replacing it by ωabs(q₀) s leaves an O(s²) error. The residual
therefore shrinks quadratically with the interval and is independent of the step.
"""
function orbits()
    println("\nOrbit equivalence: GK over s against GC over t = ωabs(q₀) s\n")
    @printf("  %-12s %-14s %s\n", "s interval", "endpoint diff", "ratio")

    q₀, par = GK.initial_conditions_deeply_passing()
    params  = (μ = par.μ,)
    ωfac    = GK.ωabs(0.0, q₀)

    previous = NaN
    for s_end in (4E-3, 2E-3, 1E-3, 5E-4)
        gk = integrate(GK.odeproblem(q₀; parameters = params, timestep = s_end / 2000, timespan = (0.0, s_end)), Gauss(2))
        gc = integrate(GC.odeproblem(q₀; parameters = params, timestep = ωfac * s_end / 2000, timespan = (0.0, ωfac * s_end)), Gauss(2))
        d  = maximum(abs.(gk.q[end] .- gc.q[end]))
        @printf("  %-12.0e %-14.2e %s\n", s_end, d, isnan(previous) ? "" : @sprintf("%.1f", previous / d))
        previous = d
    end
end

pointwise()
splitting()
orbits()
