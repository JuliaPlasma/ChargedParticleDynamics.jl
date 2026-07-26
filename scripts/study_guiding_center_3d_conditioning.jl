#
# How well conditioned is the constrained canonical guiding centre model?
#
# The 3D model replaces the parallel-velocity constraint by two of the three independent components
# of v × b = 0. Li, Zhang & Liu label them
#
#     g¹ = b₁ v₃ - b₃ v₁,    g² = b₂ v₃ - b₃ v₂,    g³ = b₁ v₂ - b₂ v₁
#
# and show that the Poisson bracket in the denominator of both Lagrange multipliers is
# b₃ [B + (p-A)·(∇×b)] for the pair (g¹, g²) and b₁ [ ... ] for the pair (g¹, g³).
#
# This package uses (g³, g¹), so the multipliers are singular where b₁ = 0. In cylindrical
# coordinates b₁ = b_R, which vanishes on the midplane of an axisymmetric equilibrium — exactly
# where the supplied initial conditions sit. This script measures how close to that the shipped
# initial conditions are, and what it does to the multipliers.
#
# See docs/src/findings.md, "Conditioning of the 3D guiding centre constraint pair".
#
# Run with:  julia --project=scripts scripts/study_guiding_center_3d_conditioning.jl
#

using ChargedParticleDynamics
using Printf

const G3 = ChargedParticleDynamics.GuidingCenter3d

const CASES = (("Dipole3d",                "cartesian",   G3.Dipole3d,                :initial_conditions_dipole),
               ("QuadraticPotentials3d",   "cartesian",   G3.QuadraticPotentials3d,   :initial_conditions_quadratic),
               ("TokamakSmallCartesian",   "cartesian",   G3.TokamakSmallCartesian,   :initial_conditions_barely_passing),
               ("TokamakMediumCartesian",  "cartesian",   G3.TokamakMediumCartesian,  :initial_conditions_barely_passing),
               ("TokamakSmallCylindrical", "cylindrical", G3.TokamakSmallCylindrical, :initial_conditions_barely_passing),
               ("TokamakMediumCylindrical","cylindrical", G3.TokamakMediumCylindrical,:initial_conditions_barely_passing),
               ("TokamakIterCylindrical",  "cylindrical", G3.TokamakIterCylindrical,  :initial_conditions_barely_passing),
               ("SolovevIter",             "cylindrical", G3.SolovevIter,             :initial_conditions_barely_passing),
               ("SolovevIterXpoint",       "cylindrical", G3.SolovevIterXpoint,       :initial_conditions_barely_passing),
               ("SolovevSymmetricField",   "cylindrical", G3.SolovevSymmetricField,   :initial_conditions_barely_passing),
               ("TokamakSmallToroidal",    "toroidal",    G3.TokamakSmallToroidal,    :initial_conditions_barely_passing))

function main()
    println("""
    At the shipped initial condition of each equilibrium:

      b₁    the component whose vanishing makes the active constraint pair (g³, g¹) singular
      λₒ    the Poisson bracket {g₃, g₁} that divides both Lagrange multipliers
      λ₁,λ₂ the multipliers themselves
    """)

    @printf("  %-27s %-12s %-12s %-13s %-12s %-12s\n",
            "equilibrium", "coords", "b₁", "λₒ", "λ₁", "λ₂")

    for (name, coords, M, icsname) in CASES
        isdefined(M, icsname) || continue
        # every GuidingCenter3d `initial_conditions_*` returns a NamedTuple (q, p, params)
        ic = getfield(M, icsname)()
        q, p, par = ic.q, ic.p, ic.params
        t = 0.0

        b₁ = M.b₁(t, q)
        λₒ = M.λₒ(t, q, p)
        λ₁ = M.λ₁(t, q, p, par)
        λ₂ = M.λ₂(t, q, p, par)

        @printf("  %-27s %-12s %-12.3e %-13.3e %-12.3e %-12.3e\n", name, coords, b₁, λₒ, λ₁, λ₂)
    end

    println("""

    Seven of the eleven equilibria have b₁ = 0 *exactly* at their shipped initial condition, so λₒ
    vanishes and both multipliers are ±Inf or NaN. The model cannot be started there at all: `hodeproblem`
    fails immediately with a NaN in the Newton direction.

    This is not confined to the curvilinear equilibria. TokamakSmallCartesian fails too — its
    initial condition sits at y = z = 0, where B_x = 0 — so the deciding factor is where the
    initial condition is placed relative to the field, not the coordinate system. What the four
    survivors have in common is simply that their initial condition is off the symmetry plane.

    That the survivors are exactly the equilibria the test suite exercises, plus the two toy
    problems that have their own scripts, is not a coincidence: the remaining 3D test blocks were
    commented out because of this, and `scripts/guiding_center_3d.jl` displaces its initial
    condition off the midplane by sqrt(eps()) to work around it.

    Choosing the pair (g¹, g²) instead would put b₃ = b_φ in the denominator, which is the largest
    component of a tokamak field and does not vanish on the midplane. The choice is currently
    hard-coded in `guiding_center_3d_equations.jl`, with the alternative commented out beside it;
    making it selectable would make most of these equilibria usable.""")
end

main()
