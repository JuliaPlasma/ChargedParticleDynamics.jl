#
# How well conditioned is the constrained canonical guiding centre model, and does the choice of
# constraint formulation matter?
#
# The 3D model replaces the parallel-velocity constraint by two of the three independent components
# of v × b = 0. Li, Zhang & Liu label them
#
#     g¹ = b₁ v₃ - b₃ v₁,    g² = b₂ v₃ - b₃ v₂,    g³ = b₁ v₂ - b₂ v₁
#
# and show that the Poisson bracket {g₁, g₂} in the denominator of both Lagrange multipliers carries
# a factor of the b-component belonging to the constraint the pair leaves out:
#
#     pair (g³, g¹) → +b₁ [B + (p-A)·(∇×b)]    selected by  constraints = :g31
#     pair (g¹, g²) → +b₃ [ ... ]                           constraints = :g12
#     pair (g², g³) → -b₂ [ ... ]                           constraints = :g23
#
# The minus on the last of those comes from the ordering of that pair rather than from anything
# physical — reversing a pair flips λₒ and both multipliers together — so only |λₒ| is meaningful in
# the table below. It is spelled out because λₒ and bₘ then do not always share a sign.
#
# The package used to implement only (g³, g¹), so it was singular wherever b₁ = 0. In cylindrical
# coordinates b₁ = b_R, which vanishes on the midplane of an axisymmetric equilibrium — exactly where
# the supplied initial conditions sit — and seven of the eleven equilibria below could not be started
# at all. All three pairs are now available, and each equilibrium declares a `default_constraints`
# that is regular at its own initial conditions.
#
# On top of the choice of pair there is the choice of formulation: the Hamilton-Dirac equations
# (`hodeproblem`, Eq. 11), the canonicalised system with the extended Hamiltonian
# (`hodeproblem_canonical`, Eq. 13), and the compact form that drops the constraint-proportional
# terms from the right-hand side (`hodeproblem_compact`, Eq. 29). The last of these is regular
# everywhere and does not depend on the pair at all; see the derivation at the top of
# `src/guiding_center_3d/guiding_center_3d_compact.jl`.
#
# Section 1 measures the conditioning of each pair at each equilibrium's initial condition, section 2
# checks that the formulations agree wherever they are all usable, and section 3 measures what each
# one costs.
#
# See docs/src/findings.md, "Conditioning of the 3D guiding centre constraint formulations".
#
# Sections 2 and 3 integrate every equilibrium with every variant, and each variant is a distinct
# specialisation of the integrator, so most of the wall clock is compilation: expect several minutes.
# The script flushes after each equilibrium so progress is visible when its output is redirected.
#
# Run with:  julia --project=scripts scripts/study_guiding_center_3d_conditioning.jl
#

using ChargedParticleDynamics
using GeometricIntegrators
using LinearAlgebra
using Logging
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

const PAIRS = (:g31, :g12, :g23)

# The tolerance the test suite settled on; see the comment on `options` in
# `test/guiding_center_3d_tests.jl`.
const OPTIONS = (f_abstol = 1E-12, max_iterations = 50)

const METHOD = PartitionedGauss(2)

# A hundred steps at the smaller of the equilibrium's own `DEFAULT_TIMESTEP` and 0.1, which is the
# workload the test suite uses. Neither bound can be dropped. A single step size across all eleven is
# not usable — forcing Δt = 0.1 on the dipole, whose natural step is 0.03, makes every formulation
# diverge, which says something about that step size and nothing about the formulations. But taking
# each equilibrium's own step is not usable either: `TokamakSmallCartesian` integrates at Δt = 200,
# where every solve is stiff enough that its six variants alone outlast the whole rest of this script.
const NSTEPS = 100

workload(M) = (timespan = (0.0, NSTEPS * min(M.DEFAULT_TIMESTEP, 0.1)),
               timestep = min(M.DEFAULT_TIMESTEP, 0.1))

mean(x) = isempty(x) ? 0.0 : sum(x) / length(x)

# `DataSeries` is indexable but not an iterator, so `maximum(abs, ds)` does not work on one.
mx(ds) = maximum(abs(ds[i]) for i in eachindex(ds))

# every GuidingCenter3d `initial_conditions_*` returns a NamedTuple (q, p, params)
initial(M, icsname) = getfield(M, icsname)()


# ---------------------------------------------------------------------------------------------
# A logger that harvests the iteration count out of the "Solver took N iterations." warnings, as in
# `study_solver_tolerances.jl`. `SimpleSolvers` emits one per solve that reaches `warn_iterations`,
# so `warn_iterations=1` turns the warning into a per-solve report of the work a step cost.
# ---------------------------------------------------------------------------------------------

mutable struct IterationLogger <: AbstractLogger
    iterations::Vector{Int}
end

Logging.shouldlog(::IterationLogger, level, _module, group, id) = true
Logging.min_enabled_level(::IterationLogger) = Logging.Debug
Logging.catch_exceptions(::IterationLogger) = false

function Logging.handle_message(logger::IterationLogger, level, message, _module, group, id,
                                file, line; kwargs...)
    m = match(r"Solver took (\d+) iterations", string(message))
    m === nothing || push!(logger.iterations, parse(Int, m.captures[1]))
    nothing
end

const ITERLOG = IterationLogger(Int[])


# ---------------------------------------------------------------------------------------------
# 1. Conditioning of the three pairs
# ---------------------------------------------------------------------------------------------

function conditioning()
    println("""
    1. Conditioning of the three constraint pairs
    =============================================

    At the shipped initial condition of each equilibrium:

      b₁ b₂ b₃   the components of b in that equilibrium's chart; the pair omitting gᵏ divides by bₖ
      λₒ         the Poisson bracket {g₁, g₂} that divides both Lagrange multipliers, per pair
      D          the same bracket as the compact form computes it, Σₘ bₘ {cᵢ,cⱼ}(m) / Σₘ bₘ bₘ, which
                 never divides by a single component and so stays finite where the pairs do not
      default    the pair the equilibrium ships with
    """)

    @printf("  %-26s %-12s %10s %10s %10s   %11s %11s %11s   %10s %8s\n",
            "equilibrium", "coords", "b₁", "b₂", "b₃",
            "λₒ(:g31)", "λₒ(:g12)", "λₒ(:g23)", "D", "default")

    for (name, coords, M, icsname) in CASES
        isdefined(M, icsname) || continue
        ic = initial(M, icsname)
        q, p, t = ic.q, ic.p, 0.0

        λ = [M.λₒ(t, q, p, M.constraint_pair(s)) for s in PAIRS]

        @printf("  %-26s %-12s %10.3e %10.3e %10.3e   %11.3e %11.3e %11.3e   %10.3e %8s\n",
                name, coords, M.b₁(t, q), M.b₂(t, q), M.b₃(t, q), λ...,
                M.compact_denominator(t, q, p), M.default_constraints())
    end
    flush(stdout)

    println("""

    A vanishing λₒ means both multipliers are ±Inf or NaN and the problem cannot be started: the
    Newton solver hits the NaN in its direction vector on the first step. Seven of the eleven
    equilibria are in that position for the pair (g³, g¹) that the package used to hard-code, and
    `SolovevSymmetricField` is in it for two of the three pairs at once — b = e₂ there, so only
    (g², g³) survives.

    Nothing about this is confined to the curvilinear equilibria. `TokamakSmallCartesian` fails for
    (g³, g¹) too, because its initial condition sits at y = z = 0 where B_x = 0. What decides the
    outcome is where the initial condition lies relative to the field, not the coordinate system.

    The D column is the point of the compact formulation: it is finite for every equilibrium here,
    including the ones where two of the three pairs are singular. The three ratios {cᵢ,cⱼ}(m) / bₘ
    that it averages agree to machine precision wherever they are all defined, in every chart, which
    is what makes the average a faithful stand-in rather than a fudge.""")
end


# ---------------------------------------------------------------------------------------------
# 2. Do the formulations agree?
# ---------------------------------------------------------------------------------------------

"""
    variants(M, ic)

The variants to compare for equilibrium module `M` at initial condition `ic`, as
`(label, regular, constructor)` triples: the Hamilton-Dirac form for each of the three pairs, the
canonicalised form and the compact form for the equilibrium's own default pair, and the compact form
in its regularised, pair-independent variant.

`regular` says whether the variant's denominator survives the initial condition. Attempting a
singular one is not merely useless but slow: the multipliers are infinite from the first stage, and
the backtracking line search grinds through its thousand-iteration budget on every step before the
solver finally reports the NaN.
"""
function variants(M, ic)
    d = M.default_constraints()
    w = workload(M)
    q, p, t = ic.q, ic.p, 0.0

    # a pair is usable where the bracket it divides by does not vanish …
    pairok(s) = (λ = M.λₒ(t, q, p, M.constraint_pair(s)); isfinite(λ) && !iszero(λ))
    # … and the literal compact form where the single component of b it divides by does not
    compactok(s) = !iszero(M.bᵢ(M.compact_index(s), t, q))

    (("hode :g31", pairok(:g31), ic -> M.hodeproblem(ic; constraints = :g31, w...)),
     ("hode :g12", pairok(:g12), ic -> M.hodeproblem(ic; constraints = :g12, w...)),
     ("hode :g23", pairok(:g23), ic -> M.hodeproblem(ic; constraints = :g23, w...)),
     ("canonical :$d", pairok(d), ic -> M.hodeproblem_canonical(ic; w...)),
     ("compact :$d", compactok(d), ic -> M.hodeproblem_compact(ic; constraints = d, w...)),
     ("compact :parallel", true, ic -> M.hodeproblem_compact(ic; w...)))
end

"""
    solve(M, ic, makeproblem)

Integrate one variant and return the final state together with the invariant errors along it, or
`nothing` if it could not be integrated at all.
"""
function solve(M, ic, makeproblem)
    try
        problem = makeproblem(ic)
        sol = integrate(problem, METHOD; OPTIONS...)
        y = vcat(sol.q[end], sol.p[end])
        any(!isfinite, y) && return nothing

        c = M.compute_constraints(sol)
        _, e = M.compute_energy_error(sol.t, sol.q, sol.p, parameters(problem))

        return (y = y,
                energy = mx(e),
                constraints = maximum(mx(getproperty(c, k)) for k in (:g₁, :g₂, :g₃)))
    catch
        return nothing
    end
end

function agreement()
    println("""


    2. Do the formulations agree?
    =============================

    Each variant integrated for $NSTEPS steps at min(`DEFAULT_TIMESTEP`, 0.1) with
    $(nameof(typeof(METHOD))):

      Δy       relative deviation of the final state (q, p) from the reference, ‖Δy‖/‖y‖
      |ΔH/H|   largest relative energy error along the orbit
      max|gᵏ|  largest of all three constraints along the orbit, whichever pair the variant retained

    The reference is the compact form in its `:parallel` variant, which is the only one of the six
    defined for every equilibrium — it depends on no constraint pair and divides by no component of b.
    """)

    for (name, _, M, icsname) in CASES
        isdefined(M, icsname) || continue
        ic = initial(M, icsname)

        results = [(label, regular, regular ? solve(M, ic, makeproblem) : nothing)
                   for (label, regular, makeproblem) in variants(M, ic)]

        reference = let r = last(results)[3]
            r === nothing ? nothing : r.y
        end

        println("  $name")
        @printf("    %-20s %12s %12s %12s\n", "variant", "Δy", "|ΔH/H|", "max|gᵏ|")

        for (label, regular, r) in results
            if r === nothing
                @printf("    %-20s %12s %12s %12s\n", label, regular ? "diverged" : "singular", "-", "-")
            else
                Δy = reference === nothing ? NaN : norm(r.y - reference) / norm(reference)
                @printf("    %-20s %12.3e %12.3e %12.3e\n", label, Δy, r.energy, r.constraints)
            end
        end
        println()
        flush(stdout)
    end

    println("""
    Wherever two variants are both well conditioned they track each other to the accuracy of the
    nonlinear solve, which is the substantive check on the new formulations: the three pairs
    parameterise the same constrained system, and the compact form agrees with the Hamilton-Dirac one
    it was derived from. `TokamakMediumCartesian` is the sharpest case, being the one equilibrium
    where no component of b vanishes at the initial condition and all three pairs are therefore
    usable at once.

    Being *usable* and being *well conditioned* are not the same thing, though, and `Dipole3d` is
    where that shows. All three pairs are regular at its initial condition, but the orbit takes b₁
    and b₂ through zero, and the two pairs that divide by them hold the constraints only to 1E-3 and
    the energy to 4E-5 — four orders of magnitude worse than (g¹, g²), which divides by b₃ and holds
    them to 2E-8 and 2E-9. Nothing is singular there; the multipliers are simply large. Conditioning
    along the orbit, not merely at the initial condition, is what the choice of pair should follow.""")
end


# ---------------------------------------------------------------------------------------------
# 3. What does each formulation cost?
# ---------------------------------------------------------------------------------------------

"""
    cost(makeproblem, ic)

Per-step cost of one variant and the Newton iteration statistics behind it. A first integration is
discarded so that compilation does not land in the timing. Returns `nothing` for a variant that
cannot be integrated.
"""
function cost(makeproblem, ic)
    try
        problem = makeproblem(ic)
        nsteps = GeometricIntegrators.GeometricBase.ntime(problem)
        integrate(problem, METHOD; warn_iterations = 1, OPTIONS...)
        empty!(ITERLOG.iterations)
        elapsed = @elapsed integrate(problem, METHOD; warn_iterations = 1, OPTIONS...)
        its = copy(ITERLOG.iterations)
        return (ms_per_step = elapsed / nsteps * 1E3,
                mean_iterations = mean(its),
                peak_iterations = isempty(its) ? 0 : maximum(its))
    catch
        return nothing
    end
end

function costs()
    println("""


    3. What does each formulation cost?
    ===================================

    Per-step wall clock and Newton iteration counts, same workload as section 2. Single sample from
    `@elapsed` after a discarded warm-up run, so treat differences below ~20 % as noise.
    """)

    @printf("  %-26s %-20s %10s %10s %8s\n",
            "equilibrium", "variant", "ms/step", "mean its", "peak")

    for (name, _, M, icsname) in CASES
        isdefined(M, icsname) || continue
        ic = initial(M, icsname)

        for (label, regular, makeproblem) in variants(M, ic)
            w = if regular
                with_logger(ITERLOG) do
                    cost(makeproblem, ic)
                end
            end
            if w === nothing
                @printf("  %-26s %-20s %10s %10s %8s\n", name, label,
                        regular ? "diverged" : "singular", "-", "-")
            else
                @printf("  %-26s %-20s %10.3f %10.1f %8d\n",
                        name, label, w.ms_per_step, w.mean_iterations, w.peak_iterations)
            end
        end
        flush(stdout)
    end

    println("""

    Between the three pairs there is nothing to choose on cost — they are the same expressions with
    permuted indices, and they come out within a few percent of each other everywhere — so the pair
    should be picked on conditioning alone.

    Between the formulations there is a great deal to choose. Taking the Hamilton-Dirac form as unity:

      compact, literal    0.86-0.95   it drops the multipliers for a single bracket ratio
      compact, :parallel  1.21-1.44   the regularised denominator averages all three brackets, so it
                                      evaluates nine bracket terms where the literal form evaluates three
      canonicalised       8.9-18.8    ∂λ/∂q and ∂λ/∂p need every second derivative of the field

    The order of magnitude on the canonicalised form is the number worth remembering: it buys exact
    symplecticity rather than approximate, and section 2 shows it giving up as much as two or three
    orders of magnitude of energy conservation and four or five of constraint conservation for it at
    these step sizes. Those are worst cases and not typical ones — on five of the eleven it is within a
    factor of two on both, and on TokamakSmallCartesian it conserves energy slightly better.

    Its cost spread is wider than the others' because part of the cost is the harder nonlinear solve
    rather than the right-hand side — Dipole3d, the worst case, is also the only equilibrium whose mean
    iteration count exceeds one. The regularisation of the compact form costs about half as much again
    as the literal one, which is a cheap price for never dividing by a component of b.""")
end


function main()
    conditioning()
    agreement()
    costs()
end

main()
