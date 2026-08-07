#
# Why did the ITER-like Solov'ev test blocks take minutes when everything else took seconds?
#
# The test modules used to request `f_abstol = 1E-15, f_reltol = 1E-15` from the nonlinear solver.
# `SimpleSolvers` accepts an iterate when
#
#     rfₐ ≤ f_abstol + f_reltol * ‖F(x₀)‖
#
# (`assess_convergence` in `nonlinear/nonlinear_solver_status.jl`), i.e. the usual `atol + rtol‖F₀‖`
# residual test. Setting *both* terms to 1E-15 collapses the right-hand side to ~1E-15 no matter how
# large the problem is, and in particular removes the relative term — which exists precisely so that
# a large-magnitude solve can converge once its residual has fallen by `f_reltol` from the initial
# one. On the ITER-scale problems nothing could then satisfy the criterion, and the solver ran to
# `max_iterations` (1000 by default) on nearly half of all steps, each iteration costing a
# ForwardDiff Jacobian and a backtracking line search.
#
# Two distinct things go wrong, and the sections below separate them:
#
#   * `f_abstol = 1E-15` is on its own too small for the 4D guiding centre in an ITER-scale field.
#     Its residual carries the momentum equation p = ϑ(q), and ‖ϑ‖ ≈ 15 there, so no residual can be
#     driven below `‖ϑ‖ eps ≈ 3.3E-15`. Section 1 tabulates that floor.
#
#   * The Pauli model is *not* explained that way — its momentum is small (‖p‖ ≈ 0.14, floor
#     3E-17) and `f_abstol = 1E-15` alone is perfectly fine for it. It only degrades when `f_reltol`
#     is pinned as well, i.e. it is the loss of the relative term that hurts, not the absolute
#     bound. Section 2 shows both effects side by side.
#
# Section 3 measures the per-step cost of each model/formulation/integrator once the tolerances are
# sound, which shows that the two slowest 3D blocks in CI were never an ITER effect at all.
#
# See the comment on `options` in `test/guiding_center_3d_tests.jl`.
#
# Run with:  julia --project=scripts scripts/study_solver_tolerances.jl
#

using ChargedParticleDynamics
using GeometricIntegrators
using LinearAlgebra
using Logging
using Printf

const GC3 = ChargedParticleDynamics.GuidingCenter3d
const GC4 = ChargedParticleDynamics.GuidingCenter4d
const P3 = ChargedParticleDynamics.PauliParticle3d

# The tolerance the test suite settled on, for reference in the tables below.
const F_ABSTOL = 1E-12

mean(x) = isempty(x) ? 0.0 : sum(x) / length(x)

# ---------------------------------------------------------------------------------------------
# A logger that harvests the iteration count out of the "Solver took N iterations." warnings.
#
# `SimpleSolvers` emits one per solve that reaches `warn_iterations`, so setting `warn_iterations=1`
# turns that warning into a per-solve report of how much work the step actually cost. The message is
# unchanged in `SimpleSolvers` 0.10, so this still works; the alternative is the qualified status API
# (`SimpleSolvers.status`, `SimpleSolvers.isstalled`), which is deliberately *not* exported and needs
# the solver to be driven a step at a time rather than through `integrate`.
#
# This is not why the test suite has `test/quiet_solver_warnings.jl`: that exists for the one block
# whose solves genuinely do not converge. See its header.
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

"""
    solver_work(problem, method; options...)

Integrate `problem` once and return the per-step cost together with the Newton iteration statistics
behind it. A first integration is discarded so that compilation does not land in the timing.
"""
function solver_work(problem, method; kwargs...)
    nsteps = GeometricIntegrators.GeometricBase.ntime(problem)
    integrate(problem, method; warn_iterations=1, kwargs...)
    empty!(ITERLOG.iterations)
    elapsed = @elapsed integrate(problem, method; warn_iterations=1, kwargs...)
    its = ITERLOG.iterations
    (ms_per_step = elapsed / nsteps * 1E3,
     nsteps = nsteps,
     mean_iterations = mean(its),
     peak_iterations = isempty(its) ? 0 : maximum(its),
     capped_steps = count(>=(1000), its))
end

# ---------------------------------------------------------------------------------------------
# 1. Magnitudes entering the residual, and the round-off floor they imply
# ---------------------------------------------------------------------------------------------

const SCALE_CASES = (("GC4d", GC4.SolovevIterXpoint),
                     ("GC4d", GC4.SolovevIter),
                     ("GC4d", GC4.SolovevSymmetricField),
                     ("GC4d", GC4.TokamakIterCylindrical),
                     ("GC4d", GC4.TokamakMediumCartesian),
                     ("GC4d", GC4.TokamakSmallCartesian),
                     ("Pauli3d", P3.SolovevIter),
                     ("Pauli3d", P3.SolovevIterXpoint),
                     ("Pauli3d", P3.TokamakIterCylindrical))

"""
    residual_scale(family, M)

The magnitude of the momentum appearing in the residual at `M`'s shipped
`initial_conditions_barely_passing`. For the 4D guiding centre that is the one-form ϑ = A + u b
evaluated at q₀; for the Pauli model the momentum p₀ is supplied directly. Returns `nothing` if the
equilibrium does not ship that initial condition.
"""
function residual_scale(family, M)
    isdefined(M, :initial_conditions_barely_passing) || return nothing
    ic = M.initial_conditions_barely_passing()

    if family == "GC4d"
        q = ic[1]
        x, u = q[1:3], q[4]
        A = [M.A₁(0.0, x), M.A₂(0.0, x), M.A₃(0.0, x)]
        b = [M.b₁(0.0, x), M.b₂(0.0, x), M.b₃(0.0, x)]
        return (x = x, normA = norm(A), normp = norm(A .+ u .* b), what = "‖ϑ‖")
    else
        q, p = ic[1], ic[2]
        A = [M.A₁(0.0, q), M.A₂(0.0, q), M.A₃(0.0, q)]
        return (x = q, normA = norm(A), normp = norm(p), what = "‖p‖")
    end
end

function residual_scales()
    println("""
    1. Magnitudes entering the residual
    ===================================

    Evaluated at each equilibrium's shipped `initial_conditions_barely_passing`. Note that the
    Solov'ev modules work in coordinates normalised by R₀, so their q₀ is ~0.4 rather than 2.5.

    An *absolute* tolerance below the last column is unreachable for a residual of that magnitude.
    """)

    @printf("  %-8s %-26s %-20s %10s %12s %12s  %s\n",
            "family", "equilibrium", "q₀", "‖A‖", "‖ϑ‖ or ‖p‖", "·eps", "f_abstol=1E-15")

    for (family, M) in SCALE_CASES
        s = residual_scale(family, M)
        s === nothing && continue
        floor = s.normp * eps(Float64)
        @printf("  %-8s %-26s %-20s %10.4g %12.4g %12.3e  %s\n",
                family, string(nameof(M)), string(round.(s.x; digits = 3)),
                s.normA, s.normp, floor, floor < 1E-15 ? "reachable" : "UNREACHABLE")
    end

    println("""

    For the 4D guiding centre this is the whole story: the ITER-like Solov'ev equilibria and
    `TokamakIterCylindrical` all sit above 1E-15 because ‖A‖ ~ B₀R₀² and ITER has B₀ = 5.3,
    R₀ = 6.2, while the medium and small tokamaks sit comfortably below it.
    `SolovevSymmetricField` is by far the worst at ‖ϑ‖ ≈ 256, and its difficulty is a step-size
    limit rather than a tolerance one: at f_abstol = 1E-12, `Gauss(2)` on its `odeproblem` averages
    29.4 iterations and reaches the iteration cap on 45 of 100 steps at Δt = 1E2, against 2.2
    iterations and no capped step at the Δt = 1E0 it runs at.

    The Pauli rows are the counterexample that stops this being a complete explanation: their
    momentum is small, so this floor says `f_abstol = 1E-15` should have been fine — and measured on
    its own, it is. What broke the Pauli blocks was pinning `f_reltol` as well, which section 2
    separates out.

    Either way, this rules out simply dropping the `options` from the test suite.
    `GeometricIntegratorsBase.default_options` scales its tolerance with the stage system,
    `f_abstol = max(8, solversize(method, problem)) * eps(datatype(problem))`, which is 1.8E-15 for
    the 4D guiding centre and Pauli solves and 2.7E-15 for the 3D `hodeproblem` — both still below
    the 4D guiding centre's ITER floor.""")
end

# ---------------------------------------------------------------------------------------------
# 2. Which knob is responsible
# ---------------------------------------------------------------------------------------------

# (label, options) — `f_reltol` unset means it keeps its default of √eps.
const TOLERANCE_SETTINGS = (("f_abstol=1E-15, f_reltol=1E-15", (f_abstol = 1E-15, f_reltol = 1E-15)),
                            ("f_abstol=1E-15", (f_abstol = 1E-15,)),
                            ("f_abstol=1E-14", (f_abstol = 1E-14,)),
                            ("f_abstol=1E-12", (f_abstol = 1E-12,)),
                            ("f_abstol=1E-12, f_reltol=1E-12", (f_abstol = 1E-12, f_reltol = 1E-12)))

const SWEEP_CASES = (("GC4d SolovevIterXpoint iode",
                      () -> GC4.SolovevIterXpoint.iodeproblem(
                          GC4.SolovevIterXpoint.initial_conditions_barely_passing();
                          timespan = (0.0, 20.0), timestep = 0.1), VPRKGauss(2)),
                     ("Pauli3d SolovevIter iode",
                      () -> P3.SolovevIter.iodeproblem(
                          P3.SolovevIter.initial_conditions_barely_passing();
                          timespan = (0.0, 20.0), timestep = 0.1), VPRKGauss(2)),
                     ("GC4d TokamakMediumCartesian iode",
                      () -> GC4.TokamakMediumCartesian.iodeproblem(
                          GC4.TokamakMediumCartesian.initial_conditions_barely_passing();
                          timespan = (0.0, 20.0), timestep = 0.1), VPRKGauss(2)))

function tolerance_sweep()
    println("""


    2. Which knob is responsible
    ============================
    """)

    for (label, makeproblem, method) in SWEEP_CASES
        println("  $label")
        @printf("    %-34s %12s %12s %10s %14s\n",
                "setting", "ms/step", "mean iters", "peak", "steps at cap")
        for (setting, opts) in TOLERANCE_SETTINGS
            w = solver_work(makeproblem(), method; opts...)
            @printf("    %-34s %12.3f %12.1f %10d %10d/%d\n",
                    setting, w.ms_per_step, w.mean_iterations, w.peak_iterations,
                    w.capped_steps, w.nsteps)
        end
        println()
    end

    println("""    The 4D guiding centre degrades from `f_abstol = 1E-15` alone, and recovers between
    1E-15 and 1E-14 — exactly where section 1 puts its floor. The Pauli model is untouched by
    `f_abstol = 1E-15` and only degrades once `f_reltol = 1E-15` is added on top, so for it the
    damage is the loss of the relative term. The medium tokamak is flat across every setting: its
    floor is 2.6E-16, so even the tightest combination was reachable.

    `f_abstol = 1E-12` with `f_reltol` left alone is the only setting that is comfortable for all
    three, which is what the test suite now uses. The final state is unchanged to ~1E-14 across
    every row — the extra hundreds of iterations per step buy nothing at all.""")
end

# ---------------------------------------------------------------------------------------------
# 3. Per-step cost by model, formulation and integrator
# ---------------------------------------------------------------------------------------------

const COST_CASES = (("GC3d SolovevIterXpoint hodeproblem +extrapolation",
                     () -> GC3.SolovevIterXpoint.hodeproblem(
                         GC3.SolovevIterXpoint.initial_conditions_barely_passing();
                         timespan = (0.0, 10.0), timestep = 0.1), PartitionedGauss(2),
                     (initialguess = MidpointExtrapolation(5),)),
                    ("GC3d SolovevIterXpoint hodeproblem",
                     () -> GC3.SolovevIterXpoint.hodeproblem(
                         GC3.SolovevIterXpoint.initial_conditions_barely_passing();
                         timespan = (0.0, 10.0), timestep = 0.1), PartitionedGauss(2), NamedTuple()),
                    ("GC3d TokamakMediumCartesian hodeproblem +extrapolation",
                     () -> GC3.TokamakMediumCartesian.hodeproblem(
                         GC3.TokamakMediumCartesian.initial_conditions_barely_passing();
                         timespan = (0.0, 10.0), timestep = 0.1), PartitionedGauss(2),
                     (initialguess = MidpointExtrapolation(5),)),
                    ("GC3d TokamakMediumCartesian hodeproblem",
                     () -> GC3.TokamakMediumCartesian.hodeproblem(
                         GC3.TokamakMediumCartesian.initial_conditions_barely_passing();
                         timespan = (0.0, 10.0), timestep = 0.1), PartitionedGauss(2), NamedTuple()),
                    ("GC4d SolovevIterXpoint ode",
                     () -> GC4.SolovevIterXpoint.odeproblem(
                         GC4.SolovevIterXpoint.initial_conditions_barely_passing();
                         timespan = (0.0, 200.0), timestep = 1.0), Gauss(2), NamedTuple()),
                    ("GC4d SolovevIterXpoint iode",
                     () -> GC4.SolovevIterXpoint.iodeproblem(
                         GC4.SolovevIterXpoint.initial_conditions_barely_passing();
                         timespan = (0.0, 20.0), timestep = 0.1), VPRKGauss(2), NamedTuple()),
                    ("Pauli3d SolovevIter hodeproblem",
                     () -> P3.SolovevIter.hodeproblem(
                         P3.SolovevIter.initial_conditions_barely_passing();
                         timespan = (0.0, 10.0), timestep = 0.1), PartitionedGauss(2), NamedTuple()),
                    ("Pauli3d SolovevIter iode",
                     () -> P3.SolovevIter.iodeproblem(
                         P3.SolovevIter.initial_conditions_barely_passing();
                         timespan = (0.0, 10.0), timestep = 0.1), VPRKGauss(2), NamedTuple()))

function step_costs()
    println("""


    3. Per-step cost at the tolerance the suite now uses (f_abstol = $(F_ABSTOL))
    ============================================================================
    """)

    @printf("    %-50s %12s %12s\n", "workload", "ms/step", "mean iters")
    for (label, makeproblem, method, extra) in COST_CASES
        try
            w = solver_work(makeproblem(), method; f_abstol = F_ABSTOL, extra...)
            @printf("    %-50s %12.3f %12.1f\n", label, w.ms_per_step, w.mean_iterations)
        catch e
            @printf("    %-50s %s\n", label, first(sprint(showerror, e), 40))
        end
    end

    println("""

    Once the tolerance is sound, the remaining spread has nothing to do with ITER. The 3D `hodeproblem`
    path costs orders of magnitude more per step than the Pauli `hodeproblem` even though Newton converges
    in one iteration, and `TokamakMediumCartesian` is *slower* than the ITER Solov'ev — so the two
    ten-minute 3D blocks in CI were a step-count and per-step-cost problem, not a convergence one.

    This used to say the cost was nested ForwardDiff over the second derivatives of the Hamiltonian,
    re-evaluated by the line search. None of that was true of this path. `ElectromagneticFields`
    generates its field functions symbolically with SymEngine at precompile time and does not depend on
    ForwardDiff; `hodeproblem` never evaluates a second derivative of the Hamiltonian, only
    `hodeproblem_canonical` does; and of the ~98 right-hand side evaluations per step exactly 4 are on
    `Dual`, the line search being 8% of wall clock against the Jacobian's 13%.

    The cost was `MidpointExtrapolation(5)` at 70% of wall clock — not the default for either method, and
    dropped — over a right-hand side that re-entered the generated field code once per bracket term. On
    ElectromagneticFields 0.6.2, `db₁dx₁` was 1905 statements with 108 separate evaluations of `log(x₁)`
    for the ITER Solov'ev X-point, and was called seventy times per evaluation; it is now called once.

    The third cost was the expression swell in that generated code, and lived upstream. It is delivered:
    ElectromagneticFields 0.6.3 eliminates common subexpressions in the bodies it emits, taking `db₁dx₁`
    from 549 ns to 57 and `d²b₁dx₁dx₁` from 1733 to 80, and changing no value — every field function
    this package injects is bit-identical between 0.6.2 and 0.6.3.

    `Project.toml` now requires 0.7.0, which unlike 0.6.3 is *not* value-preserving: it reverses `b` in
    the four left-handed charts, so six of this package's equilibria — every cylindrical, toroidal and
    Solov'ev one bar `SolovevSymmetricField`, which is cartesian despite the name — integrate a
    different orbit than they did. It changes the sign of terms in the generated code rather than the
    amount of arithmetic, so the timings above re-measure inside their previous spread. See
    docs/src/findings.md, "The reversed field of the left-handed charts".

    `docs/src/findings.md` measures all four combinations of the two changes, because they are not the
    same factor: for `hodeproblem` they are separable and multiplicative, 6.9x here and 2.8x upstream,
    while for `hodeproblem_canonical` they reinforce each other, the second derivatives that form needs
    having gained the most from the elimination.""")
end

function main()
    logger = global_logger(ITERLOG)
    try
        residual_scales()
        tolerance_sweep()
        step_costs()
    finally
        global_logger(logger)
    end
end

main()
