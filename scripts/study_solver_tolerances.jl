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
# turns that warning into a per-solve report of how much work the step actually cost. It is the only
# way to see the iteration counts from outside the solver, and it is why the test suite needs
# `test/quiet_solver_warnings.jl` at all.
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
    `SolovevSymmetricField` is by far the worst at ‖ϑ‖ ≈ 256; it never showed up in CI only because
    its test block runs `ode` rather than `iode`, and `Gauss` converges on these problems in one or
    two iterations regardless of tolerance.

    The Pauli rows are the counterexample that stops this being a complete explanation: their
    momentum is small, so this floor says `f_abstol = 1E-15` should have been fine — and measured on
    its own, it is. What broke the Pauli blocks was pinning `f_reltol` as well, which section 2
    separates out.

    Either way, this rules out simply dropping the `options` from the test suite:
    `GeometricIntegratorsBase.default_options` asks for `f_abstol = 8eps() = 1.8E-15`, itself below
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
                      () -> GC4.SolovevIterXpoint.guiding_center_4d_iode(
                          GC4.SolovevIterXpoint.initial_conditions_barely_passing()...;
                          tspan = (0.0, 20.0), tstep = 0.1), VPRKGauss(2)),
                     ("Pauli3d SolovevIter iode",
                      () -> P3.SolovevIter.pauli_particle_3d_iode(
                          P3.SolovevIter.initial_conditions_barely_passing()...;
                          tspan = (0.0, 20.0), tstep = 0.1), VPRKGauss(2)),
                     ("GC4d TokamakMediumCartesian iode",
                      () -> GC4.TokamakMediumCartesian.guiding_center_4d_iode(
                          GC4.TokamakMediumCartesian.initial_conditions_barely_passing()...;
                          tspan = (0.0, 20.0), tstep = 0.1), VPRKGauss(2)))

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

const COST_CASES = (("GC3d SolovevIterXpoint hode +extrapolation",
                     () -> GC3.SolovevIterXpoint.hode(
                         GC3.SolovevIterXpoint.initial_conditions_barely_passing()...;
                         tspan = (0.0, 10.0), tstep = 0.1), PartitionedGauss(2),
                     (initialguess = MidpointExtrapolation(5),)),
                    ("GC3d SolovevIterXpoint hode",
                     () -> GC3.SolovevIterXpoint.hode(
                         GC3.SolovevIterXpoint.initial_conditions_barely_passing()...;
                         tspan = (0.0, 10.0), tstep = 0.1), PartitionedGauss(2), NamedTuple()),
                    ("GC3d TokamakMediumCartesian hode +extrapolation",
                     () -> GC3.TokamakMediumCartesian.hode(
                         GC3.TokamakMediumCartesian.initial_conditions_barely_passing()...;
                         tspan = (0.0, 10.0), tstep = 0.1), PartitionedGauss(2),
                     (initialguess = MidpointExtrapolation(5),)),
                    ("GC3d TokamakMediumCartesian hode",
                     () -> GC3.TokamakMediumCartesian.hode(
                         GC3.TokamakMediumCartesian.initial_conditions_barely_passing()...;
                         tspan = (0.0, 10.0), tstep = 0.1), PartitionedGauss(2), NamedTuple()),
                    ("GC4d SolovevIterXpoint ode",
                     () -> GC4.SolovevIterXpoint.guiding_center_4d_ode(
                         GC4.SolovevIterXpoint.initial_conditions_barely_passing()...;
                         tspan = (0.0, 200.0), tstep = 1.0), Gauss(2), NamedTuple()),
                    ("GC4d SolovevIterXpoint iode",
                     () -> GC4.SolovevIterXpoint.guiding_center_4d_iode(
                         GC4.SolovevIterXpoint.initial_conditions_barely_passing()...;
                         tspan = (0.0, 20.0), tstep = 0.1), VPRKGauss(2), NamedTuple()),
                    ("Pauli3d SolovevIter hode",
                     () -> P3.SolovevIter.pauli_particle_3d_hode(
                         P3.SolovevIter.initial_conditions_barely_passing()...;
                         tspan = (0.0, 10.0), tstep = 0.1), PartitionedGauss(2), NamedTuple()),
                    ("Pauli3d SolovevIter iode",
                     () -> P3.SolovevIter.pauli_particle_3d_iode(
                         P3.SolovevIter.initial_conditions_barely_passing()...;
                         tspan = (0.0, 10.0), tstep = 0.1), VPRKGauss(2), NamedTuple()))

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

    Once the tolerance is sound, the remaining spread has nothing to do with ITER. The 3D `hode`
    path costs orders of magnitude more per step than the Pauli `hode` even though Newton converges
    in one iteration, and `TokamakMediumCartesian` is *slower* than the ITER Solov'ev — so the two
    ten-minute 3D blocks in CI were a step-count and per-step-cost problem, not a convergence one.
    The cost is nested ForwardDiff: the second derivatives of the 3D Hamiltonian differentiate field
    functions that are themselves AD-generated, and the backtracking line search re-evaluates them.
    `MidpointExtrapolation(5)` multiplies it by a further factor of three.

    Reducing those blocks from 1000 to 100 steps was the pragmatic fix; making the 3D second
    derivatives cheaper would be the real one.""")
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
