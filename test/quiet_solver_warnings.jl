#
# Count the nonlinear-solver warnings the test suite provokes, per test file, and keep them out of
# the log. There should be none, and `runtests.jl` asserts it.
#
# `SimpleSolvers` warns when a solve reaches `warn_iterations`, and once per line search that cannot
# find a step satisfying the sufficient-decrease condition. Both are dropped here and counted per
# test file; `runtests.jl` then *asserts* that the count is zero for every file. The suppression is
# what keeps a regression from flooding the log on its way to failing that assertion.
#
# That assertion is the point, and it rests on `SimpleSolvers` 0.10: a converging solve is silent,
# and a solve that cannot converge says so once, naming the residual it achieved. Anything that
# appears is therefore attributable, and a file that starts warning is a real change rather than
# noise — a tolerance that has slipped below a residual floor, or a step at which the integrators
# no longer converge.
#
# `SimpleSolvers` is not a dependency of this environment — it reaches the suite through
# `GeometricIntegrators` — so its messages are identified by the name of the module they originate
# from rather than by importing it.
#
# ## Why the count is zero, and what would make it rise again
#
# Three things keep it there, and each is a place a regression could appear.
#
# `GeometricIntegratorsBase` scales its default tolerance with the stage system,
# `f_abstol = max(8, solversize(method, problem)) * eps(datatype(problem))`. That is what keeps the
# calls which deliberately pass no options at all — `integrate(prob, Gauss(2))` in
# `structure_tests.jl` and the three in `plots_tests.jl` — quiet while still exercising the library
# defaults.
#
# Every equilibrium declares a step at which its own model integrates. Where one does not, the
# solver reaches its iteration cap on a large fraction of steps; `GuidingCenter4d`'s and
# `PauliParticle3d`'s `TokamakSmallCartesian` are the two that are easiest to get wrong, both being
# unusable at `Δt = 500`.
#
# The 4D variational problems are integrated with `SymmetricProjection(VPRKGauss(2))`. Plain
# `VPRKGauss(2)` does not control the parasitic mode of the degenerate discretisation, and a solve
# asked to continue a trajectory that has left the device caps rather than converges. See
# `guiding_center_4d_tests.jl` and `docs/src/findings.md`.
#
# This is deliberately not a way of hiding failures: the tests assert that `integrate` returns a
# solution, which is unaffected by the warnings. Note also that the converse does not hold — an
# orbit can be lost with the solver perfectly quiet — so silence here is a necessary condition for a
# healthy suite, not a sufficient one.
#
module QuietSolverWarnings

using Logging

export quiet_solver_warnings!, suppressed_warning_count, suppressed_warning_counts
export current_test_file!

const QUIET_MODULES = (:SimpleSolvers,)
const COUNTS = Dict{String,Int}()
const CURRENT = Ref("<startup>")

struct QuietLogger{L<:AbstractLogger} <: AbstractLogger
    parent::L
end

# `shouldlog` returns `true` for the quiet modules so that `handle_message` is reached; it drops the
# record there. Counting cannot be done here, because `shouldlog` is consulted once per *call site*
# for loggers that cache it, whereas `handle_message` runs once per message.
function Logging.shouldlog(logger::QuietLogger, level, _module, group, id)
    level < Logging.Error && nameof(_module) ∈ QUIET_MODULES && return true
    Logging.shouldlog(logger.parent, level, _module, group, id)
end

Logging.min_enabled_level(logger::QuietLogger) = Logging.min_enabled_level(logger.parent)
Logging.catch_exceptions(logger::QuietLogger) = Logging.catch_exceptions(logger.parent)

function Logging.handle_message(logger::QuietLogger, level, message, _module, group, id,
                                file, line; kwargs...)
    if level < Logging.Error && nameof(_module) ∈ QUIET_MODULES
        COUNTS[CURRENT[]] = get(COUNTS, CURRENT[], 0) + 1
        return nothing
    end
    Logging.handle_message(logger.parent, level, message, _module, group, id, file, line; kwargs...)
end

quiet_solver_warnings!() = global_logger(QuietLogger(global_logger()))

"""Attribute subsequent suppressed messages to test file `f`."""
current_test_file!(f) = (CURRENT[] = f)

"""Suppressed messages per test file."""
suppressed_warning_counts() = COUNTS

"""Total number of suppressed messages, over all test files."""
suppressed_warning_count() = sum(values(COUNTS); init = 0)

end

using .QuietSolverWarnings
