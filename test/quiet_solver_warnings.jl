#
# Suppress the repetitive nonlinear-solver warnings for the duration of the test suite.
#
# The guiding-centre and Pauli-particle models are stiff enough that the Newton solver routinely
# exhausts its iteration budget at the accuracy the tests ask for, and `SimpleSolvers` emits one
# warning per occurrence: a single testset can produce well over thirty thousand of them, which
# buries every other line of a CI log. Only the count is reported, at the end of the run.
#
# This is deliberately *not* a way of hiding failures: the tests assert that `integrate` returns a
# solution, which is unaffected by the warnings. See the comment on the assertions in the
# individual test files for why "no warnings" is not a usable criterion here.
#
module QuietSolverWarnings

using Logging

export quiet_solver_warnings!, suppressed_warning_count

const QUIET_MODULES = (:SimpleSolvers,)
const COUNT = Ref(0)

struct QuietLogger{L<:AbstractLogger} <: AbstractLogger
    parent::L
end

function Logging.shouldlog(logger::QuietLogger, level, _module, group, id)
    if level < Logging.Error && nameof(_module) ∈ QUIET_MODULES
        COUNT[] += 1
        return false
    end
    Logging.shouldlog(logger.parent, level, _module, group, id)
end

Logging.min_enabled_level(logger::QuietLogger) = Logging.min_enabled_level(logger.parent)
Logging.catch_exceptions(logger::QuietLogger) = Logging.catch_exceptions(logger.parent)
Logging.handle_message(logger::QuietLogger, args...; kwargs...) =
    Logging.handle_message(logger.parent, args...; kwargs...)

quiet_solver_warnings!() = global_logger(QuietLogger(global_logger()))

suppressed_warning_count() = COUNT[]

end

using .QuietSolverWarnings
