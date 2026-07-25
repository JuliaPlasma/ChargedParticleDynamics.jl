#
# Suppress the repetitive nonlinear-solver warnings for the duration of the test suite.
#
# `SimpleSolvers` emits one warning per solve that exhausts its iteration budget, and which orbits
# trip it depends on the floating-point details of the platform. Only the count is reported, at the
# end of the run.
#
# That count is also a tripwire. It used to sit above fifty thousand because the tests requested an
# `f_abstol` below the round-off floor of the ITER-scale residuals, so Newton could not converge and
# ran to its iteration limit on nearly half of all steps; see the comment on `options` in
# `guiding_center_3d_tests.jl`. With the tolerances now in use the blocks that pass those options
# contribute none of them, so a sharp rise means a tolerance has slipped below a residual floor
# again.
#
# The residue that remains comes from `structure_tests.jl` and `plots_tests.jl`, which deliberately
# integrate without passing `options` in order to exercise the library defaults — and those defaults
# ask for `f_abstol = 8eps()`, which is itself below the ITER-scale floor.
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
