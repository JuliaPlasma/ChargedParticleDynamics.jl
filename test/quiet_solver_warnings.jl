#
# Suppress the repetitive nonlinear-solver warnings for the duration of the test suite.
#
# `SimpleSolvers` emits one warning per solve that exhausts its iteration budget, and which orbits
# trip it depends on the floating-point details of the platform. Only the count is reported, at the
# end of the run.
#
# The count is worth reading, but it is *not* a tripwire and nothing asserts on it. It used to sit
# above fifty thousand because the tests requested an `f_abstol` below the round-off floor of the
# ITER-scale residuals, so Newton could not converge and ran to its iteration limit on nearly half
# of all steps; see the comment on `options` in `guiding_center_3d_tests.jl`. With the tolerances
# now in use the blocks that pass those options contribute none of them, so a sharp rise means a
# tolerance has slipped below a residual floor again — but that has to be noticed by reading the
# log, not by a failing test.
#
# Asserting a bound was considered and rejected: which orbits exhaust the budget depends on the
# floating-point details of the platform, which is the whole reason the warnings are filtered here,
# so any threshold tight enough to catch a regression would also flake across platforms.
#
# The residue that remains comes from the calls in `structure_tests.jl` and `plots_tests.jl` that
# deliberately integrate without passing `options` in order to exercise the library defaults — and
# those defaults ask for `f_abstol = 8eps()`, which is itself below the ITER-scale floor.
#
# Every call that *does* pass options now asks for `1E-12`. The 3D guiding centre diagnostics block
# of `structure_tests.jl` used to restate `8eps()` by hand while still passing options, which is the
# worst of both: it configured the solver and then asked it for a residual below the round-off floor,
# and four steps of the Solov'ev X-point ran to the iteration cap chasing it. See the comment on
# `options` there.
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
