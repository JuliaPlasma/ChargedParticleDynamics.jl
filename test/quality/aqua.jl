using Aqua
using ChargedParticleDynamics
using Test

# Package-level quality assurance: method ambiguities, unbound type parameters, undefined exports,
# the agreement between `Project.toml` and `test/Project.toml`, stale dependencies, `[compat]`
# bounds, type piracy and persistent tasks.
#
# On Julia 1.11 the persistent-tasks check precompiles a wrapper package whose manifest mirrors the
# test environment, and with CairoMakie in that environment the wrapper also builds the Makie
# extensions `ChargedParticlePlots` and `PoincareInvariantsMakieExt`. The check took 53 s and 56 s
# in two runs, more than Aqua's default `tmax` of 30 s. A real persistent task blocks forever, so a
# larger `tmax` hides none.
Aqua.test_all(ChargedParticleDynamics; persistent_tasks = (; tmax = 300))
