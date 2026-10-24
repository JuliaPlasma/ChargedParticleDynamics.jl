using Aqua
using ChargedParticleDynamics
using Test

# Package-level quality assurance: method ambiguities, unbound type parameters, undefined exports,
# the agreement between `Project.toml` and `test/Project.toml`, stale dependencies, `[compat]`
# bounds, type piracy and persistent tasks.
Aqua.test_all(ChargedParticleDynamics)
