# Changelog

All notable changes to ChargedParticleDynamics.jl are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [Unreleased]

### Added

- **Guiding centre dynamics in 3D.** A new `GuidingCenter3d` module implements the guiding
  centre equations as a constrained canonical Hamiltonian system in the position `r` and its
  conjugate momentum `p = A(r) + u b(r)`, with the parallel velocity recovered as a dependent
  variable. The two independent components of `v × b = 0` are imposed as constraints `g₁` and
  `g₂`, whose Lagrange multipliers `λ₁` and `λ₂` are determined algebraically. Constructed by
  `hode`. Available for thirteen equilibria: dipole, quadratic potentials, symmetric quadratic,
  theta pinch, symmetric Solov'ev, ITER-like Solov'ev with and without X-point, the small tokamak
  in cartesian, cylindrical and toroidal coordinates, the medium tokamak in cartesian and
  cylindrical coordinates, and the ITER-size tokamak in cylindrical coordinates.
- **Canonicalised 3D guiding centre formulation.** `hode_canonical` absorbs the multipliers into
  the Hamiltonian, `H̄ = H + λ₁ g₁ + λ₂ g₂`, which agrees with `H` on the constraint manifold but
  yields a genuinely canonical system. It requires the second derivatives of the magnetic field
  and lives in `src/guiding_center_3d/guiding_center_3d_canonical.jl`, included from all thirteen
  equilibria. `hamiltonian_canonical` is now exported alongside it.
- **Makie plotting.** A `ChargedParticlePlots` package extension, loaded when `Makie` and
  `LaTeXStrings` are available, provides `plot_fieldlines`, `plot_invariant`,
  `plot_invariant_error`, `plot_components`, `plot_trajectory_cartesian`,
  `plot_trajectory_cylindrical`, `plot_trajectory_poloidal`, `plot_trajectory_projection`,
  `plot_trajectory_3d`, `plot_poincare_invariant_error`, `plot_poincare_loop`,
  `plot_poincare_surface` and `plot_poincare_trajectories`, together with their mutating
  variants. These are plain functions returning `(figure, axis)` and taking plain coordinate
  vectors, rather than `Plots.jl` recipes taking solution objects.
- `cartesian_solution` converts a solution to cartesian coordinates, returning the time series
  and the `R`, `X`, `Y` and `Z` components.
- Documentation page for the 3D guiding centre model, and a worked example integrating the
  charged particle, Pauli particle and 3D and 4D guiding centre models in an ITER-like
  equilibrium in cylindrical coordinates.
- `initial_conditions_default` for the small tokamak in cartesian coordinates.
- Example script for the Pauli particle.

### Changed

- Updated to the current `GeometricEquations` problem interface and the current
  `GeometricSolutions` diagnostics interface. `compute_energy`, `compute_toroidal_momentum` and
  their error variants now take the problem parameters, following the new `compute_invariant`
  and `compute_invariant_error` signatures; the convenience methods taking a `GeometricSolution`
  obtain them from the problem. `compute_momentum_error`, `compute_one_form` and
  `compute_error_drift` are re-exported from `GeometricSolutions`.
- Initial conditions and coordinate helpers moved from `src/initialization` to `src/utils`.
- Minimum Julia version raised from 1.6 to 1.10.
- CI now tests Julia 1.10, 1.11, 1.12 and nightly.
- The documentation is built from its own environment (`docs/Project.toml`) rather than the
  package environment.
- Various initial conditions, time steps and integration intervals of the 3D guiding centre
  equilibria were adjusted.

### Fixed

- `compute_energy` and friends threw a `MethodError` for every model, because the imports had
  been switched from `GeometricProblems.Diagnostics` to `GeometricSolutions` without updating
  the call sites for the added `parameters` argument.
- `compute_momentum_error` in the 4D guiding centre diagnostics still called into
  `GeometricProblems`, and `compute_one_form` and `compute_error_drift` were exported without
  being defined.
- `convert_coordinates_RZphi_to_xyz` used the removed `SDataSeries`; it is replaced by
  `cartesian_solution`.
- `toroidal_momentum` referred to the undefined `ϑ3` instead of `ϑ₃` in the symmetric Solov'ev
  equilibrium and in the small tokamak in cartesian, cylindrical and toroidal coordinates, in
  both the 3D and the 4D guiding centre models.
- The medium tokamak in cartesian coordinates exported `toroidal_momentum` without defining it.
- `plot_trajectory_3d_figure` declared positional instead of keyword arguments, so any keyword
  passed to `plot_trajectory_3d` raised a `MethodError`.
- `plot_invariant_error` divided by `log10(0)` for an identically vanishing error series, and
  had an unbalanced bracket in its axis label.
- The test suite aborted before the 4D guiding centre, gyrokinetic and Pauli particle tests,
  which imported `MidpointExtrapolation` from `GeometricIntegratorsBase` — no longer part of the
  test environment. Once reachable again, those tests revealed solvers that did not converge at
  the default time step of some equilibria, and a degenerate Hermite extrapolation of the initial
  guess for the Pauli theta pinch, whose momentum is constant along the initial trajectory.
- The example block plotting the ITER magnetic field was missing its `@` and was never executed;
  the initial conditions figure of the initialization page is reactivated using Makie.
- Numerous fixes to the Pauli particle Solov'ev equilibria, the 3D guiding centre equations and
  the scripts and examples.

### Removed

- The `Plots`, `RecipesBase`, `GeometricProblems`, `GeometricIntegrators`,
  `GeometricIntegratorsBase` and `ForwardDiff` dependencies, together with the `Plots.jl`
  recipes and the PyPlot-based Poincaré invariant plots they supported. `GeometricIntegrators`
  remains a test dependency.
- The unused `Documenter`, `EulerLagrange` and `Printf` dependencies, and the `SimpleSolvers`
  test dependency.


## [0.1.0] - 2024-12-13

Initial release.


[Unreleased]: https://github.com/JuliaPlasma/ChargedParticleDynamics.jl/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/JuliaPlasma/ChargedParticleDynamics.jl/releases/tag/v0.1.0
