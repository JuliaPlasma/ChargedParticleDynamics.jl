# Known issues

### K1 · The per-file warning counter separates nothing.

- **location:** `test/helpers/quiet_solver_warnings.jl`
- **evidence:** Each test file includes `test/helpers/quiet_solver_warnings.jl` into its own
  `@safetestset` module, so each has its own `COUNTS`. The `Dict` keyed by file, `CURRENT`,
  `current_test_file!` and the `<startup>` key therefore separate nothing: `current_test_file!`
  follows `quiet_solver_warnings!()` directly in every file. The two `@test` lines in each file's
  footer, `all(iszero, values(suppressed_warning_counts()))` and `suppressed_warning_count() == 0`,
  are the same statement for non-negative counts. In `quiet_solver_warnings!`,
  `logger isa QuietLogger && return logger` is covered by the `while` loop after it. The 15-line
  footer is copied into 9 files; one helper such as `assert_solver_quiet()` with one `Ref{Int}`
  counter and one `@test` would remove about 115 lines. See the footer of
  `test/ChargedParticle3d.jl`.
- **kind:** dead code
- **found:** 2026-09-26

### K2 · `TODO.md` names the old test paths.

- **location:** `TODO.md:25`
- **evidence:** `TODO.md:25, 153, 355, 381, 423, 445` name `runtests.jl` asserting per file and the
  old paths `test/*_tests.jl` and `test/quiet_solver_warnings.jl`.
  `grep -rn '_tests\.jl\|runtests' TODO.md`.
- **kind:** docs
- **found:** 2026-09-26

### K3 · Parts of the 4D guiding centre model are not reached by the tests.

- **location:** `src/guiding_center_4d/guiding_center_4d_common.jl`
- **evidence:** These mutants survive the test suite: `dH[4] → dH[5]` and
  `Ω[4, 4] = 0 → error("mutant")` in `src/guiding_center_4d/guiding_center_4d_common.jl`
  (`GuidingCenter4d.jl`); `vcat(ics.X, -ics.u)` in `src/utils/initial_conditions.jl` and
  `1.1 * μ * dBdx₁` (`integration/model_agreement.jl`); `dHdx₄ = 2 * u` (`integration/structure.jl`).
  So `dH` and `ω(Ω, t, q::FieldPoint)` of the 4D model, used by `guiding_center_4d_λ`, are not
  reached. Also a Hamiltonian mutant `μ * B(t, q)` → `2μ * B(t, q)` in
  `src/pauli_particle_3d/pauli_particle_3d.jl` survives `integration/model_agreement.jl`. Each
  mutant, applied alone to `src/`, leaves the test file named after it passing.
- **kind:** missing test
- **found:** 2026-09-26, in the test suite before and after its reorganisation into `core` and
  `slow`

### K4 · Eight comments in `src/` name test files by paths that do not exist.

- **location:** `src/guiding_center_3d/guiding_center_3d_compact.jl:66`
- **evidence:** `grep -rn '_tests\.jl' src` lists
  `src/guiding_center_3d/guiding_center_3d_compact.jl:66`,
  `src/guiding_center_3d/guiding_center_3d_canonical.jl:345`,
  `src/guiding_center_3d/tokamak_medium_cartesian.jl:21`,
  `src/guiding_center_3d/tokamak_small_cartesian.jl:19`,
  `src/guiding_center_3d/tokamak_small_cylindrical.jl:21`,
  `src/guiding_center_3d/tokamak_small_toroidal.jl:21`, `src/gyro_kinetics_4d/gc_common.jl:247` and
  `src/gyro_kinetics_4d/gc_tokamak_small_toroidal.jl:41`. They name `test/structure_tests.jl`,
  `test/guiding_center_3d_tests.jl` and `test/gyro_kinetics_4d_tests.jl`, which are
  `test/integration/structure.jl`, `test/GuidingCenter3d.jl` and `test/GyroKinetics4d.jl`.
- **kind:** docs
- **found:** 2026-09-26
