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

### K5 · The `scripts/` environment is bounded to releases that its scripts have not run against.

- **location:** `scripts/Project.toml`
- **evidence:** The floor raise to GeometricBase 0.15 sets `ElectromagneticFields = "0.9.1"`,
  `GeometricIntegrators = "0.18.6"`, `GeometricIntegratorsBase = "0.6.9"` and
  `SimpleSolvers = "0.14.1"`. The old bounds were `"0.9"`, `"0.17"`, `"0.5.1"` and `"0.10"`, so the
  skipped versions are GeometricIntegrators 0.18.0–0.18.5, GeometricIntegratorsBase 0.6.0–0.6.8
  and SimpleSolvers 0.11–0.14.0. No script source was changed and none was run. Known risks from
  the skipped CHANGELOGs: SimpleSolvers 0.11 removes `Backtracking`'s `α₀`; SimpleSolvers 0.12
  stops the line-search warnings inside `solver_step!`, so the warning string matches of
  `scripts/study_solver_tolerances.jl` and `scripts/study_guiding_center_3d_conditioning.jl` can
  silently match nothing; SimpleSolvers 0.13 changes the default linear solver for LAPACK element
  types. GeometricIntegrators 0.18.0 renames `SymplecticEulerA` and `SymplecticEulerB` to
  `SymplecticEulerARK` and `SymplecticEulerBRK`; `scripts/pauli_particle.jl:18–19` names the old
  types in two commented-out lines, which break when they are uncommented.
- **follow-up:** resolve the environment, run each script, and fix the callers in a later pull
  request.
- **kind:** not verified
- **found:** 2026-10-02

### K7 · `integration/model_agreement.jl` does not finish on Julia 1.14.0-beta1 and on nightly.

- **location:** `test/integration/model_agreement.jl`
- **evidence:** the CI jobs that run the whole suite stop at the 2 h limit after
  `Pauli particle 3D`, with no output from "Model agreement" in between. Julia pre: on
  1.14.0-alpha2 (`main`, job 125935516940, 2026-11-03) "Model agreement" passes, 279 tests in
  52.7 s; on 1.14.0-beta1 (#62, job 132086057883) the last line is `Pauli particle 3D | 36 36 55.9s`
  at 18:46:44, then "The operation was canceled." at 20:35:02. Nightly: on 1.15.0-DEV.420 (`main`,
  job 125935517025) the last line is `Pauli particle 3D` at 13:37:38, then the cancel at 15:18:09;
  #62's nightly job 132086057903 stops the same way. The test file is unchanged between the runs.
  The pre and nightly jobs are no required checks. Whether the Julia version or a dependency
  version causes the hang is not checked.
- **kind:** upstream
- **found:** #62

### K8 · The second field derivatives of the 3D guiding centre are wrong on AMD Zen 4 and Zen 5 CPUs with Julia 1.12.7.

- **location:** ElectromagneticFields `src/analytic/field_functions.jl:265`
- **evidence:** ElectromagneticFields derives `DDb` symbolically
  (`src/analytic/analytic_field.jl:510`), compiles it with `Symbolics.build_function(…; cse = true)`
  (`field_functions.jl:265`) and stores the compiled functions in its precompile image through the
  `@compile_workload` of `src/ElectromagneticFields.jl`. This package only reads the values
  (`src/guiding_center_3d/guiding_center_3d_constraints.jl:39-41`). In job 132424694448,
  "Structure" fails 36 of 2025 tests: 7 in "the field accessors are the derivatives they name"
  (`d²bᵢdxⱼdxₖ` against a finite difference, `test/integration/structure.jl:345`), 23 in
  "constraint derivatives match finite differences" (`d²gᵏdqₗdqₘ`, lines 462 and 465: values such
  as `3.3759` against `-0.1065`, `NaN` and `-Inf`, and the symmetry in `l` and `m` broken) and 6 in
  "multiplier derivatives match finite differences". The failure follows the CPU of the runner,
  which `ubuntu-latest` assigns at random, and not the versions. The diagnostic run 43696308235 of
  #64 repeats the `Julia min - ubuntu-latest` and `Downgrade - ubuntu-latest` core jobs six times
  each, with no depot cache, and then probes `DDb♭` against a central difference of `Db♭`. All 8
  jobs on AMD EPYC 9V74 (`znver4`) and 9V45 (`znver5`) fail "Structure"; all 4 jobs on Intel Xeon
  8573C (`sapphirerapids`), Intel Xeon 6975P-C (`graniterapids`) and AMD EPYC 7763 (`znver3`) pass.
  In job 132757499116 (`znver5`) the probe finds 14 wrong entries, all of the Cartesian tokamak field
  of `TokamakMediumCartesian`, such as `81.98` against `0.0076`, `NaN` against `0.0` and `-Inf`
  against `-0.171`; in job 132757499018 (`sapphirerapids`) it finds none. The generated code of
  that `DDb♭` has the same hash in both jobs, so the same Julia code gives wrong values on Zen 4
  and Zen 5 only: the defect is in the machine code that Julia 1.12.7 compiles for these CPUs,
  not in the Symbolics output. The macOS and Windows jobs, `Julia 1` on Ubuntu, and nine local
  runs on aarch64 macOS give correct values. Which compiler pass fails is not established.
- **kind:** upstream
- **found:** #62
