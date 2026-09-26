# Known issues

Defects found in review and recorded, not fixed. Each entry gives its kind and its evidence.

## KI-1 · The per-file warning counter separates nothing (dead code, docs)

Each test file includes `test/helpers/quiet_solver_warnings.jl` into its own `@safetestset`
module, so each has its own `COUNTS`. The `Dict` keyed by file, `CURRENT`, `current_test_file!`
and the `<startup>` key therefore separate nothing: `current_test_file!` follows
`quiet_solver_warnings!()` directly in every file. The two `@test` lines in each file's footer,
`all(iszero, values(suppressed_warning_counts()))` and `suppressed_warning_count() == 0`, are the
same statement for non-negative counts, and the footer comment ("anything counted before
`current_test_file!` … which the first test cannot see") is false. In `quiet_solver_warnings!`,
`logger isa QuietLogger && return logger` is covered by the `while` loop after it. The 16-line
footer is copied into 9 files; one helper such as `assert_solver_quiet()` with one `Ref{Int}`
counter and one `@test` would remove about 125 lines. Evidence: `test/helpers/quiet_solver_warnings.jl`,
`test/ChargedParticle3d.jl` footer (critics 1a and 1b of part M3).

## KI-2 · `TODO.md` names the old test paths (docs)

`TODO.md:25, 153, 355, 381, 423, 445` name `runtests.jl` asserting per file and the old paths
`test/*_tests.jl` and `test/quiet_solver_warnings.jl`. Evidence:
`grep -rn '_tests\.jl\|runtests' TODO.md`.

## KI-3 · Parts of the 4D guiding centre model are not reached by the tests (missing test)

These mutants survive, on `origin/main` as on the test-suite migration branch:
`dH[4] → dH[5]` and `Ω[4, 4] = 0 → error("mutant")` in
`src/guiding_center_4d/guiding_center_4d_common.jl` (`GuidingCenter4d.jl`); `vcat(ics.X, -ics.u)`
in `src/utils/initial_conditions.jl` and `1.1 * μ * dBdx₁` (`integration/model_agreement.jl`);
`dHdx₄ = 2 * u` (`integration/structure.jl`). So `dH` and `ω(Ω, t, q::FieldPoint)` of the 4D model,
used by `guiding_center_4d_λ`, are not reached. Also a Hamiltonian mutant `μ * B(t, q)` → `2μ * B(t, q)`
in `src/pauli_particle_3d/pauli_particle_3d.jl` survives `integration/model_agreement.jl`.
Evidence: logs `mut_gc4d.log`, `mut_gc4d2.log`, `mut_agree.log`, `mut_agree2.log`,
`mut_structure.log` of critic 1b of part M3.
