using ParallelTestRunner
using ChargedParticleDynamics

const TESTS = [
    "core" => [
        "Structure" => "integration/structure.jl",
        "Charged particle 3D" => "ChargedParticle3d.jl",
        "Guiding centre 4D" => "GuidingCenter4d.jl",
        "Gyrokinetics 4D" => "GyroKinetics4d.jl",
        "Pauli particle 3D" => "PauliParticle3d.jl",
        "Model agreement" => "integration/model_agreement.jl"
    ],
    "quality" => [
        "Aqua" => "quality/aqua.jl",
        "Plots" => "Plots.jl"
    ],
    "slow" => [
        "Guiding centre 3D" => "GuidingCenter3d.jl",
        "Poincaré invariants" => "integration/poincare_invariants.jl"
    ]
]

const SELECTED = filter(!startswith("--"), ARGS)
const GROUPS = isempty(SELECTED) ? ["core", "quality", "slow"] : SELECTED
const UNKNOWN = setdiff(SELECTED, first.(TESTS))
isempty(UNKNOWN) ||
    error("no test group ", join(UNKNOWN, ", "), " in runtests.jl; its groups are ",
        join(first.(TESTS), ", "))

runtests(ChargedParticleDynamics, filter(startswith("--"), ARGS);
    testsuite = Dict(label => :(include($(joinpath(@__DIR__, path))))
    for (group, files) in TESTS if group in GROUPS for (label, path) in files),
    serial = [label
              for (group, files) in TESTS if group in ("metal", "cuda")
              for (label, _) in files],
    serial_position = :after,
    retries = 0,
    exeflags = ["--threads=$(Threads.nthreads())"],
    max_default_jobs = max(1, Sys.CPU_THREADS ÷ Threads.nthreads()))
