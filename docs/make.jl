using Documenter
using ChargedParticleDynamics

makedocs(
    sitename = "ChargedParticleDynamics.jl",
    warnonly = Documenter.except(:autodocs_block, :cross_references, :docs_block, :doctest, :eval_block, :example_block, :footnote, :linkcheck_remotes, :linkcheck, :meta_block, :parse_error, :setup_block),
    format = Documenter.HTML(
                prettyurls = get(ENV, "CI", nothing) == "true",
                assets = [asset("assets/style.css", class=:css, islocal=true)]),
    pages = ["Overview"                      => "index.md",
             "Normalization"                 => "normalization.md",
             "Initialization"                => "initialization.md",
             "Charged Particles in 3D"       => "charged_particle_3d.md",
             "Pauli Particles in 3D"         => "pauli_particle_3d.md",
             "Guiding Center Dynamics in 4D" => "guiding_center_4d.md",
            ]
)

deploydocs(
    repo   = "github.com/JuliaPlasma/ChargedParticleDynamics.jl",
    devurl = "latest",
    devbranch = "main",
    push_preview = true,
)
