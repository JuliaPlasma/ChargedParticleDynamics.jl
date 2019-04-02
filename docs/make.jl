using Documenter
using ChargedParticleDynamics

makedocs(
    sitename = "ChargedParticleDynamics.jl",
    format = Documenter.HTML(),
    pages = ["Overview" => "index.md",
             "Charged Particles in 3D"       => "charged_particle_3d.md",
             "Guiding Center Dynamics in 4D" => "guiding_center_4d.md",
            ]
)
