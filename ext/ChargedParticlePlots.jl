module ChargedParticlePlots

using ChargedParticleDynamics
using GeometricSolutions
using LaTeXStrings
using Makie

import ChargedParticleDynamics: plot_invariant, plot_invariant_error
import ChargedParticleDynamics: plot_fieldlines_figure, plot_fieldlines, is_axisymmetric_cylindrical
import ChargedParticleDynamics: plot_poloidal_figure, plot_trajectory_poloidal, plot_trajectory_poloidal!, plot_trajectory_scatter!
import ChargedParticleDynamics: plot_components_figure, plot_components, plot_components!
import ChargedParticleDynamics: plot_trajectory_cartesian, plot_trajectory_cylindrical
import ChargedParticleDynamics: plot_trajectory_projection_figure, plot_trajectory_projection, plot_trajectory_projection!
import ChargedParticleDynamics: plot_trajectory_3d_figure, plot_trajectory_3d, plot_trajectory_3d!
import ChargedParticleDynamics: plot_poincare_invariant_error
import ChargedParticleDynamics: plot_poincare_loop, plot_poincare_surface, plot_poincare_trajectories


include("plots/common.jl")
include("plots/fieldlines.jl")
include("plots/invariants.jl")
include("plots/components.jl")
include("plots/trajectories.jl")
include("plots/poincare_invariants.jl")

end
