module ChargedParticlePlots

    using GeometricIntegrators
    # using Plots
    using PyPlot


    include("plots/plot_common.jl")

    export plot_energy, plot_energy_error, plot_energy_drift,
           plot_toroidal_momentum_error, plot_one_form, plot_momentum, plot_momentum_error,
           plot_lagrange_multiplier

    include("plots/timeseries_plots.jl")

    export plot_trajectory_2d, plot_trajectory_3d

    include("plots/trajectory_plots.jl")

    export plot_poincare_error,
           plot_poincare_loop,
           plot_poincare_surface,
           plot_poincare_trajectories

    include("plots/poincare_invariant_plots.jl")

end
