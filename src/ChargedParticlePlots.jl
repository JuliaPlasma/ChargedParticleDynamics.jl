module ChargedParticlePlots

    # using GeometricIntegrators
    using Plots


    include("plots/cartesian.jl")
    include("plots/cylindrical.jl")


    # export plot_poincare_error,
    #        plot_poincare_loop,
    #        plot_poincare_surface,
    #        plot_poincare_trajectories
    #
    # include("plots/poincare_invariant_plots.jl")

end
