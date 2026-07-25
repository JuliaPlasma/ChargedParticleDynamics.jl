
const poincare_loop_linewidth = 2
const poincare_surface_markersize = 3


plot_poincare_invariant_error(t, I; label = "I", kwargs...) = plot_error(t, compute_relative_error(I); label = label, kwargs...)

# Error of the Poincaré invariant `I` corrected by its `O(Δt²)` contribution `K`:
#     [(I(t) - I(0)) - Δt² (K(t) - K(0))] / (I(0) - Δt² K(0))
function plot_poincare_invariant_error(t, I, K, Δt; label = "I", kwargs...)
    e = ((I .- I[begin]) .- Δt^2 .* (K .- K[begin])) ./ (I[begin] - Δt^2 * K[begin])
    plot_error(t, e; label = label, kwargs...)
end


# `X`, `Y` and `Z` hold one vector of coordinates per curve: for a loop or a surface that is one
# entry per plotted time, for a bundle of trajectories one entry per trajectory.
function plot_poincare_curves!(fg, ax, X, Y, Z; labels = nothing, plottype = :lines, kwargs...)
    @assert length(X) == length(Y) == length(Z)
    for i in eachindex(X)
        label = isnothing(labels) ? nothing : labels[i]
        if plottype == :scatter
            scatter!(ax, X[i], Y[i], Z[i]; markersize = poincare_surface_markersize, label = label, kwargs...)
        else
            lines!(ax, X[i], Y[i], Z[i]; linewidth = poincare_loop_linewidth, label = label, kwargs...)
        end
    end
    fg, ax
end


function plot_poincare_loop(t, X, Y, Z; nplot = 5, plottype = :lines, sigdigits = 4, kwargs...)
    steps = select_indices(eachindex(X), nplot)
    labels = ["t = $(round(t[n]; sigdigits = sigdigits))" for n in steps]
    fg, ax = plot_3d_figure(; kwargs...)
    plot_poincare_curves!(fg, ax, X[steps], Y[steps], Z[steps]; labels = labels, plottype = plottype)
    Legend(fg[1, 2], ax)
    fg, ax
end

plot_poincare_surface(t, X, Y, Z; kwargs...) = plot_poincare_loop(t, X, Y, Z; plottype = :scatter, kwargs...)


function plot_poincare_trajectories(X, Y, Z; nplot = 5, kwargs...)
    trajectories = select_indices(eachindex(X), nplot)
    fg, ax = plot_3d_figure(; kwargs...)
    plot_poincare_curves!(fg, ax, X[trajectories], Y[trajectories], Z[trajectories])
    fg, ax
end
