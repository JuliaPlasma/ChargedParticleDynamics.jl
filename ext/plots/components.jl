
const components_linewidth = 2
const components_markersize_dense = 1
const components_markersize_sparse = 5
const components_markersize_threshold = 1000

components_markersize(n) = n ≤ components_markersize_threshold ? components_markersize_sparse : components_markersize_dense


function plot_components!(fg, axs, t, components; plottype=:lines, kwargs...)
    @assert length(axs) == length(components)
    for (ax, component) in zip(axs, components)
        if plottype == :scatter
            scatter!(ax, t, component; markersize=components_markersize(length(t)), kwargs...)
        else
            lines!(ax, t, component; linewidth=components_linewidth, kwargs...)
        end
    end
    fg, axs
end

function plot_components(t, components; labels, plottype=:lines, kwargs...)
    fg, axs = plot_components_figure(labels; kwargs...)
    plot_components!(fg, axs, t, components; plottype=plottype)
    fg, axs
end


plot_trajectory_cartesian(t, X, Y, Z, u; kwargs...) =
    plot_components(t, (X, Y, Z, u); labels=(L"X", L"Y", L"Z", L"u"), kwargs...)

plot_trajectory_cylindrical(t, R, Z, φ, u; kwargs...) =
    plot_components(t, (R, Z, φ, u); labels=(L"R", L"Z", L"\varphi", L"u"), kwargs...)
