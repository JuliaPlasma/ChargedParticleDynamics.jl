
const trajectory_poloidal_linewidth = 3
const trajectory_3d_linewidth = 1

function plot_trajectory_scatter!(fg, ax, R, Z; kwargs...)
    scatter!(ax, R, Z; kwargs...)
    fg, ax
end

function plot_trajectory_poloidal!(fg, ax, R, Z; kwargs...)
    lines!(ax, R, Z; linewidth = trajectory_poloidal_linewidth, kwargs...)
    fg, ax
end

function plot_trajectory_poloidal(R, Z, equ; label = nothing, linewidth = trajectory_poloidal_linewidth, kwargs...)
    xrange = (floor(minimum(R); digits=1), ceil(maximum(R); digits=1))
    yrange = (floor(minimum(Z); digits=1), ceil(maximum(Z); digits=1))
    fg, ax = plot_fieldlines(equ; xrange = xrange, yrange = yrange, kwargs...)
    plot_trajectory_poloidal!(fg, ax, R, Z; label=label, linewidth=linewidth)
    fg, ax
end

function plot_trajectory_poloidal(R, Z; label = nothing, linewidth = trajectory_poloidal_linewidth, kwargs...)
    fg, ax = plot_poloidal_figure(; kwargs...)
    plot_trajectory_poloidal!(fg, ax, R, Z; label=label, linewidth=linewidth)
    fg, ax
end


function plot_trajectory_projection_figure(; xlabel = "X", ylabel = "Y", size = (500, 500), aspect = 1,
        limits = (nothing, nothing), labelsize = 24, ticklabelsize = 18)
    fg = Figure(size = size, figure_padding = (5, 5, 1, 1))
    ax = Axis(fg[1, 1],
        aspect = aspect,
        xlabel = xlabel,
        ylabel = ylabel,
        limits = limits,
        xlabelsize = labelsize, xticklabelsize = ticklabelsize,
        ylabelsize = labelsize, yticklabelsize = ticklabelsize,
    )
    fg, ax
end

function plot_trajectory_projection!(fg, ax, x, y; kwargs...)
    lines!(ax, x, y; linewidth = trajectory_poloidal_linewidth, kwargs...)
    fg, ax
end

function plot_trajectory_projection(x, y; label = nothing, linewidth = trajectory_poloidal_linewidth, kwargs...)
    fg, ax = plot_trajectory_projection_figure(; kwargs...)
    plot_trajectory_projection!(fg, ax, x, y; label = label, linewidth = linewidth)
    fg, ax
end


plot_trajectory_3d_figure(; kwargs...) = plot_3d_figure(; kwargs...)

function plot_trajectory_3d!(fg, ax, X, Y, Z; kwargs...)
    lines!(ax, X, Y, Z; linewidth = trajectory_3d_linewidth, kwargs...)
    fg, ax
end

function plot_trajectory_3d(X, Y, Z; linewidth = trajectory_3d_linewidth, kwargs...)
    fg, ax = plot_trajectory_3d_figure(; kwargs...)
    plot_trajectory_3d!(fg, ax, X, Y, Z; linewidth = linewidth)
end
