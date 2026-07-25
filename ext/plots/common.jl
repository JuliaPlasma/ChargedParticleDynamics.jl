
function plot_poloidal_figure(; size=(600, 600), limits=(nothing, nothing), labelsize=24, ticklabelsize=18)
    fg = Figure(size=size, figure_padding=(5, 5, 1, 1))
    ax = Axis(fg[1, 1],
        aspect = 1,
        xlabel = "R",
        ylabel = "Z",
        limits = limits,
        xlabelsize = labelsize, xticklabelsize = ticklabelsize,
        ylabelsize = labelsize, yticklabelsize = ticklabelsize,
    )
    fg, ax
end


function plot_components_figure(labels; size=(600, 800), labelsize=24, ticklabelsize=18)
    fg = Figure(size=size, figure_padding=(5, 5, 1, 1))
    n = length(labels)
    axs = [Axis(fg[i, 1],
        xlabel = i == n ? "t" : "",
        ylabel = labels[i],
        xticklabelsvisible = i == n,
        xlabelsize = labelsize, xticklabelsize = ticklabelsize,
        ylabelsize = labelsize, yticklabelsize = ticklabelsize,
    ) for i in 1:n]
    linkxaxes!(axs...)
    fg, axs
end


function plot_3d_figure(; size=(800, 800), limits=(nothing, nothing, nothing), labelsize=24, ticklabelsize=18)
    fg = Figure(size=size, figure_padding=(1, 50, 1, 1))
    ax = Axis3(fg[1, 1],
        aspect = :equal,
        azimuth = (1.275 + 0.5) * π,
        xlabel = "X",
        ylabel = "Y",
        zlabel = "Z",
        limits = limits,
        xlabelsize = labelsize, xticklabelsize = ticklabelsize,
        ylabelsize = labelsize, yticklabelsize = ticklabelsize,
        zlabelsize = labelsize, zticklabelsize = ticklabelsize,
    )
    fg, ax
end


# Pick at most `nplot` evenly spaced indices out of `indices`, always including both ends.
function select_indices(indices, nplot)
    nplot ≥ length(indices) && return collect(indices)
    nplot ≤ 1 && return [first(indices)]
    unique(round.(Int, range(first(indices), last(indices), length=nplot)))
end
