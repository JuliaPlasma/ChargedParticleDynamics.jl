
function plot_fieldlines_figure(; xrange, yrange, kwargs...)
    plot_poloidal_figure(limits = ((xrange[begin], xrange[end]), (yrange[begin], yrange[end])), kwargs...)
end

function plot_fieldlines(equ; xrange, yrange, ngrid=(300, 200), levels=10, kwargs...)
    xgrid = LinRange(xrange..., ngrid[1])
    ygrid = LinRange(yrange..., ngrid[2])
    fieldlines = [equ.A₃(0, xgrid[i], ygrid[j], 0.0) / xgrid[i] for i in eachindex(xgrid), j in eachindex(ygrid)];

    fg, ax = plot_fieldlines_figure(; xrange = xrange, yrange = yrange, kwargs...)
    contour!(ax, xgrid, ygrid, fieldlines, levels = levels, colormap = :reds)
    fg, ax
end
