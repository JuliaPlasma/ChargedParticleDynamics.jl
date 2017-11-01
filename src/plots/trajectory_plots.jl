
function plot_trajectory_2d(q::DataSeries, filename_prefix; markersize=.5, nplot=1, XLIM=nothing, YLIM=nothing)
    fig = figure(figsize=(5,5))
    subplots_adjust(left=0.25, right=0.95, top=0.96, bottom=0.16)
    plot(q[1,0:nplot:q.nt], q[2,0:nplot:q.nt], ".", markersize=markersize)

    XLIM != nothing ? xlim(XLIM) : nothing
    YLIM != nothing ? ylim(YLIM) : nothing

    xlabel(L"$q_{1} (t)$", labelpad=10, fontsize=20)
    ylabel(L"$q_{2} (t)$", labelpad=10, fontsize=20)
    savefig(filename_prefix * "_trajectory_2d.png")
    close(fig)
end


function plot_trajectory_3d(q::DataSeries, filename_prefix; markersize=.5, ntmax=0)
    ntmax == 0 ? ntmax = q.nt : nothing

    fig = figure(figsize=(10,8))
    subplots_adjust(left=0.05, right=0.95, top=0.98, bottom=0.05)
    plot3D(q[1,0:ntmax], q[2,0:ntmax], q[3,0:ntmax], ".-", markersize=markersize)

    xlabel(L"$q_{1} (t)$", labelpad=10, fontsize=20)
    ylabel(L"$q_{2} (t)$", labelpad=10, fontsize=20)
    zlabel(L"$q_{3} (t)$", labelpad=10, fontsize=20)

    savefig(filename_prefix * "_trajectory_3d.png")
    close(fig)
end
