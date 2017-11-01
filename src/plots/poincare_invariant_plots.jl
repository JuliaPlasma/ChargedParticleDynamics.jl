

function plot_poincare_error(t::Vector, I::Vector, filename; plot_title=nothing)
    if plot_title == nothing
        plot_title=L"\Delta I"
    end

    xf, yf = get_default_formatter(t)

    I_error = (I - I[1]) / I[1]

    fig = figure(figsize=(6,5))
    subplots_adjust(left=0.15, right=0.93, top=0.94, bottom=0.12)
    plot(t, I_error, ".", markersize=2)
    xticks(linspace(t[1], t[end], 6))
    xlim(t[1], t[end])
    xlabel(L"$t$")
    ylabel(plot_title)
    ax = gca()
    ax[:xaxis][:set_major_formatter](xf)
    ax[:yaxis][:set_major_formatter](yf)
    ax[:yaxis][:set_label_coords](-0.1,0.5)
    savefig(filename)
    close(fig)
end


function plot_poincare_error(t::Vector, I::Vector, K::Vector, Δt, filename; plot_title=nothing)
    if plot_title == nothing
        plot_title=L"\Delta I"
    end

    xf, yf = get_default_formatter(t)

    I_error  = zeros(I)
    I_error .= (I .- I[1]) .- Δt^2 .* (K .- K[1])
    I_error ./= (I[1] - Δt^2 * K[1])

    fig = figure(figsize=(6,5))
    subplots_adjust(left=0.15, right=0.93, top=0.94, bottom=0.12)
    plot(t, I_error, ".", markersize=2)
    xticks(linspace(t[1], t[end], 6))
    xlim(t[1], t[end])
    xlabel(L"$t$")
    ylabel(plot_title)
    ax = gca()
    ax[:xaxis][:set_major_formatter](xf)
    ax[:yaxis][:set_major_formatter](yf)
    ax[:yaxis][:set_label_coords](-0.1,0.5)
    savefig(filename)
    close(fig)
end



function plot_poincare_trajectories(sol::Solution, nplot, filename, dpi=100; xmin=-0.6, xmax=+0.6, ymin=-0.6, ymax=+0.6)
    if nplot > 0
        fig = figure(figsize=(10,8))
        subplots_adjust(left=0.0, right=0.7, top=1.0, bottom=0.0)
        ax  = gca(projection="3d")

        plot3D(sol.q.d[1,1,:], sol.q.d[2,1,:], sol.q.d[3,1,:])

        for i in 1:sol.ni
            if mod(i, div(sol.ni, nplot)) == 0
                plot3D(sol.q.d[1,:,i], sol.q.d[2,:,i], sol.q.d[3,:,i])
            end
        end

        xlabel(L"$x$", fontsize=18, labelpad=10)
        ylabel(L"$y$", fontsize=18, labelpad=10)
        zlabel(L"$z$", fontsize=18, labelpad=10)

        xlim(xmin, xmax)
        ylim(ymin, ymax)

        savefig(filename, dpi=dpi)
        close(fig)
    end
end



function plot_poincare_loop(sol::Solution, nplot, filename, dpi=100; xmin=-0.6, xmax=+0.6, ymin=-0.6, ymax=+0.6)
    if nplot > 0
        fig = figure(figsize=(10,8))
        subplots_adjust(left=0.0, right=0.7, top=1.0, bottom=0.0)
        ax  = gca(projection="3d")

        for i in 1:sol.nt+1
            if mod(i-1, div(sol.nt, nplot)) == 0
                # tpower      = log10(sol.t.t[end])
                # exponent    = @sprintf("%2d", tpower)
                # coefficient = @sprintf("%2.2f", sol.t.t[i] / 10^tpower)
                # tstr = "$coefficient \\times 10^\{ $exponent \}"
                # plot3D(sol.q.d[1,i,:], sol.q.d[2,i,:], sol.q.d[3,i,:], label="\$ t = $tstr \$")
                # coefficient = @sprintf("%4i", sol.t.t[i])
                # tstr = "$coefficient"
                tstr = string(sol.t.t[i])
                plot3D(sol.q.d[1,i,:], sol.q.d[2,i,:], sol.q.d[3,i,:], label="\$ t = $tstr \$")
            end
        end

        xlabel(L"$x$", fontsize=18, labelpad=10)
        ylabel(L"$y$", fontsize=18, labelpad=10)
        zlabel(L"$z$", fontsize=18, labelpad=10)

        xlim(xmin, xmax)
        ylim(ymin, ymax)

        legend(loc=2, bbox_to_anchor=(1.05, 0.8))
        savefig(filename, dpi=dpi)
        close(fig)
    end
end



function plot_poincare_surface(sol::Solution, nplot, filename, dpi=100; xmin=-0.15, xmax=+0.15, ymin=-0.15, ymax=+0.15)
    if nplot > 0
        fig = figure(figsize=(10,8))
        subplots_adjust(left=0.0, right=0.7, top=1.0, bottom=0.0)
        ax  = gca(projection="3d")

        nplot > sol.nt ? nplot = sol.nt : nothing

        for i in 0:sol.nt
            if mod(i, div(sol.nt, nplot)) == 0
                tpower      = log10(sol.t[end])
                exponent    = @sprintf("%2d", floor(Int64, tpower))
                coefficient = @sprintf("%2.2f", sol.t[i] / 10^floor(tpower))
                tstr = "$coefficient \\times 10^\{ $exponent \}"
                plot3D(sol.q.d[1,i+1,:], sol.q.d[2,i+1,:], sol.q.d[3,i+1,:], ".", label="\$ t = $tstr \$")
           end
        end

        xlabel(L"$x$", fontsize=18, labelpad=10)
        ylabel(L"$y$", fontsize=18, labelpad=10)
        zlabel(L"$z$", fontsize=18, labelpad=10)

        xlim(xmin, xmax)
        ylim(ymin, ymax)

        legend(loc=2, bbox_to_anchor=(1.05, 0.8))
        savefig(filename, dpi=dpi)
        close(fig)
    end
end
