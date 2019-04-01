
function plot_energy(t::TimeSeries{DT}, h::DataSeries{DT}, filename_prefix; nplot=1, markersize=.5) where {DT}
    @assert t.n == h.nt

    fig = figure(figsize=(8,5))
    subplots_adjust(left=0.16, right=0.95, top=0.95, bottom=0.16)
    plot(t[0:nplot:end], h[1,0:nplot:h.nt], ".", markersize=markersize)
    xlim(t[0], t[end])
    xlabel("t", labelpad=10, fontsize=20)
    ylabel("H(t)", fontsize=20)

    xf, yf = get_default_formatter(t)
    ax = gca()
    ax[:xaxis][:set_major_formatter](xf)
    ax[:yaxis][:set_major_formatter](yf)

    savefig(filename_prefix * "_energy.png")
    close(fig)
end


function plot_energy_error(t::TimeSeries{DT}, h::DataSeries{DT}, filename_prefix; nplot=1, markersize=.5) where {DT}
    @assert t.n == h.nt

    fig = figure(figsize=(8,5))
    subplots_adjust(left=0.16, right=0.95, top=0.95, bottom=0.16)
    plot(t[0:nplot:end], h[2,0:nplot:h.nt], ".", markersize=markersize)
    xlim(t[0], t[end])
    xlabel("t", labelpad=10, fontsize=20)
    ylabel("[H(t) - H(0)] / H(0)", fontsize=20)

    xf, yf = get_default_formatter(t)
    ax = gca()
    ax[:xaxis][:set_major_formatter](xf)
    ax[:yaxis][:set_major_formatter](yf)

    savefig(filename_prefix * "_energy_error.png")
    close(fig)
end


function plot_energy_drift(t::TimeSeries{DT}, drift::DataSeries{DT}, filename_prefix) where {DT}
    fig = figure(figsize=(8,5))
    subplots_adjust(left=0.16, right=0.95, top=0.95, bottom=0.16)
    plot(drift[1,1:drift.nt], drift[2,1:drift.nt], ".", markersize=10)

    ylim1 = minimum(drift[2,1:drift.nt])
    ylim2 = maximum(drift[2,1:drift.nt])
    ylimΔ = ylim2 - ylim1
    ylim1 -= 0.05 * ylimΔ
    ylim2 += 0.05 * ylimΔ

    xlim(t[0], t[end])
    ylim(ylim1, ylim2)

    xlabel("t", labelpad=10, fontsize=20)
    ylabel("Energy Drift", fontsize=20)

    xf, yf = get_default_formatter(t)
    ax = gca()
    ax[:xaxis][:set_major_formatter](xf)
    ax[:yaxis][:set_major_formatter](yf)

    savefig(filename_prefix * "_energy_drift.png")
    close(fig)
end


function plot_toroidal_momentum_error(t::TimeSeries{DT}, p::DataSeries{DT}, filename_prefix; nplot=1, markersize=.5) where {DT}
    @assert t.n == p.nt

    fig = figure(figsize=(8,5))
    subplots_adjust(left=0.16, right=0.95, top=0.95, bottom=0.16)
    plot(t[0:nplot:end], p[2,0:nplot:p.nt], ".", markersize=markersize)
    xlim(t[0], t[end])
    xlabel("t", labelpad=10, fontsize=20)
    ylabel("[P(t) - P(0)] / P(0)", fontsize=20)

    xf, yf = get_default_formatter(t)
    ax = gca()
    ax[:xaxis][:set_major_formatter](xf)
    ax[:yaxis][:set_major_formatter](yf)

    savefig(filename_prefix * "_toroidal_momentum_error.png")
    close(fig)
end


function plot_one_form(t::TimeSeries{DT}, ϑ::DataSeries{DT}, filename_prefix; nplot=1, markersize=.5) where {DT}
    @assert t.n == ϑ.nt

    xf, yf = get_default_formatter(t)

    for i in 1:ϑ.nd
        fig = figure(figsize=(6,5))
        subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.16)
        plot(t[0:nplot:end], ϑ[i,0:nplot:ϑ.nt], ".", markersize=markersize)
        xlim(t[0], t[end])
        if maximum(ϑ[i]) == 0.
            ylim(-2E-16, +2E-16)
        end
        xlabel("\$ t \$", labelpad=10, fontsize=20)
        ylabel("\$ \\vartheta_{" * string(i) * "} (q(t)) \$", labelpad=6, fontsize=20)
        ax = gca()
        ax[:xaxis][:set_major_formatter](xf)
        ax[:yaxis][:set_major_formatter](yf)
        savefig(filename_prefix * "_one_form" * string(i) * ".png")
        close(fig)
    end
end


function plot_momentum(t::TimeSeries{DT}, p::DataSeries{DT}, filename_prefix; nplot=1, markersize=.5) where {DT}
    @assert t.n == p.nt

    xf, yf = get_default_formatter(t)

    for i in 1:p.nd
        fig = figure(figsize=(6,5))
        subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.16)
        plot(t[0:nplot:end], p[i,0:nplot:p.nt], ".", markersize=markersize)
        xlim(t[0], t[end])
        if maximum(p[i]) == 0.
            ylim(-2E-16, +2E-16)
        end
        xlabel("\$ t \$", labelpad=10, fontsize=20)
        ylabel("\$ p_{" * string(i) * "} (t) \$", labelpad=6, fontsize=20)
        ax = gca()
        ax[:xaxis][:set_major_formatter](xf)
        ax[:yaxis][:set_major_formatter](yf)
        savefig(filename_prefix * "_momentum" * string(i) * ".png")
        close(fig)
    end
end


function plot_momentum_error(t::TimeSeries{DT}, e::DataSeries{DT}, filename_prefix; nplot=1, markersize=.5) where {DT}
    @assert t.n == e.nt

    xf, yf = get_default_formatter(t)

    for i in 1:e.nd
        fig = figure(figsize=(6,5))
        subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.16)
        plot(t[0:nplot:end], e[i,0:nplot:e.nt], ".", markersize=markersize)
        xlim(t[0], t[end])
        if maximum(abs.(e[i])) == 0.
            ylim(-2E-16, +2E-16)
        end
        xlabel("\$ t \$", labelpad=10, fontsize=20)
        ylabel("\$ p_{" * string(i) * "} (t) - \\vartheta_{" * string(i) * "} (q(t)) \$", labelpad=6, fontsize=20)
        ax = gca()
        ax[:xaxis][:set_major_formatter](xf)
        ax[:yaxis][:set_major_formatter](yf)
        savefig(filename_prefix * "_momentum_error" * string(i) * ".png")
        close(fig)
    end
end


function plot_lagrange_multiplier(t::TimeSeries{DT}, λ::DataSeries{DT}, filename_prefix; nplot=1, markersize=.5) where {DT}
    @assert t.n == λ.nt

    xf, yf = get_default_formatter(t)

    for i in 1:λ.nd
        fig = figure(figsize=(6,5))
        subplots_adjust(left=0.16, right=0.95, top=0.95, bottom=0.16)
        plot(t[0:nplot:end], λ[i,0:nplot:λ.nt], ".", markersize=markersize)
        xlim(t[0], t[end])
        if maximum(abs.(λ[i,:])) == 0.
            ylim(-2E-16, +2E-16)
        end
        xlabel("\$t\$", labelpad=10, fontsize=20)
        ylabel("\$\\lambda_{" * string(i) * "} (t)\$", labelpad=6, fontsize=20)
        ax = gca()
        ax[:xaxis][:set_major_formatter](xf)
        ax[:yaxis][:set_major_formatter](yf)
        savefig(filename_prefix * "_lambda" * string(i) * ".png")
        close(fig)
    end
end
