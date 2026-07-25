
function plot_invariant(t, invariant; label, size=(800, 400), labelsize=24, ticklabelsize=18, padding_left=5, padding_right=5)
    fg = Figure(size=size, figure_padding=(padding_left, padding_right, 5, 5))
    ax = Axis(fg[1, 1],
        xlabel="t",
        ylabel="$(label)(t)",
        limits=((t[begin], t[end]), nothing),
        xlabelsize=labelsize, xticklabelsize=ticklabelsize,
        ylabelsize=labelsize, yticklabelsize=ticklabelsize,
    )
    lines!(ax, t, invariant)
    fg
end


# Plot a precomputed relative error series `e` for the quantity named `label`.
# The vertical axis is rescaled by the leading power of ten of the error, which is
# annotated above the axis; a vanishing error is plotted unscaled.
function plot_error(t, e; label, size=(800, 400), labelsize=24, ticklabelsize=18, padding_left=5, padding_right=5)
    fg = Figure(size=size, figure_padding=(padding_left, padding_right, 5, 5))
    ax = Axis(fg[1, 1],
        xlabel="t",
        ylabel="[$(label)(t) - $(label)(0)] / $(label)(0)",
        limits=((t[begin], t[end]), nothing),
        xlabelsize=labelsize, xticklabelsize=ticklabelsize,
        ylabelsize=labelsize, yticklabelsize=ticklabelsize,
    )
    emax = maximum(abs, e)
    if iszero(emax)
        lines!(ax, t, e)
    else
        expnt = floor(Int, log10(emax))
        Label(fg[1, 1, Top()], fontsize=18, halign=:left, L"\times 10^{%$(expnt)}")
        lines!(ax, t, e ./ 10.0^expnt)
    end
    fg
end


function plot_invariant_error(t, invariant; label, kwargs...)
    plot_error(t, compute_relative_error(invariant); label=label, kwargs...)
end
