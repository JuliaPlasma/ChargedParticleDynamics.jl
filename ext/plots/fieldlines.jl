
function plot_fieldlines_figure(; xrange, yrange, kwargs...)
    plot_poloidal_figure(limits = ((xrange[begin], xrange[end]), (yrange[begin], yrange[end])), kwargs...)
end


"""
    is_axisymmetric_cylindrical(equ)

Whether `equ` is an axisymmetric equilibrium in cylindrical coordinates ``(R, Z, \\varphi)``, which
is the only chart in which [`plot_fieldlines`](@ref) is meaningful.

Decided from the coordinate ranges the field code injects, the same `rangemin`/`rangemax` the
problem constructors use for periodicity: a coordinate is angular when its range is exactly
``[0, 2\\pi)``, and this chart is the one whose *third* coordinate is angular and whose first two
are not. Cartesian equilibria have no angular coordinate; toroidal ``(r, \\theta, \\varphi)`` has
two, so both are excluded.
"""
function is_axisymmetric_cylindrical(equ)
    isdefined(equ, :rangemin) && isdefined(equ, :rangemax) || return false

    lo = equ.rangemin(zeros(3))
    hi = equ.rangemax(zeros(3))

    angular(i) = isfinite(lo[i]) && isfinite(hi[i]) && hi[i] - lo[i] ≈ 2π

    angular(3) && !angular(1) && !angular(2)
end


@doc raw"""
    plot_fieldlines(equ; xrange, yrange, ngrid = (300, 200), levels = 10, kwargs...)

Contour the poloidal flux ``\psi`` of `equ` over the poloidal plane, i.e. the projections of the
magnetic field lines.

The flux is the *covariant* toroidal component of the vector potential, `equ.A₃`, which the field
code injects as ``A_3 = R \, A_{\varphi} = \psi``. Its contours are the flux surfaces, since
``B \cdot \nabla A_3 = 0``. The physical component ``A_{\varphi} = A_3 / R`` is a different
function, and ``B \cdot \nabla (A_3 / R) \neq 0``, so contouring it does not draw field lines — the
same distinction `ElectromagneticFields` 0.8.0 draws in its own equilibrium plots.

!!! warning "Axisymmetric cylindrical equilibria only"
    `equ.A₃` is the flux only where the first two coordinates are ``(R, Z)`` and the third is the
    toroidal angle. That identification is meaningless for a cartesian or a toroidal equilibrium,
    and the contours it draws there are a plausible-looking wrong picture rather than an error.
    `plot_fieldlines` therefore rejects any `equ` that is not axisymmetric cylindrical — see
    [`is_axisymmetric_cylindrical`](@ref).

    This also guards [`plot_trajectory_poloidal`](@ref)`(R, Z, equ)`, which forwards its `equ` here.
    To plot a trajectory over a chart this does not cover, drop the argument and call the
    two-argument `plot_trajectory_poloidal(R, Z)`, which draws the trajectory without the field
    lines.
"""
function plot_fieldlines(equ; xrange, yrange, ngrid=(300, 200), levels=10, kwargs...)
    is_axisymmetric_cylindrical(equ) || throw(ArgumentError(
        "plot_fieldlines contours ψ = A₃ and is only meaningful for an axisymmetric " *
        "equilibrium in cylindrical coordinates (R, Z, φ); $(equ) is not one. Call " *
        "plot_trajectory_poloidal(R, Z) without the equilibrium to plot without field lines."))

    xgrid = LinRange(xrange..., ngrid[1])
    ygrid = LinRange(yrange..., ngrid[2])
    fieldlines = [equ.A₃(0, xgrid[i], ygrid[j], 0.0) for i in eachindex(xgrid), j in eachindex(ygrid)];

    fg, ax = plot_fieldlines_figure(; xrange = xrange, yrange = yrange, kwargs...)
    contour!(ax, xgrid, ygrid, fieldlines, levels = levels, colormap = :reds)
    fg, ax
end
