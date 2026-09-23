using ElectromagneticFields: periodic, rangemin, rangemax

@doc raw"""
    periodic_domain(field, T, n)

The periodic domain of an `n`-component state whose first three components are the coordinates of
`field`'s chart, as the `(xmin, xmax)` tuple that `GeometricEquations` takes as `periodicity`.

A component is periodic where `periodic(field)` says so, and wraps on the chart's own
`rangemin`/`rangemax`. Every other component, including every one past the third, gets
``(-\infty, +\infty)``.

`periodic(field)` answers one `Bool` per coordinate and is not itself a `periodicity`: passed in its
place, a three-element vector destructures to `(false, false)` without any error. A bounded range
does not imply periodicity either, which is why the bounds are taken only where `periodic` says so.
"""
function periodic_domain(field, ::Type{T}, n) where {T}
    per = periodic(field)
    lo = rangemin(field, zeros(T, 3))
    hi = rangemax(field, zeros(T, 3))

    xmin = fill(-T(Inf), n)
    xmax = fill(+T(Inf), n)

    for i in 1:3
        if per[i]
            xmin[i] = lo[i]
            xmax[i] = hi[i]
        end
    end

    (xmin, xmax)
end
