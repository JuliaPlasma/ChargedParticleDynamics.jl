using GeometricSolutions
import GeometricEquations
using ElectromagneticFields: FieldFunctions, coordinates, equilibrium

export cartesian_solution

"""
    check_chart(field, reference)

Throw an `ArgumentError` unless `field` holds an equilibrium of the same kind as `reference`, and
so one written in the same chart. A formula that an equilibrium module writes in its own
coordinates — its `toroidal_momentum`, which is `ϑ₃` where the third coordinate is the toroidal
angle and `x ϑ₂ - y ϑ₁` in a cartesian chart — is right only for a field in that chart, while the
right-hand sides take a field in any orthogonal chart.
"""
function check_chart(field, reference)
    a = nameof(typeof(equilibrium(field)))
    b = nameof(typeof(equilibrium(reference)))
    a === b ||
        throw(ArgumentError("this formula is written in the chart of a $b; the field holds a $a"))
    nothing
end

"""
    cartesian_solution(sol)
    cartesian_solution(sol, field)

The solution `sol` in cartesian coordinates, together with the major radius: a named tuple of the
time series and the `R`, `X`, `Y` and `Z` data series. The coordinate helpers come from
`coordinates(field)`, so the field must be axisymmetric; the one-argument form takes it from the
solution's own `parameters`.
"""
function cartesian_solution(sol)
    cartesian_solution(sol, GeometricEquations.parameters(sol.problem).field)
end

function cartesian_solution(sol, field::FieldFunctions)
    crd = coordinates(field)

    R = [crd.R.(sol.t[i], sol.q[i, 1], sol.q[i, 2], sol.q[i, 3]) for i in eachindex(sol.t)]
    X = [crd.X.(sol.t[i], sol.q[i, 1], sol.q[i, 2], sol.q[i, 3]) for i in eachindex(sol.t)]
    Y = [crd.Y.(sol.t[i], sol.q[i, 1], sol.q[i, 2], sol.q[i, 3]) for i in eachindex(sol.t)]
    Z = [crd.Z.(sol.t[i], sol.q[i, 1], sol.q[i, 2], sol.q[i, 3]) for i in eachindex(sol.t)]

    (
        t = sol.t,
        R = DataSeries(R),
        X = DataSeries(X),
        Y = DataSeries(Y),
        Z = DataSeries(Z)
    )
end
