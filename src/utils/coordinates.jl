using GeometricSolutions

export cartesian_solution


function cartesian_solution(sol, equ)
    R = [equ.R.(sol.t[i], sol.q[i,1], sol.q[i,2], sol.q[i,3]) for i in eachindex(sol.t)]
    X = [equ.X.(sol.t[i], sol.q[i,1], sol.q[i,2], sol.q[i,3]) for i in eachindex(sol.t)]
    Y = [equ.Y.(sol.t[i], sol.q[i,1], sol.q[i,2], sol.q[i,3]) for i in eachindex(sol.t)]
    Z = [equ.Z.(sol.t[i], sol.q[i,1], sol.q[i,2], sol.q[i,3]) for i in eachindex(sol.t)]

    (
        t = sol.t,
        R = DataSeries(R),
        X = DataSeries(X),
        Y = DataSeries(Y),
        Z = DataSeries(Z),
    )
end
