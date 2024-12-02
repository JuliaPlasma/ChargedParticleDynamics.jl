module Cartesian

    using LaTeXStrings
    using Plots
    using Plots.PlotMeasures
    using RecipesBase

    using GeometricSolutions: GeometricSolution

    include("common.jl")

    @userplot PlotTrajectory
    @recipe function f(p::PlotTrajectory; plottype=nothing)
        if length(p.args) != 1 || !(typeof(p.args[1]) <: Solution)
            error("Guiding center trajectory plots should be given a solution. Got: $(typeof(p.args))")
        end

        sol = p.args[1]

        if sol.nt ≤ 200
            markersize=5
        else
            markersize=1
        end

        seriestype := :scatter
        legend := :none

        if plottype == nothing
            size := (600,800)

            layout := @layout [Xplot
                               Yplot
                               Zplot
                               Uplot]

            ylabels = (L"X", L"Y", L"Z", L"u")

            for i in 1:4
                @series begin
                    xlabel := L"t"
                    ylabel := ylabels[i]
                    subplot := i
                    sol.t.t[:], sol.q.d[i,:]
                end
            end

        elseif plottype==:trajectoryXY
            xlabel --> L"X"
            ylabel --> L"Y"
            size --> (500,500)

            sol.q.d[1,:], sol.q.d[2,:]

        elseif plottype==:trajectoryYZ
            xlabel --> L"Y"
            ylabel --> L"Z"
            size --> (500,500)

            sol.q.d[2,:], sol.q.d[3,:]

        elseif plottype==:trajectoryXZ
            xlabel --> L"X"
            ylabel --> L"Z"
            size --> (500,500)

            sol.q.d[1,:], sol.q.d[3,:]

        elseif plottype==:trajectory3d
            xlabel --> L"X"
            ylabel --> L"Y"
            zlabel --> L"Z"

            size --> (500,500)

            sol.q.d[1,:], sol.q.d[2,:], sol.q.d[3,:]
        end
    end

end
