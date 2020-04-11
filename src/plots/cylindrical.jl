module Cylindrical

    using Plots
    using Plots.PlotMeasures
    using RecipesBase
    using LaTeXStrings

    using GeometricIntegrators.Solutions: TimeSeries, DataSeries, Solution
    using GeometricProblems.PlotRecipes: subscript

    include("common.jl")

    @userplot PlotTrajectory
    @recipe function f(p::PlotTrajectory; plottype=nothing, nplot=1, xlims=:auto, ylims=:auto, latex=true)
        if length(p.args) != 1 || !(typeof(p.args[1]) <: Solution)
            error("Guiding center trajectory plots should be given a solution. Got: $(typeof(p.args))")
        end

        sol = p.args[1]

        linewidth=2
        if sol.nt ≤ 1000
            markersize := 5
        else
            markersize  := 1
            markercolor := 1
            linecolor   := 1
            markerstrokewidth := 1
            markerstrokecolor := 1
        end

        seriestype := :scatter
        guidefont  := font(18)
        tickfont   := font(12)
        legend := :none


        if plottype == nothing
            size := (600,800)

            layout := @layout [Rplot
                               Zplot
                               ϕplot
                               Uplot]

            if latex
                ylabels = (L"R", L"Z", L"\varphi", L"u")
            else
                ylabels = ("R", "Z", "φ", "u")
            end

            for i in 1:4
                @series begin
                    subplot := i
                    if i == 4
                        if latex
                            xlabel := L"t"
                        else
                            xlabel := "t"
                        end
                    else
                        xaxis := false
                    end
                    ylabel := ylabels[i]
                    sol.t[0:nplot:end], sol.q[i,0:nplot:end]
                end
            end

        elseif plottype==:trajectoryRZ
            if latex
                xlabel := L"R"
                ylabel := L"Z"
            else
                xlabel := "R"
                ylabel := "Z"
            end
            xlims  := xlims
            ylims  := ylims
            size   := (500,500)

            sol.q[1,0:nplot:end], sol.q[2,0:nplot:end]

        elseif plottype==:trajectoryXY
            if latex
                xlabel := L"X"
                ylabel := L"Y"
            else
                xlabel := "X"
                ylabel := "Y"
            end
            xlims  := xlims
            ylims  := ylims
            size   := (500,500)

            sol.q[1,0:nplot:end] .* cos.(sol.q[3,0:nplot:end]), sol.q[1,0:nplot:end] .* sin.(sol.q[3,0:nplot:end])

        elseif plottype==:trajectoryYZ
            if latex
                xlabel := L"Y"
                ylabel := L"Z"
            else
                xlabel := "Y"
                ylabel := "Z"
            end
            xlims  := xlims
            ylims  := ylims
            size   := (500,500)

            sol.q[1,0:nplot:end] .* sin.(sol.q[3,0:nplot:end]), sol.q[2,0:nplot:end]

        elseif plottype==:trajectoryXZ
            if latex
                xlabel := L"X"
                ylabel := L"Z"
            else
                xlabel := "X"
                ylabel := "Z"
            end
            xlims  := xlims
            ylims  := ylims
            size   := (500,500)

            sol.q[1,0:nplot:end] .* cos.(sol.q[3,0:nplot:end]), sol.q[2,0:nplot:end]

        elseif plottype==:trajectory3d
            if latex
                xlabel := L"X"
                ylabel := L"Y"
                zlabel := L"Z"
            else
                xlabel := "X"
                ylabel := "Y"
                zlabel := "Z"
            end
            xlims  := xlims
            ylims  := ylims
            size --> (500,500)

            sol.q[1,0:nplot:end] .* cos.(sol.q[3,0:nplot:end]), sol.q[1,0:nplot:end] .* sin.(sol.q[3,0:nplot:end]), sol.q[2,0:nplot:end]
        end
    end

end
