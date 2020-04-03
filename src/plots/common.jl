
@userplot PlotToroidalMomentumError
@recipe function f(p::PlotToroidalMomentumError; nplot=1, latex=true)
    if length(p.args) == 2 && typeof(p.args[1]) <: TimeSeries && typeof(p.args[2]) <: DataSeries
        t  = p.args[1]
        pφ = p.args[2]
        @assert t.n == pφ.nt
    else
        error("Toroidal momentum error plot should be given a timeseries and a data series. Got: $(typeof(p.args))")
    end

    legend := :none
    size := (800,400)

    @series begin
        if latex
            xlabel := L"t"
            ylabel := L"[P_\varphi (t) - P_\varphi (0)] / P_\varphi (0)"
        else
            xlabel := "t"
            ylabel := "[P(t) - P(0)] / P(0)"
        end
        xlims  := (t[0], Inf)
        yformatter := :scientific
        guidefont := font(18)
        tickfont := font(12)
        right_margin := 10mm
        t[0:nplot:end], pφ[0:nplot:end]
    end
end
