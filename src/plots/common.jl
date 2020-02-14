
@userplot PlotToroidalMomentumError
@recipe function f(p::PlotToroidalMomentumError)
    if length(p.args) == 1 && typeof(p.args[1]) <: Solution
        sol = p.args[1]
        t = sol.t
        pφ = compute_toroidal_momentum(sol)
    elseif length(p.args) == 2 && typeof(p.args[1]) <: TimeSeries && typeof(p.args[2]) <: DataSeries
        t  = p.args[1]
        pφ = p.args[2]
        @assert t.n == pφ.nt
    else
        error("Toroidal momentum error plot should be given a solution or a timeseries and a data series. Got: $(typeof(p.args))")
    end

    legend := :none
    size := (800,400)

    @series begin
        xlabel := L"t"
        ylabel := L"[P_\varphi (t) - P_\varphi (0)] / P_\varphi (0)"
        xlims  := (t[0], Inf)
        yformatter := :scientific
        guidefont := font(18)
        tickfont := font(12)
        right_margin := 10mm
        t, pφ[2,:]
    end
end
