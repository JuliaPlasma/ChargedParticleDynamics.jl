
export compute_energy, compute_toroidal_momentum, compute_momentum_error, compute_one_form, compute_error_drift


function convert_coordinates_RZphi_to_xyz(q::DataSeries{DT,2}) where {DT}
    x = SDataSeries(DT, q.nd, q.nt, q.ni)
    for k in 1:q.ni
        for j in 0:q.nt
            for i in 1:q.nd
                x[1,j,k] = q[1,j,k] * cos(q[3,j,k])
                x[2,j,k] = q[1,j,k] * sin(q[3,j,k])
                x[3,j,k] = q[2,j,k]
                x[4,j,k] = q[4,j,k]
            end
        end
    end
    return x
end

function compute_energy(sol::Solution, hamiltonian::Function)
    compute_energy(sol.t, sol.q, hamiltonian)
end

function compute_energy(t::TimeSeries{DT}, q::DataSeries{DT,2}, hamiltonian::Function) where {DT}
    h = SDataSeries(DT, 2, q.nt, 1)
    for i in 0:q.nt
        h[1,i] = hamiltonian(t[i], q[:,i])
        h[2,i] = (h[1,i] - h[1,0]) / h[1,0]
    end
    return h
end

function compute_toroidal_momentum(t::TimeSeries{DT}, q::DataSeries{DT,2}, toroidal_momentum::Function) where {DT}
    p = SDataSeries(DT, 2, q.nt, 1)
    for i in 0:q.nt
        p[1,i] = toroidal_momentum(t[i], q[:,i])
        p[2,i] = (p[1,i] - p[1,0]) / p[1,0]
    end
    return p
end

function compute_one_form(t::TimeSeries{DT}, q::DataSeries{DT,2}, ϑ::Function) where {DT}
    p = SDataSeries(DT, q.nd, q.nt, 1)
    theta = zeros(DT, q.nd)
    for j in 0:q.nt
        ϑ(t[j], q[:,j], theta)
        for i in 1:q.nd
            p[i,j] = theta[i]
        end
    end
    return p
end

function compute_momentum_error(p::DataSeries{DT}, ϑ::DataSeries{DT}) where {DT}
    @assert p.nd == ϑ.nd
    @assert p.nt == ϑ.nt
    @assert p.ni == ϑ.ni

    err = SDataSeries(DT, p.nd, p.nt, p.ni)

    for k in 1:p.ni
        for j in 0:p.nt
            for i in 1:p.nd
                err[i,j,k] = p[i,j,k] - ϑ[i,j,k]
            end
        end
    end

    return err
end

function compute_error_drift(t::TimeSeries{DT}, e::DataSeries{DT,2}, lint=100) where {DT}
    @assert t.n == e.nt

    nint  = div(t.n, lint)
    tint  = t[div(nint,2):nint:t.n]

    drift = SDataSeries(DT, 3, lint, 1)

    for i in 1:lint
        i1 = nint*(i-1)+1
        i2 = nint*i

        drift[1,i] = tint[i]
        drift[3,i] = maximum(e[2,i1:i2])
        drift[2,i] = drift[3,i] - drift[3,1]
    end

    return drift
end
