
import DecFP
import ElectromagneticFields

using GeometricIntegrators.Equations


export guiding_center_4d_ode, guiding_center_4d_iode, guiding_center_4d_iode_λ,
       guiding_center_4d_iode_dec128, guiding_center_4d_vode,
       guiding_center_4d_dg, guiding_center_4d_formal_lagrangian



function guiding_center_4d_periodicity(q)
    p = zero(q)

    try
        p .= ElectromagneticFields.periodicity(q, equ)
    catch
        @warn "No equilibrium found to determine periodicity."
    end

    return p
end


function guiding_center_4d_ode(qᵢ, μ; periodic=true)
    guiding_center_4d_v_μ(t, q, v) = guiding_center_4d_v(t, q, μ, v)

    if periodic
        ODE(guiding_center_4d_v_μ, qᵢ; periodicity=guiding_center_4d_periodicity(qᵢ))
    else
        ODE(guiding_center_4d_v_μ, qᵢ)
    end
end


function guiding_center_4d_iode(qᵢ, μ; periodic=true)
    pᵢ = guiding_center_4d_pᵢ(qᵢ)

    guiding_center_4d_v_μ(t, q, p, v) = guiding_center_4d_v(t, q, p, μ, v)
    guiding_center_4d_f_μ(t, q, p, v) = guiding_center_4d_f(t, q, p, μ, v)

    if periodic
        IODE(guiding_center_4d_ϑ, guiding_center_4d_f_μ,
             guiding_center_4d_g, qᵢ, pᵢ;
             v=guiding_center_4d_v_μ, periodicity=guiding_center_4d_periodicity(qᵢ))
    else
        IODE(guiding_center_4d_ϑ, guiding_center_4d_f_μ,
             guiding_center_4d_g, qᵢ, pᵢ;
             v = guiding_center_4d_v_μ)
    end
end

function guiding_center_4d_iode_dec128(qᵢ, μ; periodic=true)
    guiding_center_4d_iode(Dec128.(qᵢ), Dec128.(μ); periodic=periodic)
end

function guiding_center_4d_iode_λ(qᵢ, μ; periodic=true)
    pᵢ = guiding_center_4d_pᵢ(qᵢ)
    λ₀ = guiding_center_4d_λ₀(qᵢ)

    guiding_center_4d_v_μ(t, q, p, v) = guiding_center_4d_v(t, q, p, μ, v)
    guiding_center_4d_f_μ(t, q, p, v) = guiding_center_4d_f(t, q, p, μ, v)

    if periodic
        IODE(guiding_center_4d_ϑ, guiding_center_4d_f_μ,
             guiding_center_4d_g, qᵢ, pᵢ, λ₀;
             v=guiding_center_4d_v_μ, periodicity=guiding_center_4d_periodicity(qᵢ))
    else
        IODE(guiding_center_4d_ϑ, guiding_center_4d_f_μ,
             guiding_center_4d_g, qᵢ, pᵢ, λ₀;
             v=guiding_center_4d_v_μ)
    end
end

function guiding_center_4d_vode(qᵢ, μ; periodic=true)
    pᵢ = guiding_center_4d_pᵢ(qᵢ)

    guiding_center_4d_v_μ(t, q, p, v) = guiding_center_4d_v(t, q, p, μ, v)
    guiding_center_4d_f_μ(t, q, p, v) = guiding_center_4d_f(t, q, p, μ, v)

    if periodic
        VODE(guiding_center_4d_ϑ, guiding_center_4d_f_μ,
             guiding_center_4d_g, guiding_center_4d_v_μ,
             ω, dH, qᵢ, pᵢ; periodicity=guiding_center_4d_periodicity(qᵢ))
    else
        VODE(guiding_center_4d_ϑ, guiding_center_4d_f_μ,
             guiding_center_4d_g, guiding_center_4d_v_μ,
             ω, dH, qᵢ, pᵢ)
    end
end

function guiding_center_4d_dg(qᵢ, μ; κ=0.0, periodic=true)
    if periodic
        periodicity = []
    else
        periodicity=guiding_center_4d_periodicity(qᵢ)
    end

    guiding_center_4d_ϑ_κ(t, q, v, p) = guiding_center_4d_ϑ(κ, t, q, v, p)
    guiding_center_4d_f_κ(t, q, v, f) = guiding_center_4d_f(κ, t, q, v, μ, f)
    guiding_center_4d_g_κ(t, q, λ, g) = guiding_center_4d_g(κ, t, q, λ, g)

    IODE(guiding_center_4d_ϑ_κ, guiding_center_4d_f_κ,
         guiding_center_4d_g_κ, guiding_center_4d_v,
         qᵢ, qᵢ; periodicity=periodicity)
end

function guiding_center_4d_formal_lagrangian(qᵢ, μ; periodic=true)
    pᵢ = guiding_center_4d_pᵢ(qᵢ)

    guiding_center_4d_v_μ(t, q, p, v) = guiding_center_4d_v(t, q, p, μ, v)
    guiding_center_4d_f_μ(t, q, p, v) = guiding_center_4d_f(t, q, p, μ, v)

    if periodic
        VODE(guiding_center_4d_ϑ, guiding_center_4d_f_μ, guiding_center_4d_g, guiding_center_4d_v_μ,
             ω, dH, qᵢ, pᵢ; periodicity=guiding_center_4d_periodicity(qᵢ))
    else
        VODE(guiding_center_4d_ϑ, guiding_center_4d_f_μ, guiding_center_4d_g, guiding_center_4d_v_μ,
             ω, dH, qᵢ, pᵢ)
    end
end
