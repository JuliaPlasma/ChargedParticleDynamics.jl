
import DecFP
import ElectromagneticFields

using GeometricIntegrators.Equations


export guiding_center_4d_ode, guiding_center_4d_iode, guiding_center_4d_iode_λ,
       guiding_center_4d_iode_dec128, guiding_center_4d_vode,
       guiding_center_4d_dg, guiding_center_4d_formal_lagrangian



function guiding_center_4d_periodicity(q, periodic=true)
    periodicity = zeros(eltype(q), size(q,1))

    if periodic
        try
            periodicity .= ElectromagneticFields.periodicity(q, equ)
        catch
            @warn "No equilibrium found to determine periodicity."
        end
    end

    return periodicity
end


function guiding_center_4d_ode(qᵢ, params; periodic=true)
    ODE(guiding_center_4d_v, qᵢ;
            parameters=params, h=hamiltonian,
            periodicity=guiding_center_4d_periodicity(qᵢ, periodic))
end


function guiding_center_4d_iode(qᵢ, params; periodic=true)
    pᵢ = guiding_center_4d_pᵢ(qᵢ)

    IODE(guiding_center_4d_ϑ, guiding_center_4d_f, guiding_center_4d_g, qᵢ, pᵢ;
            parameters=params, h=hamiltonian, v=guiding_center_4d_v,
            periodicity=guiding_center_4d_periodicity(qᵢ, periodic))
end

function guiding_center_4d_iode_dec128(qᵢ, params; periodic=true)
    guiding_center_4d_iode(Dec128.(qᵢ), Dec128.(μ); periodic=periodic)
end

function guiding_center_4d_iode_λ(qᵢ, params; periodic=true)
    pᵢ = guiding_center_4d_pᵢ(qᵢ)
    λ₀ = guiding_center_4d_λ₀(qᵢ)

    IODE(guiding_center_4d_ϑ, guiding_center_4d_f, guiding_center_4d_g, qᵢ, pᵢ, λ₀;
            parameters=params, h=hamiltonian, v=guiding_center_4d_v,
            periodicity=guiding_center_4d_periodicity(qᵢ, periodic))
end

function guiding_center_4d_vode(qᵢ, params; periodic=true)
    pᵢ = guiding_center_4d_pᵢ(qᵢ)

    VODE(guiding_center_4d_ϑ, guiding_center_4d_f, guiding_center_4d_g, qᵢ, pᵢ;
            parameters=params, h=hamiltonian, v=guiding_center_4d_v,
            Ω=guiding_center_4d_ω, ∇H=guiding_center_4d_dH,
            periodicity=guiding_center_4d_periodicity(qᵢ, periodic))
end

function guiding_center_4d_dg(qᵢ, params; κ=0.0, periodic=true)
    guiding_center_4d_ϑ_κ(t, q, v, p, params) = guiding_center_4d_ϑ(t, q, v, p, params, κ)
    guiding_center_4d_f_κ(t, q, v, f, params) = guiding_center_4d_f(t, q, v, f, params, κ)
    guiding_center_4d_g_κ(t, q, λ, g, params) = guiding_center_4d_g(t, q, λ, g, params, κ)

    IODE(guiding_center_4d_ϑ_κ, guiding_center_4d_f_κ, guiding_center_4d_g_κ, qᵢ, qᵢ;
            parameters=params, h=hamiltonian, v=guiding_center_4d_v,
            periodicity=guiding_center_4d_periodicity(qᵢ, periodic))
end

function guiding_center_4d_formal_lagrangian(qᵢ, params; periodic=true)
    pᵢ = guiding_center_4d_pᵢ(qᵢ)

    VODE(guiding_center_4d_ϑ, guiding_center_4d_f, guiding_center_4d_g, qᵢ, pᵢ;
            parameters=params, h=hamiltonian, v=guiding_center_4d_v,
            Ω=guiding_center_4d_ω, ∇H=guiding_center_4d_dH,
            periodicity=guiding_center_4d_periodicity(qᵢ, periodic))
end
