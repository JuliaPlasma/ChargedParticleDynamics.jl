
import ElectromagneticFields

import GeometricIntegrators.Equations: ODE, SODE


export guiding_center_4d_ode, guiding_center_4d_sode


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
    ODE(v, qᵢ; parameters=params, h=hamiltonian,
            periodicity=guiding_center_4d_periodicity(qᵢ, periodic))
end


function guiding_center_4d_sode(qᵢ, params; periodic=true)
    SODE((v₁, v₂, v₃, v₄, v₅, v₆), qᵢ; parameters=params,
            periodicity=guiding_center_4d_periodicity(qᵢ, periodic))
end
