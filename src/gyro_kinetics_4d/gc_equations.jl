
import ElectromagneticFields

import GeometricIntegrators.Equations: ODE, SODE


export guiding_center_4d_ode, guiding_center_4d_sode


function guiding_center_4d_periodicity(q, periodic=true)
    period = zeros(eltype(q), size(q,1))

    if periodic
        period[1:3] .= periodicity(q)
    end

    return period
end


transform_v(V) = (t,q,q̇,parameters) -> V(t, transform_q̃_to_q(t, q, parameters), q̇, parameters)

function guiding_center_4d_ode(qᵢ, params; periodic=true)
    ODE(transform_v(v), qᵢ; parameters=params, h=hamiltonian,
            periodicity=guiding_center_4d_periodicity(qᵢ, periodic))
end

function guiding_center_4d_sode(qᵢ, params; periodic=true)
    SODE(Tuple(transform_v(V) for V in (v₁, v₂, v₃, v₄, v₅, v₆)), qᵢ; parameters=params,
            periodicity=guiding_center_4d_periodicity(qᵢ, periodic))
end
