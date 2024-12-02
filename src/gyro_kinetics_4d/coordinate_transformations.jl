

export transform_q_to_q̃!, transform_q_to_q̃,
       transform_q̃_to_q!, transform_q̃_to_q


function transform_q_to_q̃!(q̃, t, q, params)
    q̃ .= q .* params.ω₀
end

function transform_q̃_to_q!(q, t, q̃, params)
    q .= q̃ ./ params.ω₀
end

transform_q_to_q̃(t, q, params) = transform_q_to_q̃!(zero(q), t, q, params)
transform_q̃_to_q(t, q̃, params) = transform_q̃_to_q!(zero(q̃), t, q̃, params)
