
export transform_q_to_q̃!, transform_q_to_q̃,
       transform_q̃_to_q!, transform_q̃_to_q


@doc raw"""
    transform_q_to_q̃(t, q, params)
    transform_q_to_q̃!(q̃, t, q, params)

Uniform rescaling of the state by the reference value ``\omega_{0}`` of ``B^{\star}_{\parallel}``
carried in `params`, ``\tilde{q} = \omega_{0} \, q``, and its inverse
[`transform_q̃_to_q!`](@ref).

!!! note "Not applied by the equations of motion"
    These are utilities, kept because a coordinate-transforming integrator was once written
    against them. The problem constructors work in the untransformed state ``q = (R, u)``.

    The notes this module follows write the rescaled characteristics as
    ``\tilde{R} = B^{\star}_{\parallel} R``, ``\tilde{V}_{\parallel} = B^{\star}_{\parallel}
    V_{\parallel}``, which reads like a change of variables but is not one: what the derivation
    actually does is absorb the phasespace Jacobian into the distribution function,
    ``\tilde{f} = B^{\star}_{\parallel} f``, leaving the *same trajectories* traversed in the
    reparametrised time ``dt = B^{\star}_{\parallel} \, ds``. That factor is already carried by
    the vector field [`v`](@ref); rescaling the state on top of it describes a different system.

    Note also that ``B^{\star}_{\parallel}`` is a function of ``(R, u)``, so freezing it at the
    initial condition — which is what `params.ω₀` is — would in any case not reproduce the
    substitution above.
"""
function transform_q_to_q̃!(q̃, t, q, params)
    q̃ .= q .* params.ω₀
end

"""
    transform_q̃_to_q(t, q̃, params)
    transform_q̃_to_q!(q, t, q̃, params)

Inverse of [`transform_q_to_q̃!`](@ref).
"""
function transform_q̃_to_q!(q, t, q̃, params)
    q .= q̃ ./ params.ω₀
end

"Allocating form of [`transform_q_to_q̃!`](@ref)."
transform_q_to_q̃(t, q, params) = transform_q_to_q̃!(zero(q), t, q, params)

"Allocating form of [`transform_q̃_to_q!`](@ref)."
transform_q̃_to_q(t, q̃, params) = transform_q̃_to_q!(zero(q̃), t, q̃, params)
