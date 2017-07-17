module TokamakBarelyPassing

    import MagneticEquilibria

    # export load_equilibrium
    export guiding_center_4d_ode, guiding_center_4d_iode,
           hamiltonian, toroidal_momentum, u, α, α1, α2, α3, α4, dα, β, β1, β2, β3, b1, b2, b3, dH

    # function load_equilibrium(equ)
    #     MagneticEquilibria.load_equilibrium(equ; target_module=TokamakBarelyPassing)
    # end

    const R0 = 1.
    const B0 = 1.
    const q  = 2.

    const μ  = 2.448E-6
    const q₀ = [1.05, 0., 0., 8.117E-4]

    MagneticEquilibria.load_equilibrium(MagneticEquilibria.AxisymmetricTokamak(R0, B0, q); target_module=TokamakBarelyPassing)

    include("guiding_center_4d_RZphi.jl")

end
