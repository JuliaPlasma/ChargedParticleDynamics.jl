"""
Fast deeply trapped particle in analytic, axisymmetric tokamak equilibrium.
"""
module TokamakFastDeeplyTrapped

    import ElectromagneticFields

    export guiding_center_4d_ode, guiding_center_4d_iode, guiding_center_4d_iode_λ,
           hamiltonian, toroidal_momentum, u, ω, α, α1, α2, α3, α4, dα, β, β1, β2, β3, b1, b2, b3, dH

    # Δt=5, nt=50

    const R₀ = 2
    const B₀ = 5
    const q  = 2

    const μ  = 1E-2
    const q₀ = [2.5, 0., 0., 1E-1]

    ElectromagneticFields.load_equilibrium(ElectromagneticFields.AxisymmetricTokamakCylindrical(R₀, B₀, q); target_module=TokamakFastDeeplyTrapped)

    include("guiding_center_4d_coords_RZphi.jl")

end
