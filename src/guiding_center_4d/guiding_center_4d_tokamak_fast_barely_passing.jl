"""
Fast barely passing particle in analytic, axisymmetric tokamak equilibrium.
"""
module TokamakFastBarelyPassing

    import ElectromagneticFields

    export guiding_center_4d_ode, guiding_center_4d_iode, guiding_center_4d_iode_λ, guiding_center_4d_iode_dec128,
           hamiltonian, toroidal_momentum, u, ω, α, α1, α2, α3, α4, dα, β, β1, β2, β3, B, B₁, B₂, B₃, b₁, b₂, b₃, dH

    # Δt=2.5, nt=50

    const R₀ = 2
    const B₀ = 5
    const q  = 2

    const μ  = 1E-2
    const q₀ = [2.5, 0., 0., 3.425E-1]

    ElectromagneticFields.load_equilibrium(ElectromagneticFields.AxisymmetricTokamakCylindrical(R₀, B₀, q); target_module=TokamakFastBarelyPassing)

    include("guiding_center_4d_coords_RZphi.jl")

end
