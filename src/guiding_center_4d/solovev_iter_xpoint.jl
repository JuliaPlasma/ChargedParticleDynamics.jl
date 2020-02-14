"""
Analytic ITER-like Solov'ev equilibrium with X-point.
"""
module GuidingCenter4dSolovevIterXpoint

    using ElectromagneticFields: load_equilibrium, SolovevXpointITER

    export initial_conditions_barely_passing, initial_conditions_barely_trapped,
           initial_conditions_deeply_passing, initial_conditions_deeply_trapped

    export hamiltonian, toroidal_momentum


    equ = SolovevXpointITER()
    load_equilibrium(equ; target_module=GuidingCenter4dSolovevIterXpoint)

    initial_conditions_barely_passing() = ([2.5, 0., 0., 3.425E-1], 1E-2)
    initial_conditions_barely_trapped() = ([2.5, 0., 0., 3.375E-1], 1E-2)
    initial_conditions_deeply_passing() = ([2.5, 0., 0., 5E-1], 1E-2)
    initial_conditions_deeply_trapped() = ([2.5, 0., 0., 1E-1], 1E-2)

    include("guiding_center_4d_common.jl")
    include("guiding_center_4d_equations.jl")

    function toroidal_momentum(t,q)
        R(t,q) * ϑ₃(t,q)
    end

    include("guiding_center_4d_diagnostics.jl")

end
