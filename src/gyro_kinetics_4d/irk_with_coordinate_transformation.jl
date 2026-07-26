#
# REFERENCE MATERIAL — NOT PART OF THE PACKAGE.
#
# This file is not `include`d by any module and does not compile. It is the 0.x-era sketch of the
# coordinate-transforming Runge-Kutta integrator described in `TODO.md`, kept as the starting point
# for the port onto the current `IRK`. It is written against a `GeometricIntegrators` API removed
# several major versions ago — `Parameters`, `ODEIntegratorCache`, `CacheDict`,
# `create_nonlinear_solver`, `function_stages!`, `AtomicSolutionODE`,
# `GeometricIntegrators.Integrators` — and the transformation itself was never finished: all three
# `coordinate_transformation_*` functions below are commented out and the one that would apply the
# transform is an empty stub.
#
# What is worth carrying over:
#
#   * the `(Q, Q̃)` stage-pair cache layout of `IntegratorCacheFIRKwCT`, with the vector field
#     evaluated on the untransformed stage and the solve carried out on the transformed one;
#   * the `compute_stages!` / `function_stages!` split, which is spelled `components!` / `residual!`
#     in the current integrator;
#   * the commented-out `coordinate_transformation_rhs!`, which spells out the implicit relation
#     `b = q̃ - B*∥(q) q` that has to be solved for `q` inside each step.
#
# The internal names still use the old `FIRK` spelling; the current scheme is `IRK`. Renaming them
# belongs with the port, not with keeping the sketch on file. See `TODO.md`, *A coordinate-
# transforming IRK integrator*.
#

using ForwardDiff

# using GeometricIntegrators
# using GeometricIntegrators.Integrators
# using SimpleSolvers

# using GeometricIntegrators.Integrators: Parameters, CacheDict, ODEIntegratorCache
# using GeometricIntegrators.Integrators: equation, equations, tableau, timestep, eachdim, eachstage
# using GeometricIntegrators.Integrators: create_nonlinear_solver, create_nonlinear_solver_with_jacobian, create_internal_stage_vector
# using GeometricIntegrators.Integrators: update_solution!, update_vector_fields!

using .GuidingCenter4dSolovevIterXpoint


"Parameters for right-hand side function of fully implicit Runge-Kutta methods."
mutable struct ParametersFIRKwCT{DT, TT, D, S, ET <: NamedTuple, PT <: NamedTuple, FT, JT} <: Parameters{DT,TT}
    equs::ET
    params::PT
    tab::Tableau{TT}
    Δt::TT

    F::FT
    Jconfig::JT

    t::TT
    q::Vector{DT}

    # t̃::TT
    # x₀::Vector{DT}

    function ParametersFIRKwCT{DT,D}(equs::ET, params::PT, tab::Tableau{TT}, Δt::TT) where {DT, TT, D, ET <: NamedTuple, PT <: NamedTuple}
        F = (v,q) -> equs[:v](zero(TT), q, v)
        tq = zeros(DT,D)
        tv = zeros(DT,D)
        Jconfig = ForwardDiff.JacobianConfig(F, tv, tq)

        new{DT, TT, D, tab.s, ET, PT, typeof(F), typeof(Jconfig)}(equs, params, tab, Δt, F, Jconfig, zero(TT), zeros(DT,D))#, zero(TT), x₀)
    end
end


@doc raw"""
Fully implicit Runge-Kutta integrator cache.

### Fields

* `q̃`: initial guess of solution
* `ṽ`: initial guess of vector field
* `s̃`: holds shift due to periodicity of solution
* `Q`: internal stages of solution
* `V`: internal stages of vector field
* `Y`: vector field of internal stages
"""
struct IntegratorCacheFIRKwCT{DT,D,S} <: ODEIntegratorCache{DT,D}
    q::Vector{DT}
    q̅::Vector{DT}

    q̃::Vector{DT}
    ṽ::Vector{DT}
    s̃::Vector{DT}

    Q::Vector{Vector{DT}}
    Q̃::Vector{Vector{DT}}
    Ṽ::Vector{Vector{DT}}
    Ỹ::Vector{Vector{DT}}
    
    J::Vector{Matrix{DT}}

    function IntegratorCacheFIRKwCT{DT,D,S}() where {DT,D,S}
        q = zeros(DT,D)
        q̅ = zeros(DT,D)

        q̃ = zeros(DT,D)
        ṽ = zeros(DT,D)
        s̃ = zeros(DT,D)

        Q = create_internal_stage_vector(DT, D, S)
        Q̃ = create_internal_stage_vector(DT, D, S)
        Ṽ = create_internal_stage_vector(DT, D, S)
        Ỹ = create_internal_stage_vector(DT, D, S)
        
        J = [zeros(DT,D,D) for i in 1:S]

        new(q, q̅, q̃, ṽ, s̃, Q, Q̃, Ṽ, Ỹ, J)
    end
end

function IntegratorCache{ST}(params::ParametersFIRKwCT{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCacheFIRKwCT{ST,D,S}(; kwargs...)
end

@inline GeometricIntegrators.Integrators.CacheType(ST, params::ParametersFIRKwCT{DT,TT,D,S}) where {DT,TT,D,S} = IntegratorCacheFIRKwCT{ST,D,S}


"Fully implicit Runge-Kutta integrator with coordinate transformation."
struct IntegratorFIRKwCT{DT, TT, D, S,
                PT <: ParametersFIRKwCT{DT,TT},
                ST <: NonlinearSolver{DT},
                IT <: InitialGuessODE{TT}}# <: IntegratorRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorFIRKwCT(params::ParametersFIRKwCT{DT,TT,D,S}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,ST,IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorFIRKwCT{DT,D}(equations::NamedTuple, parameters::NamedTuple, tableau::Tableau{TT}, Δt::TT; exact_jacobian=false) where {DT,TT,D}
        # get number of stages
        S = tableau.s

        # create solver for coordinate transformation
        # ct_solver = get_config(:nls_solver)(zeros(DT,3), (x,b) -> coordinate_transformation_rhs!(x, b, parameters))

        # create params
        params = ParametersFIRKwCT{DT,D}(equations, parameters, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create solver
        if exact_jacobian
            solver = create_nonlinear_solver_with_jacobian(DT, D*S, params, caches)
        else
            solver = create_nonlinear_solver(DT, D*S, params, caches)
        end

        # create initial guess
        iguess = InitialGuessODE{DT,D}(get_config(:ig_interpolation), equations[:v], Δt)

        # create integrator
        IntegratorFIRKwCT(params, solver, iguess, caches)
    end

    function IntegratorFIRKwCT{DT,D}(v::Function, ωabs::Function, parameters::NamedTuple, tableau::Tableau{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorFIRKwCT{DT,D}(NamedTuple{(:v,:ωabs)}((v, (t,q) -> ωabs(t,q,parameters))), parameters, tableau, Δt; kwargs...)
    end

    function IntegratorFIRKwCT(equation::ODE{DT,TT}, ωabs::Function, tableau::Tableau{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorFIRKwCT{DT, ndims(equation)}(equation.v, ωabs, equation.parameters, tableau, Δt; kwargs...)
    end
end


@inline Base.ndims(int::IntegratorFIRKwCT{DT,TT,D,S}) where {DT,TT,D,S} = D


function update_params!(int::IntegratorFIRKwCT, sol::AtomicSolutionODE)
    # set time for nonlinear solver and copy previous solution
    int.params.t  = sol.t
    int.params.q .= sol.q
end


# function coordinate_transformation_ig!(x, params::ParametersFIRKwCT)
#     # compute B⋆∥(t,q)
#     B = params.equs[:ωabs](params.t̃, [params.x₀, params.u], params.params)

#     # compute initial guess for transformed coordinate
#     for k in eachindex(x)
#         x[k] = params.q̃[k] / B
#     end
# end

# function coordinate_transformation_rhs!(x::Vector{ST}, b::Vector{ST}, params::ParametersFIRKwCT{DT,TT,D},
#                           caches::CacheDict) where {ST,DT,TT,D}

#     # get cache for internal stages
#     cache = caches[ST]

#     # copy nonlinear solver solution x to Q
#     for k in eachindex(cache.q)
#         cache.q[k] = x[k]
#     end

#     # compute B⋆∥(t,q)
#     B = params.equs[:ωabs](params.t̃, cache.q, params.params)

#     # compute b = q̃-q*B(q)
#     for k in eachindex(b)
#         b[k] = cache.q̃[k] - B * cache.q[k]
#     end
# end

# function coordinate_transformation!(q̃::Vector{ST}, q::Vector{ST}, params::ParametersFIRKwCT{DT,TT,D}) where {ST,DT,TT,D}

# end


function GeometricIntegrators.Integrators.initialize!(int::IntegratorFIRKwCT{DT}, sol::AtomicSolutionODE{DT},
                         cache::IntegratorCacheFIRKwCT{DT}=int.caches[DT]) where {DT}

    sol.t̅ = sol.t - timestep(int)

    # transform_q̃_to_q!(sol.q, cache.q, int.params.params)

    equations(int)[:v](sol.t, sol.q, sol.v)

    initialize!(int.iguess, sol.t, sol.q, sol.v,
                            sol.t̅, sol.q̅, sol.v̅)
end


function initial_guess!(int::IntegratorFIRKwCT{DT}, sol::AtomicSolutionODE{DT},
                        cache::IntegratorCacheFIRKwCT{DT}=int.caches[DT]) where {DT}

    local offset::Int

    # transform_q̃_to_q!(sol.q, cache.q, int.params.params)
    # transform_q̃_to_q!(sol.q̅, cache.q̅, int.params.params)

    # compute initial guess for internal stages
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q, sol.v, sol.q̅, sol.v̅, cache.Q̃[i], cache.Ṽ[i], tableau(int).q.c[i])
        # transform_q_to_q̃!(cache.Q[i], cache.Q̃[i], int.params.params)
    end
    for i in eachstage(int)
        offset = ndims(int)*(i-1)
        for k in eachdim(int)
            int.solver.x[offset+k] = 0
            for j in eachstage(int)
                int.solver.x[offset+k] += timestep(int) * tableau(int).q.a[i,j] * cache.Ṽ[j][k]
            end
        end
    end
end


function compute_stages!(x::Vector{ST}, Q::Vector{Vector{ST}}, Q̃::Vector{Vector{ST}}, Ṽ::Vector{Vector{ST}}, Ỹ::Vector{Vector{ST}},
                         params::ParametersFIRKwCT{DT,TT,D}) where {ST,DT,TT,D}

    local tᵢ::TT

    # copy x to Y and compute Q = q + Δt Y
    for i in eachindex(Q̃,Ỹ)
        for k in eachindex(Q̃[i],Ỹ[i])
            Ỹ[i][k] = x[D*(i-1)+k]
            Q̃[i][k] = params.q[k] + Ỹ[i][k]
        end
        # transform_q̃_to_q!(Q̃[i], Q[i], params.params)
    end

    # compute V = v(Q)
    for i in eachindex(Q̃,Ṽ)
        tᵢ = params.t + params.Δt * params.tab.q.c[i]
        params.equs[:v](tᵢ, Q̃[i], Ṽ[i])
    end
end

"Compute stages of fully implicit Runge-Kutta methods."
function GeometricIntegrators.Integrators.function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersFIRKwCT{DT,TT,D},
                          caches::CacheDict) where {ST,DT,TT,D}
    # temporary variables
    local y::ST

    # get cache for internal stages
    cache = caches[ST]

    # compute stages from nonlinear solver solution x
    compute_stages!(x, cache.Q, cache.Q̃, cache.Ṽ, cache.Ỹ, params)

    # compute b = - (Y-AV)
    for i in eachindex(cache.Ỹ)
        for k in eachindex(cache.Ỹ[i])
            y = 0
            for j in eachindex(cache.Ṽ)
                y += params.tab.q.a[i,j] * cache.Ṽ[j][k]
            end
            b[D*(i-1)+k] = - cache.Ỹ[i][k] + params.Δt * y
        end
    end
end

# "Compute stages of fully implicit Runge-Kutta methods."
# function jacobian!(x::Vector{DT}, jac::Matrix{DT}, cache::IntegratorCacheFIRKwCT{DT,D,S},
#                    params::ParametersFIRKwCT{DT,TT,D,S}) where {DT,TT,D,S}
#     local tᵢ::TT

#     for i in 1:S
#         for k in 1:D
#             cache.Ỹ[i][k] = x[D*(i-1)+k]
#             cache.Q̃[i][k] = params.q[k] + cache.Ỹ[i][k]
#         end
#         transform_q̃_to_q!(cache.Q̃[i], cache.Q[i], params.params)
#         tᵢ = params.t + params.Δt * params.tab.q.c[i]
#         F = (v,q) -> params.equs[:v](tᵢ, q, v)
#         ForwardDiff.jacobian!(cache.J[i], params.F, cache.ṽ, cache.Q[i], params.Jconfig)
#     end

#     jac .= 0

#     @inbounds for j in 1:S
#         for i in 1:S
#             for l in 1:D
#                 for k in 1:D
#                     m = D*(i-1)+k
#                     n = D*(j-1)+l

#                     if i == j && k == l
#                         jac[m,n] += -1
#                     end

#                     jac[m,n] += params.Δt * params.tab.q.a[i,j] * cache.J[j][k,l]
#                     jac[m,n] += params.Δt * params.tab.q.â[i,j] * cache.J[j][k,l]
#                 end
#             end
#         end
#     end
# end


"Integrate ODE with fully implicit Runge-Kutta integrator."
function GeometricIntegrators.Integrators.integrate_step!(int::IntegratorFIRKwCT{DT,TT}, sol::AtomicSolutionODE{DT,TT},
                         cache::IntegratorCacheFIRKwCT{DT}=int.caches[DT]) where {DT,TT}

    # update nonlinear solver parameters from atomic solution
    update_params!(int, sol)

    # compute initial guess
    initial_guess!(int, sol, cache)

    # reset atomic solution
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector field at internal stages
    compute_stages!(int.solver.x, cache.Q, cache.Q̃, cache.Ṽ, cache.Ỹ, int.params)

    # compute final update
    update_solution!(sol.q, sol.q̃, cache.Ṽ, tableau(int).q.b, tableau(int).q.b̂, timestep(int))

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.q)
end
