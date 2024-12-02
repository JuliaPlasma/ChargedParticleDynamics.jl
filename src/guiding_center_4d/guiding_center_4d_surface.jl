
using PoincareInvariants

export guiding_center_4d_surface_ode,
       guiding_center_4d_surface_iode,
       guiding_center_4d_surface_iode_λ

export guiding_center_4d_ode_poincare_invariant_2nd,
       guiding_center_4d_ode_poincare_invariant_2nd_trapezoidal,
       guiding_center_4d_iode_poincare_invariant_2nd,
       guiding_center_4d_iode_poincare_invariant_2nd_trapezoidal,
       guiding_center_4d_iode_poincare_invariant_2nd_λ


function f_surface(i, j, nx, ny)
    f_surface(i/nx, j/ny)
end

function initial_conditions_surface(nx, ny)
    q₀ = [zeros(4) for _ in 1:nx*ny]

    for j in 1:ny
        for i in 1:nx
            q₀[nx*(j-1)+i] = f_surface(i, j, nx, ny)
        end
    end

    return q₀
end


guiding_center_4d_surface_ode_init(q₀) = guiding_center_4d_ode(q₀, (μ = μ_surface(),); periodic=false)
guiding_center_4d_surface_iode_init(q₀) = guiding_center_4d_iode(q₀, (μ = μ_surface(),); periodic=false)
guiding_center_4d_surface_iode_λ_init(q₀) = guiding_center_4d_iode_λ(q₀, (μ = μ_surface(),); periodic=false)


function guiding_center_4d_surface_ode(nx, ny)
    guiding_center_4d_surface_ode_init(initial_conditions_surface(nx, ny))
end

function guiding_center_4d_surface_iode(nx, ny)
    guiding_center_4d_surface_iode_init(initial_conditions_surface(nx, ny))
end

function guiding_center_4d_surface_iode_λ(nx, ny)
    guiding_center_4d_surface_iode_λ_init(initial_conditions_surface(nx, ny))
end


function guiding_center_4d_ode_poincare_invariant_2nd(Δt, nx, ny, ntime, nsave, DT=Float64)
    PoincareInvariant2nd(guiding_center_4d_surface_ode_init, f_surface, ω, Δt, 4, nx, ny, ntime, nsave, DT)
end

function guiding_center_4d_ode_poincare_invariant_2nd_trapezoidal(Δt, nx, ny, ntime, nsave, DT=Float64)
    PoincareInvariant2ndTrapezoidal(guiding_center_4d_surface_ode_init, f_surface, ω, Δt, 4, nx, ny, ntime, nsave, DT)
end

function guiding_center_4d_iode_poincare_invariant_2nd(Δt, nx, ny, ntime, nsave, DT=Float64)
    PoincareInvariant2nd(guiding_center_4d_surface_iode_init, f_surface, ω, Δt, 4, nx, ny, ntime, nsave, DT)
end

function guiding_center_4d_iode_poincare_invariant_2nd_trapezoidal(Δt, nx, ny, ntime, nsave, DT=Float64)
    PoincareInvariant2ndTrapezoidal(guiding_center_4d_surface_iode_init, f_surface, ω, (D²ϑ₁, D²ϑ₂, D²ϑ₃, D²ϑ₄), Δt, 4, nx, ny, ntime, nsave)
end

function guiding_center_4d_iode_poincare_invariant_2nd_λ(Δt, nx, ny, ntime, nsave, DT=Float64)
    PoincareInvariant2nd(guiding_center_4d_surface_iode_init, f_surface, ω, (D²ϑ₁, D²ϑ₂, D²ϑ₃, D²ϑ₄), Δt, 4, nx, ny, ntime, nsave, DT)
end
