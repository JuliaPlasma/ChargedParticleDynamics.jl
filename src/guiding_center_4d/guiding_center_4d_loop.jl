
using PoincareInvariants

export guiding_center_4d_loop_ode,
       guiding_center_4d_loop_iode,
       guiding_center_4d_loop_vode

export guiding_center_4d_ode_poincare_invariant_1st,
       guiding_center_4d_iode_poincare_invariant_1st,
       guiding_center_4d_vode_poincare_invariant_1st


function f_loop(i, n)
   f_loop(i/n)
end

function initial_conditions_loop(n)
   [f_loop(i,n) for i in 1:n]
end


guiding_center_4d_loop_ode_init(q₀) = guiding_center_4d_ode(q₀; parameters = (μ = μ_loop(),), periodic=false)
guiding_center_4d_loop_iode_init(q₀) = guiding_center_4d_iode(q₀; parameters = (μ = μ_loop(),), periodic=false)
guiding_center_4d_loop_vode_init(q₀) = guiding_center_4d_vode_formal_lagrangian(q₀; parameters = (μ = μ_loop(),), periodic=false)


function guiding_center_4d_loop_ode(n)
   guiding_center_4d_loop_ode_init(initial_conditions_loop(n))
end

function guiding_center_4d_loop_iode(n)
   guiding_center_4d_loop_iode_init(initial_conditions_loop(n))
end

function guiding_center_4d_loop_vode(n)
   guiding_center_4d_loop_iode_init(initial_conditions_loop(n))
end


function guiding_center_4d_ode_poincare_invariant_1st(Δt, nloop, ntime, nsave, DT=Float64)
   PoincareInvariant1st(guiding_center_4d_loop_ode_init, f_loop, ϑ, Δt, 4, nloop, ntime, nsave, DT)
end

function guiding_center_4d_iode_poincare_invariant_1st(Δt, nloop, ntime, nsave, DT=Float64)
   PoincareInvariant1st(guiding_center_4d_loop_iode_init, f_loop, ϑ, Δt, 4, nloop, ntime, nsave, DT)
end

function guiding_center_4d_vode_poincare_invariant_1st(Δt, nloop, ntime, nsave, DT=Float64)
   PoincareInvariant1st(guiding_center_4d_loop_vode_init, f_loop, ϑ, Δt, 4, nloop, ntime, nsave, DT)
end
