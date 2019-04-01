
using PoincareInvariants

export guiding_center_4d_loop_ode,
       guiding_center_4d_loop_iode

export guiding_center_4d_ode_poincare_invariant_1st,
       guiding_center_4d_iode_poincare_invariant_1st,
       guiding_center_4d_vode_poincare_invariant_1st


function f_loop(i, n)
   f_loop(i/n)
end

function get_initial_conditions(n)
   q₀ = zeros(4, n)

   for i in 1:n
       q₀[:,i] = f_loop(i, n)
   end

   return q₀
end


function guiding_center_4d_loop_ode(n)
   q₀ = get_initial_conditions(n)
   guiding_center_4d_ode(q₀; periodic=false)
end

function guiding_center_4d_loop_iode(n)
   q₀ = get_initial_conditions(n)
   guiding_center_4d_iode(q₀; periodic=false)
end


function guiding_center_4d_ode_poincare_invariant_1st(Δt, nloop, ntime, nsave, DT=Float64)
   guiding_center_4d_ode_init(q₀) = guiding_center_4d_ode(q₀; periodic=false)
   PoincareInvariant1st(guiding_center_4d_ode_init, f_loop, α, Δt, 4, nloop, ntime, nsave, DT)
end

function guiding_center_4d_iode_poincare_invariant_1st(Δt, nloop, ntime, nsave, DT=Float64)
   guiding_center_4d_iode_init(q₀) = guiding_center_4d_iode(q₀; periodic=false)
   PoincareInvariant1st(guiding_center_4d_iode_init, f_loop, α, Δt, 4, nloop, ntime, nsave, DT)
end

function guiding_center_4d_vode_poincare_invariant_1st(Δt, nloop, ntime, nsave, DT=Float64)
   guiding_center_4d_vode_init(q₀) = guiding_center_4d_vode_formal_lagrangian(q₀; periodic=false)
   PoincareInvariant1st(guiding_center_4d_vode_init, f_loop, α, Δt, 4, nloop, ntime, nsave, DT)
end
