
using DecFP
# using DoubleDouble


function α(t, q, p)
    p[1] = α1(t,q)
    p[2] = α2(t,q)
    p[3] = α3(t,q)
    p[4] = α4(t,q)
    nothing
end


function ω(t, q, Ω)
    Ω[1,1] = 0
    Ω[1,2] = dα1d2(t,q) - dα2d1(t,q)
    Ω[1,3] = dα1d3(t,q) - dα3d1(t,q)
    Ω[1,4] = dα1d4(t,q)

    Ω[2,1] = dα2d1(t,q) - dα1d2(t,q)
    Ω[2,2] = 0
    Ω[2,3] = dα2d3(t,q) - dα3d2(t,q)
    Ω[2,4] = dα2d4(t,q)

    Ω[3,1] = dα3d1(t,q) - dα1d3(t,q)
    Ω[3,2] = dα3d2(t,q) - dα2d3(t,q)
    Ω[3,3] = 0
    Ω[3,4] = dα3d4(t,q)

    Ω[4,1] = dα4d1(t,q) - dα1d4(t,q)
    Ω[4,2] = dα4d2(t,q) - dα2d4(t,q)
    Ω[4,3] = dα4d3(t,q) - dα3d4(t,q)
    Ω[4,4] = 0

    nothing
end


function D²ϑ₁(t, q, D²ϑ)
    D²ϑ[1,1] = d²A₁dx₁dx₁(t,q) + u(t,q) * d²b₁dx₁dx₁(t,q)
    D²ϑ[1,2] = d²A₁dx₁dx₂(t,q) + u(t,q) * d²b₁dx₁dx₂(t,q)
    D²ϑ[1,3] = d²A₁dx₁dx₃(t,q) + u(t,q) * d²b₁dx₁dx₃(t,q)
    D²ϑ[1,4] = db₁dx₁(t,q)

    D²ϑ[2,1] = d²A₁dx₂dx₁(t,q) + u(t,q) * d²b₁dx₂dx₁(t,q)
    D²ϑ[2,2] = d²A₁dx₂dx₂(t,q) + u(t,q) * d²b₁dx₂dx₂(t,q)
    D²ϑ[2,3] = d²A₁dx₂dx₃(t,q) + u(t,q) * d²b₁dx₂dx₃(t,q)
    D²ϑ[2,4] = db₁dx₂(t,q)

    D²ϑ[2,1] = d²A₁dx₃dx₁(t,q) + u(t,q) * d²b₁dx₃dx₁(t,q)
    D²ϑ[2,2] = d²A₁dx₃dx₂(t,q) + u(t,q) * d²b₁dx₃dx₂(t,q)
    D²ϑ[2,3] = d²A₁dx₃dx₃(t,q) + u(t,q) * d²b₁dx₃dx₃(t,q)
    D²ϑ[2,4] = db₁dx₃(t,q)

    D²ϑ[4,1] = db₁dx₁(t,q)
    D²ϑ[4,2] = db₁dx₂(t,q)
    D²ϑ[4,3] = db₁dx₃(t,q)
    D²ϑ[4,4] = 0

    nothing
end

function D²ϑ₂(t, q, D²ϑ)
    D²ϑ[1,1] = d²A₂dx₁dx₁(t,q) + u(t,q) * d²b₂dx₁dx₁(t,q)
    D²ϑ[1,2] = d²A₂dx₁dx₂(t,q) + u(t,q) * d²b₂dx₁dx₂(t,q)
    D²ϑ[1,3] = d²A₂dx₁dx₃(t,q) + u(t,q) * d²b₂dx₁dx₃(t,q)
    D²ϑ[1,4] = db₂dx₁(t,q)

    D²ϑ[2,1] = d²A₂dx₂dx₁(t,q) + u(t,q) * d²b₂dx₂dx₁(t,q)
    D²ϑ[2,2] = d²A₂dx₂dx₂(t,q) + u(t,q) * d²b₂dx₂dx₂(t,q)
    D²ϑ[2,3] = d²A₂dx₂dx₃(t,q) + u(t,q) * d²b₂dx₂dx₃(t,q)
    D²ϑ[2,4] = db₂dx₂(t,q)

    D²ϑ[2,1] = d²A₂dx₃dx₁(t,q) + u(t,q) * d²b₂dx₃dx₁(t,q)
    D²ϑ[2,2] = d²A₂dx₃dx₂(t,q) + u(t,q) * d²b₂dx₃dx₂(t,q)
    D²ϑ[2,3] = d²A₂dx₃dx₃(t,q) + u(t,q) * d²b₂dx₃dx₃(t,q)
    D²ϑ[2,4] = db₂dx₃(t,q)

    D²ϑ[4,1] = db₂dx₁(t,q)
    D²ϑ[4,2] = db₂dx₂(t,q)
    D²ϑ[4,3] = db₂dx₃(t,q)
    D²ϑ[4,4] = 0

    nothing
end

function D²ϑ₃(t, q, D²ϑ)
    D²ϑ[1,1] = d²A₃dx₁dx₁(t,q) + u(t,q) * d²b₃dx₁dx₁(t,q)
    D²ϑ[1,2] = d²A₃dx₁dx₂(t,q) + u(t,q) * d²b₃dx₁dx₂(t,q)
    D²ϑ[1,3] = d²A₃dx₁dx₃(t,q) + u(t,q) * d²b₃dx₁dx₃(t,q)
    D²ϑ[1,4] = db₃dx₁(t,q)

    D²ϑ[2,1] = d²A₃dx₂dx₁(t,q) + u(t,q) * d²b₃dx₂dx₁(t,q)
    D²ϑ[2,2] = d²A₃dx₂dx₂(t,q) + u(t,q) * d²b₃dx₂dx₂(t,q)
    D²ϑ[2,3] = d²A₃dx₂dx₃(t,q) + u(t,q) * d²b₃dx₂dx₃(t,q)
    D²ϑ[2,4] = db₃dx₂(t,q)

    D²ϑ[2,1] = d²A₃dx₃dx₁(t,q) + u(t,q) * d²b₃dx₃dx₁(t,q)
    D²ϑ[2,2] = d²A₃dx₃dx₂(t,q) + u(t,q) * d²b₃dx₃dx₂(t,q)
    D²ϑ[2,3] = d²A₃dx₃dx₃(t,q) + u(t,q) * d²b₃dx₃dx₃(t,q)
    D²ϑ[2,4] = db₃dx₃(t,q)

    D²ϑ[4,1] = db₃dx₁(t,q)
    D²ϑ[4,2] = db₃dx₂(t,q)
    D²ϑ[4,3] = db₃dx₃(t,q)
    D²ϑ[4,4] = 0

    nothing
end

function D²ϑ₄(t, q, D²ϑ)
    D²ϑ .= 0
    nothing
end


function dα(t, q, dα)
    dα[1,1] = dα1d1(t,q)
    dα[1,2] = dα1d2(t,q)
    dα[1,3] = dα1d3(t,q)
    dα[1,4] = dα1d4(t,q)

    dα[2,1] = dα2d1(t,q)
    dα[2,2] = dα2d2(t,q)
    dα[2,3] = dα2d3(t,q)
    dα[2,4] = dα2d4(t,q)

    dα[3,1] = dα3d1(t,q)
    dα[3,2] = dα3d2(t,q)
    dα[3,3] = dα3d3(t,q)
    dα[3,4] = dα3d4(t,q)

    dα[4,1] = 0
    dα[4,2] = 0
    dα[4,3] = 0
    dα[4,4] = 0

    nothing
end


function β1(t,q)
   return dα3d2(t,q) - dα2d3(t,q)
end

function β2(t,q)
   return dα1d3(t,q) - dα3d1(t,q)
end

function β3(t,q)
   return dα2d1(t,q) - dα1d2(t,q)
end


function β(t,q)
   return sqrt(β1(t,q)^2 + β2(t,q)^2 + β3(t,q)^2)
end


function toroidal_momentum(t,q)
    α3(t,q)
end


function hamiltonian(t,q)
    0.5 * u(t,q)^2 + μ*B(t,q)
end


function dHd1(t, q)
    μ * dBdx₁(t,q)
end

function dHd2(t, q)
    μ * dBdx₂(t,q)
end

function dHd3(t, q)
    μ * dBdx₃(t,q)
end

function dHd4(t, q)
    u(t,q)
end


function dH(t, q, dH)
    dH[1] = dHd1(t, q)
    dH[2] = dHd2(t, q)
    dH[3] = dHd3(t, q)
    dH[4] = dHd4(t, q)
    nothing
end


function f1(t, q, v)
    dα1d1(t,q) * v[1] + dα2d1(t,q) * v[2] + dα3d1(t,q) * v[3]
end

function f2(t, q, v)
    dα1d2(t,q) * v[1] + dα2d2(t,q) * v[2] + dα3d2(t,q) * v[3]
end

function f3(t, q, v)
    dα1d3(t,q) * v[1] + dα2d3(t,q) * v[2] + dα3d3(t,q) * v[3]
end

function f4(t, q, v)
    dα1d4(t,q) * v[1] + dα2d4(t,q) * v[2] + dα3d4(t,q) * v[3]
end


function guiding_center_4d_periodicity(q₀)
    periodicity = zeros(eltype(q₀), size(q₀,1))
    periodicity[3] = 2π
    return periodicity
end


function guiding_center_4d_ode_v(t, q, v)
    local lB₁ = B₁(t,q)
    local lB₂ = B₂(t,q)
    local lB₃ = B₃(t,q)

    local lβ₁ = β1(t,q)
    local lβ₂ = β2(t,q)
    local lβ₃ = β3(t,q)
    local lβ  = lβ₁ * dα1d4(t,q) + lβ₂ * dα2d4(t,q) + lβ₃ * dα3d4(t,q)

    local ∇₁B = dBdx₁(t,q)
    local ∇₂B = dBdx₂(t,q)
    local ∇₃B = dBdx₃(t,q)

    v[1] = + ( u(t,q) * lβ₁ - μ * ( ∇₂B * dα3d4(t,q) - ∇₃B * dα2d4(t,q) ) ) / lβ
    v[2] = + ( u(t,q) * lβ₂ - μ * ( ∇₃B * dα1d4(t,q) - ∇₁B * dα3d4(t,q) ) ) / lβ
    v[3] = + ( u(t,q) * lβ₃ - μ * ( ∇₁B * dα2d4(t,q) - ∇₂B * dα1d4(t,q) ) ) / lβ
    v[4] = - μ * ( ∇₁B * lβ₁
                 + ∇₂B * lβ₂
                 + ∇₃B * lβ₃ ) / lβ

    nothing
end

function guiding_center_4d_ode(q₀=q₀; periodic=true)
    if periodic
        ODE(guiding_center_4d_ode_v, q₀; periodicity=guiding_center_4d_periodicity(q₀))
    else
        ODE(guiding_center_4d_ode_v, q₀)
    end
end


function guiding_center_4d_iode_α(t, q, v, p)
    α(t, q, p)
end

function guiding_center_4d_iode_f(t, q, v, f)
    f[1] = f1(t,q,v) - dHd1(t,q)
    f[2] = f2(t,q,v) - dHd2(t,q)
    f[3] = f3(t,q,v) - dHd3(t,q)
    f[4] = f4(t,q,v) - dHd4(t,q)
    nothing
end

function guiding_center_4d_iode_g(t, q, λ, g)
    g[1] = f1(t,q,λ)
    g[2] = f2(t,q,λ)
    g[3] = f3(t,q,λ)
    g[4] = f4(t,q,λ)
    nothing
end

function guiding_center_4d_iode_v(t, q, p, v)
    guiding_center_4d_ode_v(t, q, v)
end


function guiding_center_4d_p₀(q₀, t₀=0)
    p₀ = zero(q₀)

    if ndims(q₀) == 1
        α(t₀, q₀, p₀)
    else
        for i in 1:size(q₀,2)
            tq = zeros(eltype(q₀), size(q₀,1))
            tp = zeros(eltype(p₀), size(p₀,1))
            simd_copy_xy_first!(tq, q₀, i)
            α(t₀, tq, tp)
            simd_copy_yx_first!(tp, p₀, i)
        end
    end
    p₀
end


function guiding_center_4d_λ(q::Vector{DT}, λ::Vector{DT}, Ω::Matrix{DT}, dh::Vector{DT}, Δt::DT, t::DT) where {DT}
    dH(t, q, dh)
    ω(t, q, Ω)
    simd_mult!(λ, inv(Ω), dh)
    # simd_scale!(λ, 1/Δt)
    nothing
end


function guiding_center_4d_λ₀(q₀::Array{DT}, Δt::DT=1., t₀::DT=0.) where {DT}
    D  = size(q₀, 1)
    λ₀ = zero(q₀)
    q  = zeros(DT, D)
    λ  = zeros(DT, D)
    Ω  = zeros(DT, D, D)
    dh = zeros(DT, D)

    if ndims(q₀) == 1
        guiding_center_4d_λ(q₀, λ₀, Ω, dh, Δt, t₀)
    else
        for i in 1:size(q₀,2)
            simd_copy_xy_first!(q, q₀, i)
            guiding_center_4d_λ(q, λ, Ω, dh, Δt, t₀)
            simd_copy_yx_first!(λ, λ₀, i)
        end
    end
    λ₀
end


# function guiding_center_4d_iode_double(q₀=q₀; periodic=true)
#     guiding_center_4d_iode(Double.(q₀); periodic=periodic)
# end

function guiding_center_4d_iode_dec128(q₀=q₀; periodic=true)
    guiding_center_4d_iode(Dec128.(q₀); periodic=periodic)
end

function guiding_center_4d_iode(q₀=q₀; periodic=true)
    p₀ = guiding_center_4d_p₀(q₀)
    if periodic
        IODE(guiding_center_4d_iode_α, guiding_center_4d_iode_f,
             guiding_center_4d_iode_g, guiding_center_4d_iode_v,
             q₀, p₀; periodicity=guiding_center_4d_periodicity(q₀))
    else
        IODE(guiding_center_4d_iode_α, guiding_center_4d_iode_f,
             guiding_center_4d_iode_g, guiding_center_4d_iode_v,
             q₀, p₀)
    end
end

function guiding_center_4d_iode_λ(q₀=q₀; periodic=true)
    p₀ = guiding_center_4d_p₀(q₀)
    λ₀ = guiding_center_4d_λ₀(q₀)
    if periodic
        IODE(guiding_center_4d_iode_α, guiding_center_4d_iode_f,
             guiding_center_4d_iode_g, guiding_center_4d_iode_v,
             q₀, p₀, λ₀; periodicity=guiding_center_4d_periodicity(q₀))
    else
        IODE(guiding_center_4d_iode_α, guiding_center_4d_iode_f,
             guiding_center_4d_iode_g, guiding_center_4d_iode_v,
             q₀, p₀, λ₀)
    end
end

function guiding_center_4d_vode(q₀=q₀; periodic=true)
    p₀ = guiding_center_4d_p₀(q₀)
    if periodic
        VODE(guiding_center_4d_iode_α, guiding_center_4d_iode_f,
             guiding_center_4d_iode_g, guiding_center_4d_iode_v,
             ω, dH, q₀, p₀; periodicity=guiding_center_4d_periodicity(q₀))
    else
        VODE(guiding_center_4d_iode_α, guiding_center_4d_iode_f,
             guiding_center_4d_iode_g, guiding_center_4d_iode_v,
             ω, dH, q₀, p₀)
    end
end
