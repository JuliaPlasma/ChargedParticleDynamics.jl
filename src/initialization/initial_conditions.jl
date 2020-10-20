
using LinearAlgebra

using ElectromagneticFields: crossproduct

export InitialConditions, InitialConditionsGC
export charged_particle, guiding_center, pauli_particle
export e, me, mp, md, mα

const e  = 1.602176634E-19    # electron charge
const me = 9.1093837015E-31   # electron mass
const mp = 1.6726219237E-27   # proton mass
const md = 3.3435837724E-27   # deuteron mass
const mα = 6.6446573357E-27   # α-particle mass


"""
Store initial conditions for charged particle, Pauli particle and guiding center models.

Fields:
* `x`: particle position
* `X`: gyro center position
* `ρ`: gyro radius vector
* `vvec`: velocity vector
* `vpar`: parallel velocity vector
* `vper`: perpendicular velocity vector
* `v`: absolute value of velocity
* `u`: absolute value of parallel velocity
* `μ`: magnetic moment
* `θ`: gyro angle ∈ [0,2π]
* `α`: pitch angle ∈ [0,π/2]
* `ω`: gyro frequency
* `mass`: mass in kg
* `energy`: energy in eV
* `charge`: charge number
"""
struct InitialConditions{T <: Real}
    x::Vector{T}
    X::Vector{T}
    ρ::Vector{T}

    vvec::Vector{T}
    vpar::Vector{T}
    vper::Vector{T}

    v::T
    u::T
    μ::T

    θ::T
    α::T
    ω::T

    mass::T
    energy::T
    charge::Int
end


"""
Compute initial conditions from the following arguments:
* `X`: gyro center position
* `θ`: gyro angle
* `α`: pitch angle
* `E`: energy
* `M`: mass
* `C`: charge number
* `a`, b, c: magnetic field unit vectors in physical coordinates
* `b`: magnetic field unit vector in contravariant coordinates
* `B`: amplitude of magnetic field
* `g̅`: inverse metric coefficients
* `DF̄`: inverse Jacobian matrix
* `J`: Jacobian determinant
* `l=1`: length normalization
"""
function InitialConditions(X::AbstractVector{T}, θ::T, α::T, Etot::T, M::T, C, â::Function, b̂::Function, ĉ::Function, b::Function, B::Function, g̅::Function, DF̄::Function, J::Function; l₀=1) where {T}
    # computet parallel and perpedicular energy according to pitch angle
    Epar = Etot * (1 - sin(α))
    Eper = Etot * sin(α)

    # compute absolute value of velocity v from energy Etot,
    # parallel velocity u, perpendicular velocity w, and
    # magnetic moment μ
    v = sqrt(2 * Etot * M / e) / l₀
    u = sqrt(2 * Epar * M / e) / l₀
    w = sqrt(2 * Eper * M / e) / l₀
    μ = w^2 / 2B(0,X)

    # compute v̂par and v̂per in physical coordinates
    v̂par = u * b̂(0,X)
    v̂per = w * ( - â(0,X) * sin(θ) - ĉ(0,X) * cos(θ) )

    # transform vpar and vper to contravariant coordinates
    vpar = DF̄(0,X) * v̂par
    vper = DF̄(0,X) * v̂per

    # compute total velocity vector field
    vvec = vpar .+ vper

    # compute gyro radius vector ρ and particle position x
    ω = e * B(0,X) / M
    ρ = [crossproduct(b(0,X), vvec, g̅(0,X), J(0,X), l) for l in 1:3] / B(0,X)
    x = X .+ ρ

    InitialConditions{T}(x, X, ρ, vvec, vpar, vper, v, u, μ, θ, α, ω, M, Etot, C)
end


"""
Compute initial conditions from the following arguments:
* `X`: gyro center position
* `θ`: gyro angle
* `u`: absolute value of parallel velocity
* `μ`: magnetic moment
* `M`: mass
* `C`: charge number
* `a`, b, c: magnetic field unit vectors in physical coordinates
* `b`: magnetic field unit vector in contravariant coordinates
* `B`: amplitude of magnetic field
* `g̅`: inverse metric coefficients
* `DF̄`: inverse Jacobian matrix
* `J`: Jacobian determinant
* `l=1`: length normalization
"""
function InitialConditionsGC(X::AbstractVector{T}, θ::T, u::T, μ::T, M::T, C, â::Function, b̂::Function, ĉ::Function, b::Function, B::Function, g̅::Function, DF̄::Function, J::Function; l₀=1) where {T}
    # compute absolute value of perpendicular and total velocity
    w = sqrt(2B(0,X) * μ)
    v = sqrt(u^2 + w^2)

    # compute perpendicular and total energy
    Eper = e * l₀^2 * w^2 / 2 / M 
    Etot = e * l₀^2 * v^2 / 2 / M

    # compute pitch angle
    α = asin(Eper / Etot)

    # compute v̂par and v̂per in physical coordinates
    v̂par = u * b̂(0,X)
    v̂per = w * ( - â(0,X) * sin(θ) - ĉ(0,X) * cos(θ) )

    # transform vpar and vper to contravariant coordinates
    vpar = DF̄(0,X) * v̂par
    vper = DF̄(0,X) * v̂per

    # compute total velocity vector field
    vvec = vpar .+ vper

    # compute gyro radius vector ρ and particle position x
    ω = e * B(0,X) / M
    ρ = [crossproduct(b(0,X), vvec, g̅(0,X), J(0,X), l) for l in 1:3] / B(0,X)
    x = X .+ ρ

    InitialConditions{T}(x, X, ρ, vvec, vpar, vper, v, u, μ, θ, α, ω, M, Etot, C)
end


function Base.show(io::IO, ics::InitialConditions)
    print(io, "Charged Particle Initial Conditions with\n")
    print(io, "  x   = ", ics.x, "\n")
    print(io, "  X   = ", ics.X, "\n")
    print(io, "  ρ   = ", ics.ρ, "\n")
    print(io, "  v   = ", ics.vvec, "\n")
    print(io, "  v∥  = ", ics.vpar, "\n")
    print(io, "  v⟂  = ", ics.vper, "\n")
    print(io, "  |v| = ", ics.v, "\n")
    print(io, "  u   = ", ics.u, "\n")
    print(io, "  μ   = ", ics.μ, "\n")
    print(io, "  θ   = ", ics.θ, "\n")
    print(io, "  α   = ", ics.α, "\n")
    print(io, "  ω   = ", ics.ω, "\n")
    print(io, "  M   = ", ics.mass, "\n")
    print(io, "  E   = ", ics.energy, "\n")
    print(io, "  C   = ", ics.charge)
end


"""
Extracts the charged particle initial conditions and returns the tuple `(x,v)`.
If the keyword argument `noncanonical` is set to `true`, the functions returns
the vector `vcat(x,v)`.
"""
function charged_particle(ics::InitialConditions{T}; noncanonical=false) where {T}
    if noncanonical
        return vcat(ics.x, ics.v)
    else
        return (ics.x, ics.vvec)
    end
end


"""
Extracts the guiding center initial conditions and returns the tuple `(vcat(X,u),μ)`.
"""
function guiding_center(ics::InitialConditions{T}) where {T}
    (vcat(ics.X, ics.u), ics.μ)
end


"""
Extracts the Pauli particle initial conditions and returns the tuple `(X,vpar,μ)`.
If the keyword argument `noncanonical` is set to `true`, the functions returns
the tuple `(vcat(X,vpar),μ)`.
"""
function pauli_particle(ics::InitialConditions{T}; noncanonical=false) where {T}
    if noncanonical
        return (vcat(ics.X, ics.vpar), ics.μ)
    else
        return (ics.X, ics.vpar, ics.μ)
    end
end
