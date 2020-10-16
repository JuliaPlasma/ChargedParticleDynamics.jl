
using LinearAlgebra


const e  = 1.602176634E-19    # electron charge
const me = 9.1093837015E-31   # electron mass
const mp = 1.6726219237E-27   # proton mass
const md = 3.3435837724E-27   # deuteron mass
const mα = 6.6446573357E-27   # α-particle mass


"""
Store initial conditions for charged particle, Pauli particle and guiding center models.

Fields:
* x: particle position
* X: gyro center position
* ρ: gyro radius vector
* vvec: velocity vector
* vpar: parallel velocity vector
* vper: perpendicular velocity vector
* v: absolute value of velocity
* u: absolute value of parallel velocity
* μ: magnetic moment
* θ: gyro angle ∈ [0,2π]
* α: pitch angle ∈ [0,π/2]
* mass: mass in kg
* energy: energy in eV
* charge: charge number
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
* X: gyro center position
* θ: gyro angle
* α: pitch angle
* E: energy
* M: mass
* C: charge number
* a, b, c: magnetic field unit vectors
* B: amplitude of magnetic field
* l=1: length normalization
"""
function InitialConditions(X::AbstractVector{T}, θ::T, α::T, E::T, M::T, C, a::Function, b::Function, c::Function, B::Function; l=1) where {T}
    v = sqrt(2 * E * M / e) / l
    u = v * (1 - sin(α))
    μ = v^2 * sin(α) / 2B(0,X)

    vpar = u * b(0,X)
    vper = v * sin(α) * ( - a(0,X) * sin(θ) - c(0,X) * cos(θ) )
    vvec = vpar .+ vper

    ω = e * B(0,X) / M 
    ρ = b(0,X) × vvec
    x = X .+ ρ

    InitialConditions{T}(x, X, ρ, vvec, vpar, vper, v, u, μ, θ, α, ω, M, E, C)
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
