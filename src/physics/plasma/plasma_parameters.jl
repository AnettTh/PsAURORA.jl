using AURORA
using AURORA; mₑ, qₑ, ε₀
using LinearAlgebra

# TODO: Document, document, document!

struct PlasmaParameters
    λ    :: Vector{Float64}     # Magnetic latitude grid [rad]
    n_e  :: Vector{Float64}     # Cold electron density [m⁻³]
    Ω_e  :: Vector{Float64}     # Gyrofrequency [rad/s]
    ω_pe :: Vector{Float64}     # Plasma frequency [rad/s]
end


# NOTE: This will not be compatible with Tsyganenko
function Ωe_at_λ(λ, L, magnetic_field)

    r = L * RE * cos(λ)^2
    B = magnetic_field(r * cos(λ), 0.0, r * sin(λ))
    Ω_e = gyro_frequency(norm(B), qₑ, mₑ)

    return Ω_e
end

function PlasmaParameters(λ_grid::AbstractVector, n_e0::Float64, magnetic_field::Function, L::Float64)
    n_e = fill(n_e0, length(λ_grid))
    Ω_e = Ωe_at_λ.(λ_grid, L, magnetic_field)
    ω_pe = @. sqrt(n_e * qₑ^2 / (mₑ * ε₀))

    return PlasmaParameters(collect(λ_grid), n_e, Ω_e, ω_pe)
end
