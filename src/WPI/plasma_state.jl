using AURORA
using AURORA; mₑ, qₑ, ε₀
using LinearAlgebra

# TODO: Figure out how to define λ_grid from the plasma state itself?
"""
    PlasmaState

Cold plasma parameters along a magnetic field line, as a function of magnetic latitude.

# Fields

- `λ`: Magnetic latitude grid [rad].
- `n_e`: Cold electron number density grid [m⁻³].
- `Ω_e`: Electron gyrofrequency grid [rad/s].
- `ω_pe`: Electron plasma frequency grid [rad/s].
- `ω_lb`: Lower-band chorus wave frequency range [rad/s].
"""
struct PlasmaState
    λ    :: Vector{Float64}
    n_e  :: Vector{Float64}
    Ω_e  :: Vector{Float64}
    ω_pe :: Vector{Float64}
    ω_lb :: Vector{Float64}
end


# NOTE: This will not be compatible with Tsyganenko, figure it out!
"""
    Ωe_at_λ(λ, L, magnetic_field)

Calculate the electron gyrofreequency `Ω_e` as a function of magnetic latitude `λ`.

# Arguments

- `λ`: Magnetic latitude [rad].
- `L`: L-shell position [RE].
- `magnetic_field`: Function that describes the magnetic field based on `L` amd `λ`.

# Returns

- The gyrofrequency as a fuction of magnetic latitude [rad/s].
"""
function Ωe_at_λ(λ, L, magnetic_field)

    B = magnetic_field(L, λ)

    return gyro_frequency(norm(B), qₑ, mₑ)
end


"""
    PlasmaState(
    λ_grid::AbstractVector,
    n_e0::Float64,
    magnetic_field::Function,
    L::Float64
)

Construct a `PlasmaState` along at a magnetic field line on a given L-shell.

Computes the electron gyrofrequency `Ω_e` from the magnetic field model, the plasma
frequency `ω_pe` from the electron density `n_e0` and the lower-band chorus frequency range
`ω_lb` as `lb_low*Ω_e` to `lb_high*Ω_e` at the equator. The default range is from 0.25-0.5
of Ω_e.

# Arguments

- `λ_grid`: Magnetic latitude grid [rad].
- `n_e0`: Cold electron number density, assumed constant along the field line [m⁻³].
- `magnetic_field`: Magnetic field function `f(L, λ)`.
- `L`: L-shell number [RE].

# Keyword Arguments

- `lb_low`: Fraction of `Ω_e` for the lower boundary of the lower-band chorus wave.
- `lb_high`: Fraction of `Ω_e` for the upper boundary of the lower-band chorus wave.

# Returns

- A `PlasmaState` with all quantities evaluated on `λ_grid`.
"""
function PlasmaState(
    λ_grid::AbstractVector,
    n_e0::Float64,
    magnetic_field::Function,
    L::Float64;
    lb_low=0.25,
    lb_high=0.5
)
    n_λ = length(λ_grid)

    # Electron plasma frequency
    n_e = fill(n_e0, n_λ)
    ω_pe = @. sqrt(n_e * qₑ^2 / (mₑ * ε₀))

    # Electron cyclotron frequency
    Ω_e = Ωe_at_λ.(λ_grid, L, magnetic_field)

    # Whistler-mode chorus wave frequencies
    Ω_eq = Ω_e[1]
    ω_lb = range(Ω_eq*lb_low, Ω_eq*lb_high, length=n_λ)

    return PlasmaState(collect(λ_grid), n_e, Ω_e, ω_pe, ω_lb)
end
