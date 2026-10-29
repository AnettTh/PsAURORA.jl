using AURORA
using AURORA; mₑ, qₑ, ε₀
using LinearAlgebra

# IDEA: Figure out how to define λ_grid from the plasma state itself? Maybe use find_R_max to locate where the fieldline hits the ionosphere?
"""
    PlasmaState

Cold plasma parameters along a magnetic field line, as a function of magnetic latitude.

# Fields

- `λ`: Magnetic latitude grid [rad].
- `n_e`: Cold electron number density grid [m⁻³].
- `Ω_e`: Electron gyrofrequency grid [rad/s].
- `ω_pe`: Electron plasma frequency grid [rad/s].
- `ω_lb`: Lower-band chorus wave frequency range [rad/s].
- `ϕ`: Longitude [rad].
- `ne_model`: Model used to make density distribution.
"""
struct PlasmaState{M}
    λ        :: Vector{Float64}
    n_e      :: Vector{Float64}
    Ω_e      :: Vector{Float64}
    ω_pe     :: Vector{Float64}
    ω_lb     :: Vector{Float64}
    ϕ        :: Float64
    ne_model :: M
end

# TODO: See if these functions really are needed anymore
"""
    Ωe_at_λ(λ, L, ϕ, magnetic_field::AbstractMagneticField)

Calculate the electron gyrofreequency `Ω_e` as a function of magnetic latitude `λ`.

# Arguments

- `λ`: Magnetic latitude [rad].
- `L`: L-shell position [RE].
- `ϕ`: Longitude [rad].
- `magnetic_field`: Function that describes the magnetic field based on position.

# Returns

- The gyrofrequency as a fuction of magnetic latitude [rad/s].
"""
function Ωe_at_λ(λ, L, ϕ, magnetic_field::AbstractMagneticField)
    B = magnetic_field(Spherical(L, λ, ϕ))
    return abs(qₑ) * norm(B) / mₑ
end


"""
    ωpe_at_λ(ne_model::Function, L, λ)

Calculate the plasma frequency of electrons for some latitude `λ` and L-shell `L`, given
some density model.

# Arguments

- `ne_model`: Function describing the density model in use, dependent on `L` and `λ`.
- `L`: Position L-shell.
- `λ`: Position latitude [rad].

# Returns

- Electron plasma frequency [rad/s].
"""
function ωpe_at_λ(ne_model::AbstractDensityModel, L, λ)

    ne = ne_model(L, λ)

    return @. sqrt(ne * qₑ^2 / (mₑ * ε₀))
end


"""
    PlasmaState(
    λ_grid::AbstractVector,
    ne_model::AbstractDensityModel,
    magnetic_field::Function,
    L;
    lb_low=0.1,
    lb_high=0.5
)
    PlasmaState(λ_grid, ϕ_eq, magnetic_field::AbstractMagneticField, L; kwargs...)

Construct a `PlasmaState` along at a magnetic field line on a given L-shell.

Computes the electron gyrofrequency `Ω_e` from the magnetic field model, the plasma
frequency `ω_pe` from the electron density `n_e0` and the lower-band chorus frequency range
`ω_lb` as `lb_low*Ω_e` to `lb_high*Ω_e` at the equator. The default range is from 0.25-0.5
of Ω_e.

# Arguments

- `λ_grid`: Magnetic latitude grid [rad].
- `ϕ`: Longitude [rad].
- `ne_model`: Function for the position-dependent density model.
- `magnetic_field`: Magnetic field function `f(L, λ)`.
- `L`: L-shell number [RE].

# Keyword Arguments

- `lb_low`: Fraction of `Ω_e` for the lower boundary of the lower-band chorus wave, default
  is 0.1.
- `lb_high`: Fraction of `Ω_e` for the upper boundary of the lower-band chorus wave, default
  is 0.5.

# Returns

- A `PlasmaState` with all quantities evaluated on `λ_grid`.
"""
function PlasmaState(
    λ_grid::AbstractVector,
    ϕ::Float64,
    ne_model::AbstractDensityModel,
    magnetic_field::AbstractMagneticField,
    L;
    lb_low=0.1,
    lb_high=0.5
)

    n_λ = length(λ_grid)

    # Electron plasma frequency
    n_e = ne_model.(L, λ_grid)
    ω_pe = @. sqrt(n_e * qₑ^2 / (mₑ * ε₀))

    # Electron cyclotron frequency
    Ω_e = Ωe_at_λ.(λ_grid, L, ϕ, Ref(magnetic_field))

    # Whistler-mode chorus wave frequencies
    Ω_eq = minimum(Ω_e)
    ω_lb = range(Ω_eq*lb_low, Ω_eq*lb_high, length=n_λ)

    return PlasmaState(collect(λ_grid), n_e, Ω_e, ω_pe, collect(ω_lb), ϕ, ne_model)
end

# TODO: Try to remove all calls like this in the script, using only the one with a λ_grid and density model
function PlasmaState(λ_grid, ϕ_eq, magnetic_field::AbstractMagneticField, L; kwargs...)
    ne_model = DentonDensity(magnetic_field, L, ϕ_eq)
    return PlasmaState(λ_grid, ϕ_eq, ne_model, magnetic_field, Float64(L); kwargs...)
end

function PlasmaState(λ_grid, ne_model::DentonDensity; kwargs...)
    return PlasmaState(
        λ_grid,
        ne_model.ϕ_eq,
        ne_model,
        ne_model.magnetic_field,
        Float64(ne_model.L);
        kwargs...
    )
end

function Base.show(io::IO, plasma::PlasmaState)
    # Plasma frequency corresponding to λ of the lowest Ω_e (not stored in the struct)
    idx_eq = argmin(plasma.Ω_e)
    Ω_e_eq = plasma.Ω_e[idx_eq]
    ω_pe_eq = plasma.ω_pe[idx_eq]

    println(io, "PlasmaState:")
    println(io, "  λ range   = $(round(rad2deg(plasma.λ[1]), digits=1))° to " *
                "$(round(rad2deg(plasma.λ[end]), digits=1))° ($(length(plasma.λ)) points)")
    println(io, "  ϕ         = $(round(rad2deg(plasma.ϕ), digits=1))° " *
                "($(round(longitude_to_MLT(plasma.ϕ, degrees=false), digits=1)) MLT)")
    println(io, "  n_e range = $(round(minimum(plasma.n_e)/1e6, digits=2)) - " *
                "$(round(maximum(plasma.n_e)/1e6, digits=2)) cm⁻³")
    println(io, "  Ω_e(eq)   = $(round(Ω_e_eq, digits=1)) rad/s")
    println(io, "  ω_pe(eq)  = $(round(ω_pe_eq, digits=1)) rad/s")
    println(io, "  ω_pe/Ω_e  = $(round(ω_pe_eq/Ω_e_eq, digits=2))")
    println(io, "  ω_lb      = $(round(plasma.ω_lb[1], digits=1)) - " *
                "$(round(plasma.ω_lb[end], digits=1)) rad/s")
    print(io,   "  ne_model  = $(typeof(plasma.ne_model))")
end
