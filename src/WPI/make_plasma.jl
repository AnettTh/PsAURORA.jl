using AURORA
using AURORA; mₑ, qₑ, ε₀
using LinearAlgebra

"""
    PlasmaState

Cold plasma parameters along a magnetic field line, as a function of magnetic latitude.

# Fields

- `λ_grid`: Magnetic latitude grid [rad].
- `n_e`: Cold electron number density grid [m⁻³].
- `Ω_e`: Electron gyrofrequency grid [rad/s].
- `ω_pe`: Electron plasma frequency grid [rad/s].
- `ω_lb`: Lower-band chorus wave frequency range [rad/s].
- `ne_model`: Model used to make density distribution.
"""
struct PlasmaState{M}
    λ_grid   :: Vector{Float64} # TODO: Remove from here and use the one in fieldline?
    n_e      :: Vector{Float64}
    Ω_e      :: Vector{Float64}
    ω_pe     :: Vector{Float64}
    ω_lb     :: Vector{Float64}
    ne_model :: M
end

# TODO: Update docstring, major revisions!
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
    fieldline::AbstractFieldLine,
    ne_model::DentonDensity;
    lb_low=0.1,
    lb_high=0.5,
    n_ω=5,
)
    issorted(fieldline.λs) || throw(ArgumentError(
        "Latitude is not monotonic along the line, so it cannot be used as a lookup grid"
    ))

    # Calculate gyrofrequency along the fieldline
    Ω_e = @. abs(qₑ) * fieldline.Bs / mₑ

    # Calculate the density and plasma frequency along the fieldline
    n_e = [
        ne_model(Cartesian(x, y, z))
        for (x, y, z) in zip(fieldline.xs, fieldline.ys, fieldline.zs)
    ]
    ω_pe = @. sqrt(n_e * qₑ^2 / (mₑ * ε₀))

    # Find the gyrofrequency at the equator, used to define the whistler frequency range
    Ω_eq = abs(qₑ) * fieldline.B_eq / mₑ
    ω_lb = collect(range(Ω_eq*lb_low, Ω_eq*lb_high, length=n_ω))

    return PlasmaState(fieldline.λs, n_e, Ω_e, ω_pe, ω_lb, ne_model)
end


function Base.show(io::IO, plasma::PlasmaState)
    # Plasma frequency corresponding to λ of the lowest Ω_e (not stored in the struct)
    idx_eq = argmin(plasma.Ω_e)
    Ω_e_eq = plasma.Ω_e[idx_eq]
    ω_pe_eq = plasma.ω_pe[idx_eq]

    println(io, "PlasmaState:")
    println(io, "  λ range   = $(round(rad2deg(plasma.λ_grid[1]), digits=1))° to " *
                "$(round(rad2deg(plasma.λ_grid[end]), digits=1))° ($(length(plasma.λ_grid)) points)")
    println(io, "  n_e range = $(round(minimum(plasma.n_e)/1e6, digits=2)) - " *
                "$(round(maximum(plasma.n_e)/1e6, digits=2)) cm⁻³")
    println(io, "  Ω_e(eq)   = $(round(Ω_e_eq, digits=1)) rad/s")
    println(io, "  ω_pe(eq)  = $(round(ω_pe_eq, digits=1)) rad/s")
    println(io, "  ω_pe/Ω_e  = $(round(ω_pe_eq/Ω_e_eq, digits=2))")
    println(io, "  ω_lb      = $(round(plasma.ω_lb[1], digits=1)) - " *
                "$(round(plasma.ω_lb[end], digits=1)) rad/s")
    print(io,   "  ne_model  = $(typeof(plasma.ne_model))")
end
