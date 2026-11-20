using AURORA
using LinearAlgebra

# Density models for plasma electron number density along field lines.
# Implements AbstractDensityModel interface with ConstantDensity and DentonDensity.

"""
    AbstractDensityModel

Abstract type for electron density models along magnetic field lines.
"""
abstract type AbstractDensityModel end

# TODO: Make this into the ∝ B²(?)-model, for comparison?
"""
    ConstantDensity(n_e)

Constant electron density model.
"""
struct ConstantDensity <: AbstractDensityModel
    n_e::Float64    # Electron number density [m⁻³]
end

"""
    DentonDensity

Density model along a magnetic field line, as described by Denton et al. (2002).
"""
struct DentonDensity{F<:AbstractMagneticField} <: AbstractDensityModel
    R_max::Float64              # Maximum radial distance along field line [m].
    ϕ_eq::Float64               # Equatorial longitude [rad].
    duct_enhancement::Float64   # Fraction of density enhancement, 0.0 = not ducted.
    R_eq::Float64               # Position at equator [RE].
    L_Denton_bin::Float64       # L for the bin used when constructing the density model.
    magnetic_field::F           # Magnetic field model.
end


"""
    DentonDensity(magnetic_field, L, ϕ_eq)
    DentonDensity(::typeof(dipole_field), L, ϕ_eq)

Constructor for running the Denton density model with either arbitrary magnetic field, or
using an analytical solution for a dipole field.

# Arguments

- `magnetic field`: Arbitrary magnetic field function that can take either
  `Dipolar(L, λ, ϕ)` or `Cartesian(x, y, z)`.
- `L`: L-shell of initial position.
- `ϕ_eq`: Longitude of position at equator [rad].

# Returns

- A `struct` holding the necessary information to run the Denton density model.

"""
function DentonDensity(fieldline::AbstractFieldLine; duct_enhancement::Float64=0.0)
    # If the max distance is above the treshold, use the threshold as the bin. Still throws
    # for lower limit
    R_eq = fieldline.r_eq / RE
    L_bin = min(R_eq, 8.2)

    return DentonDensity(
        fieldline.r_max,
        fieldline.ϕ_eq,
        duct_enhancement,
        R_eq,
        L_bin,
        fieldline.magnetic_field
        )
end


function DentonDensity(
    magnetic_field::DipoleMagneticField,
    L::Real;
    duct_enhancement=0.0,
)
    R_max = Float64(L) * RE
    return DentonDensity(R_max, 0.0, duct_enhancement, L, L, magnetic_field)
end


"""
    electron_density(model::AbstractDensityModel, L, λ)

Compute electron number density at a given L-shell and latitude, using any defined density
model.

# Arguments

- `model`: Density model used for the computation.
- `L`: L-shell of current position.
- `λ`: Magnetic latitude of current position [rad].

# Returns

- Electron number density [m⁻³].
"""
function electron_density(model::ConstantDensity, L, λ)
    return model.n_e
end

function electron_density(model::DentonDensity, position::Cartesian)
    n_e = denton_density_model(position, model.R_max, model.L_Denton_bin)
    return n_e * (1.0 + model.duct_enhancement)
end

function electron_density(model::ConstantDensity, position::Cartesian)
    return model.n_e
end

# Make both callable as functions
(model::AbstractDensityModel)(L, λ) = electron_density(model, L, λ)
(model::AbstractDensityModel)(position::Cartesian) = electron_density(model, position)

function Base.show(io::IO, d::DentonDensity)
    println(io, "DentonDensity:")
    println(io, "  L-bin:            $(d.L_Denton_bin)")
    println(io, "  Duct enhancement: $(d.duct_enhancement)")
end
