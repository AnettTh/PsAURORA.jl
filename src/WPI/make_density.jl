using AURORA

# Density models for plasma electron number density along field lines.
# Implements AbstractDensityModel interface with ConstantDensity and DentonDensity.

"""
    AbstractDensityModel

Abstract type for electron density models along magnetic field lines.
"""
abstract type AbstractDensityModel end


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
    L::Float64                  # Original equatorial L, used for initial position.
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
function DentonDensity(
    magnetic_field::AbstractMagneticField,
    L,
    ϕ_eq;
    duct_enhancement=0.0,
    hemisphere::Symbol=:north,
    adjusted_model=false,
)
    if !adjusted_model
        # Based on R_max regardless of the magnetic field model used
        R_max = find_R_max(magnetic_field, Dipolar(L, 0.0, ϕ_eq))
        return DentonDensity(R_max, ϕ_eq, duct_enhancement, L, L, magnetic_field)
    else
        # Find the footpoint of the Tsyganenko fieldline
        result = find_footpoint(
            magnetic_field,
            Dipolar(L, 0.0, ϕ_eq);
            hemisphere=hemisphere
        )

        # Use the result of the dipolar L with the same footpoint as the Tsyganenko line to
        # specify which Denton bin to use. This allows for one smooth model, avoiding jumps
        # in density around L∼7
        L_Denton_bin = result.L_dipolar

        return DentonDensity(
            result.L_tsyg_max * RE,
            ϕ_eq,
            duct_enhancement,
            L,
            L_Denton_bin,
            DipoleMagneticField(),
        )
    end
end


function DentonDensity(
    magnetic_field::DipoleMagneticField,
    L,
    ϕ_eq;
    duct_enhancement=0.0,
    kwargs...,      # Accept, but ignore, kwargs about hemispheres and adjusted model
)
    R_max = Float64(L) * RE
    return DentonDensity(R_max, ϕ_eq, duct_enhancement, L, L, magnetic_field)
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

function electron_density(model::DentonDensity, L, λ)
    n_e = denton_density_model(Dipolar(L, λ, model.ϕ_eq), model.R_max, model.L_Denton_bin)
    return n_e * (1.0 + model.duct_enhancement)
end

# Make both callable as functions
(model::AbstractDensityModel)(L, λ) = electron_density(model, L, λ)
