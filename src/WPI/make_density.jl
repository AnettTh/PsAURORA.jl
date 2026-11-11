using AURORA
using LinearAlgebra

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
    L::Float64                  # Original equatorial L, used for initial position. #TODO: Should be R_eq
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
#function DentonDensity(
#    magnetic_field::AbstractMagneticField,
#    L,
#    ϕ0;
#    duct_enhancement=0.0,
#)
#
#    # Based on R_max regardless of the magnetic field model used
#    R_max, trace = find_R_max(magnetic_field, Dipolar(L, 0.0, ϕ0); store_trace=true)
#
#    # Find ϕ_eq
#    B = [norm(magnetic_field(Cartesian(x, y, z)))
#            for (x, y, z) in zip(trace.xs, trace.ys, trace.zs)]
#    i = argmin(B)
#    ϕ_eq = mod2pi(atan(trace.ys[i], trace.xs[i]))
#
#    # If above threshold for model, use the upper limit but with extended distance,
#    # still throws below
#    L_bin = min(L, 8.2)
#
#    return DentonDensity(R_max, ϕ_eq, duct_enhancement, L, L_bin, magnetic_field)
#
#end
# TODO: Update docstring
function DentonDensity(
    magnetic_field::AbstractMagneticField,
    position_eq::Cartesian;
    duct_enhancement=0.0
)
    # Find maximal distance from Earth
    R_max = find_R_max(magnetic_field, position_eq)

    # Find positional variables needed for the Denton density model
    ϕ_eq = mod2pi(atan(position_eq.y, position_eq.x))
    R_eq = sqrt(position_eq.x^2 + position_eq.y^2 + position_eq.z^2)

    # If the max distance is above the treshold, use the threshold as the bin. Still throws
    # for lower limit
    L_bin = min(R_eq/RE, 8.2)
    return DentonDensity(R_max, ϕ_eq, duct_enhancement, R_eq/RE, L_bin, magnetic_field)
end


function DentonDensity(
    magnetic_field::DipoleMagneticField,
    L;
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

function electron_density(model::DentonDensity, L, λ)
    # TODO: Check that model.L makes sense here
    n_e = denton_density_model(Dipolar(model.L, λ, model.ϕ_eq), model.R_max, model.L_Denton_bin)
    return n_e * (1.0 + model.duct_enhancement)
end

# Make both callable as functions
(model::AbstractDensityModel)(L, λ) = electron_density(model, L, λ)
