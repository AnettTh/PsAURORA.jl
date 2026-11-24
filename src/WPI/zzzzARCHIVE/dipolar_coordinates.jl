"""
    Dipolar

Dipolar coordinates, given in metres and radians. The use of these coordinates assumes a
dipolar geometry, and hence must not be used together with the Tsyganenko magnetic field
model naively. Is constructed either by giving it the fields, or by a `Cartesian` type.

# Fields

- `L`: L-shell of the position, i.e. the point of the fieldline furthest away from the
  center of the Earth, in the equatorial plane, given in units of Earth radii [RE].
- `λ`: Latitude of the position in radians, measured from the equatorial plane and north
  (south) being positive (negative).
- `ϕ`: Longitude of the position in radians, to complete the system. Only relevant for
  book-keeping purposes, as the dipolar approximation has axisymmetry.

  # Constructors

    Dipolar(L::Float64, λ::Float64, ϕ::Float64)
    Dipolar(r::Cartesian)

  """
struct Dipolar
    L::Float64
    λ::Float64
    ϕ::Float64
end

function Cartesian(r::Dipolar)
    cλ = cos(r.λ)

    r_mag = r.L * RE * cλ^2

    x = r_mag * cλ * cos(r.ϕ)
    y = r_mag * cλ * sin(r.ϕ)
    z = r_mag * sin(r.λ)
    return Cartesian(x, y, z)
end

function Dipolar(r::Cartesian)
    r_mag = sqrt(r.x^2 + r.y^2 + r.z^2)
    λ = asin(r.z / r_mag)
    ϕ = mod2pi(atan(r.y, r.x))
    L = r_mag / (RE * cos(λ)^2)
    return Dipolar(L, λ, ϕ)
end

function SM_to_ENU(r_sm::Dipolar, ref::LTPReference, t::DateTime)
    return SM_to_ENU(Cartesian(r_sm), ref, t)
end

function SM_to_ENU(r_sm::Dipolar, ref::Symbol, t::DateTime)
    return SM_to_ENU(Cartesian(r_sm), LTPReference(ref), t)
end

function dipole_field(pos::Dipolar)
    r = Cartesian(pos)
    return dipole_field(r.x, r.y, r.z)
end

(model::DipoleMagneticField)(r::Dipolar) = dipole_field(r)
