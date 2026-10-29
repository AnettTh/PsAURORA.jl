using AURORA
using Geodesy
using PythonCall
using Dates

## Coordinate types
struct Cartesian
    x::Float64
    y::Float64
    z::Float64
end

# NOTE: This might be an issue, as it is spherical dipole coordinates, not pure spherical coordinates!
struct Spherical
    L::Float64
    λ::Float64
    ϕ::Float64
end

struct EastNorthUp
    e::Float64
    n::Float64
    u::Float64
end

struct LTPReference
    lat::Float64
    lon::Float64
    alt::Float64
end

struct MagneticCoordinate
    mlat_deg::Float64   # Magnetic apex latitude [°]
    mlon_deg::Float64   # Magnetic apex longitude [°]
    alt::Float64
    mlt::Float64
end


# Adding some locations that can be called by name
LTPReference(name::Symbol) = LTPReference(
    Dict(
        :tromso       => (69.66,  18.94, 0.0),
        :longyearbyen => (78.22,  15.63, 0.0),
        :sodankyla    => (67.37,  26.63, 0.0),
    )[name]...
)

# Define conversion between the types
function Spherical(r::Cartesian)
    r_mag = sqrt(r.x^2 + r.y^2 + r.z^2)
    λ     = asin(r.z / r_mag)
    ϕ     = atan(r.y, r.x)
    L     = r_mag / (RE * cos(λ)^2)
    return Spherical(L, λ, ϕ)
end

function Cartesian(r::Spherical)
    r_mag = r.L * RE * cos(r.λ)^2
    x     = r_mag * cos(r.λ) * cos(r.ϕ)
    y     = r_mag * cos(r.λ) * sin(r.ϕ)
    z     = r_mag * sin(r.λ)
    return Cartesian(x, y, z)
end


function EastNorthUp(r::Cartesian, ref::LTPReference)
    # Make ECEF Geodesy-element of position
    ecef = ECEF(r.x, r.y, r.z)

    # Make ECEF Geodesy-element of reference
    ref_lla = LLA(ref.lat, ref.lon, ref.alt)
    ref_ecef = ECEF(ref_lla, wgs84)

    # Convert position to ENU with reference
    enu = ENU(ecef, ref_ecef, wgs84)

    return EastNorthUp(enu.e, enu.n, enu.u)
end

function EastNorthUp(r::Cartesian, ref::Symbol)
    return EastNorthUp(r, LTPReference(ref))
end

function EastNorthUp(r::Spherical, ref::LTPReference)
    return EastNorthUp(Cartesian(r), ref)
end

function EastNorthUp(r::Spherical, ref::Symbol)
    return EastNorthUp(Cartesian(r), LTPReference(ref))
end

function Geograpic2Magnetic(ref::LTPReference; time::DateTime=DateTime(2020, 1, 1, 2, 0, 0))
    apex = pyimport("apexpy").Apex(date=year(time))
    mlat, mlon = apex.geo2apex(ref.lat, ref.lon, ref.alt / 1e3)
    mlt = apex.mlon2mlt(pyconvert(Float64, mlon), time)
    return MagneticCoordinate(
        pyconvert(Float64, mlat),
        pyconvert(Float64, mlon),
        ref.alt,
        pyconvert(Float64, mlt),
        )
end


"""
    longitude_to_MLT(ϕ; degrees::Bool=true)

Convert from longitude to Magnetic Local Time (MLT).

Takes the longitude in either degrees (default) or radians, and converts the position to
magnetic local time.

# Arguments

- `ϕ`: Longitude [° or rad].

# Keyword Arguments

- `degrees`: Indicates the units of the input longitude, default is degrees.
"""
function longitude_to_MLT(ϕ; degrees::Bool=true)
    if degrees
        return  mod(12.0 + (ϕ / 15.0), 24.0)
    else
        return mod(12.0 + (rad2deg(ϕ) / 15.0), 24.0)
    end
end
